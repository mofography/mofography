#!/usr/bin/env python3
"""Trace env_cluster for any MOF; pass cif slug or path."""
from __future__ import annotations

import copy
import sys
from pathlib import Path

import numpy as np
from ase.io import read

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import mofography as mgr
from mofography.pocket_framework_env import _cosine, fingerprint_cosine_similarity
from mofography.pocket_typing import (
    DEFAULT_MASK_THICKNESS_A,
    DEFAULT_PROBE_RADIUS_PSD_A,
    DEFAULT_SPACING_A,
    ENV_CLUSTER_INITIAL_BAND_GAP_A,
    ENV_CLUSTER_MERGE_COSINE,
    ENV_CLUSTER_MERGE_COSINE_LARGE_GAP,
    ENV_CLUSTER_MERGE_GAP_THRESHOLD_A,
    ENV_CLUSTER_MERGE_MAX_R_GAP_A,
    assign_pocket_types_env_cluster,
    assign_pocket_types_from_segmentation,
    build_pockets_from_segmentation,
    compute_pocket_fingerprints,
    compute_typing_dgrid,
    psd_bins,
)

SEG_KW = dict(
    mask_thickness=DEFAULT_MASK_THICKNESS_A,
    h=0.2,
    min_distance=1.5,
    threshold_abs=1.5,
    extra_padding=5.0,
    return_dgrid=True,
)

CIFS = {
    "cubtc": Path(REPO / "examples/TestCIFs/P1-Cu-BTC.cif"),
    "zif20": Path.home() / "repos/RASPA2/structures/mofs/cif/ZIF-20.cif",
}


def run_seg(cif: Path, dgrid_backend: str = "mic"):
    atoms = read(cif)
    mgr.update_with_covalent_bonds(atoms)
    if dgrid_backend == "rtree":
        dgrid = mgr.dgrid_from_atoms_rtree_rust(atoms, spacing=DEFAULT_SPACING_A).astype(np.float32)
    else:
        dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=DEFAULT_SPACING_A).astype(np.float32)
    results = mgr.do_periodic_segmentation_rust(dgrid, atoms, spacing=DEFAULT_SPACING_A, **SEG_KW)
    if results.get("pocket_psds") is None:
        results["pocket_psds"] = mgr.compute_psds_pocketwise_rust(
            results["dgrid"],
            results.get("regions_pbc", results["regions"]),
            atoms,
            num_samples_per_region=1000,
            bins=psd_bins(results["dgrid"]),
            probe_center_min_radius=DEFAULT_PROBE_RADIUS_PSD_A,
        )
    return atoms, results


def mean_max_r(pids, fps):
    return float(np.mean([fps[p].max_probe_radius_a for p in pids]))


def trace(name: str, cif: Path, dgrid_backend: str = "mic"):
    print(f"\n{'#' * 72}\n# {name}  ({cif.name})\n{'#' * 72}")
    atoms, results = run_seg(cif, dgrid_backend)
    pockets, pbc_groups, _, _ = build_pockets_from_segmentation(results, atoms)
    typing_dgrid = compute_typing_dgrid(atoms, results, "bondi_crc_fallback")
    fps = compute_pocket_fingerprints(
        results, atoms, pockets, pbc_groups=pbc_groups, typing_dgrid=typing_dgrid
    )

    sorted_pids = sorted([int(p["id"]) for p in pockets], key=lambda pid: fps[pid].max_probe_radius_a)
    bands: list[list[int]] = [[sorted_pids[0]]]
    for pid in sorted_pids[1:]:
        prev = fps[bands[-1][-1]].max_probe_radius_a
        curr = fps[pid].max_probe_radius_a
        if curr - prev > ENV_CLUSTER_INITIAL_BAND_GAP_A:
            bands.append([])
        bands[-1].append(pid)

    print(f"Step 1: {len(bands)} initial bands")
    for bi, grp in enumerate(bands):
        mr = [fps[p].max_probe_radius_a for p in grp]
        d1 = [next(p["d1_A"] for p in pockets if p["id"] == pid) for pid in grp]
        print(
            f"  band {bi:2d}: n={len(grp):3d}  max_r [{min(mr):.3f}, {max(mr):.3f}]  "
            f"d1 [{min(d1):.2f}, {max(d1):.2f}]  reps={grp[:8]}{'...' if len(grp)>8 else ''}"
        )

    bands_work = copy.deepcopy(bands)
    merge_log = []
    changed = True
    while changed:
        changed = False
        bands_work.sort(key=lambda g: mean_max_r(g, fps))
        for i in range(len(bands_work) - 1):
            ga, gb = bands_work[i], bands_work[i + 1]
            gap = mean_max_r(gb, fps) - mean_max_r(ga, fps)
            if gap > ENV_CLUSTER_MERGE_MAX_R_GAP_A:
                continue
            cos_vals = [
                fingerprint_cosine_similarity(fps[a], fps[b], angular_only=False)
                for a in ga
                for b in gb
            ]
            eff = (
                ENV_CLUSTER_MERGE_COSINE
                if gap <= ENV_CLUSTER_MERGE_GAP_THRESHOLD_A
                else ENV_CLUSTER_MERGE_COSINE_LARGE_GAP
            )
            mean_cos = float(np.mean(cos_vals))
            if mean_cos >= eff:
                merge_log.append((gap, mean_cos, eff, list(ga), list(gb)))
                bands_work[i] = sorted(ga + gb)
                del bands_work[i + 1]
                changed = True
                break

    print(f"\nStep 2: {len(merge_log)} adjacent-band merge(s)")
    for k, (gap, mc, eff, ga, gb) in enumerate(merge_log):
        print(f"  merge {k}: Δmean_max_r={gap:.3f}  mean_cos={mc:.3f}  thresh={eff:.2f}")
        print(f"    band A max_r ~{mean_max_r(ga, fps):.3f}  reps {ga}")
        print(f"    band B max_r ~{mean_max_r(gb, fps):.3f}  reps {gb}")

    n_types, rep_to_type, audit = assign_pocket_types_env_cluster(pockets, fps)
    typed = assign_pocket_types_from_segmentation(results, atoms, method="env_cluster")
    print(f"\nFinal: {n_types} types (expected Cu-BTC=3, ZIF-20=6)  audit merges={audit['n_merges']}")

    # type summary
    from collections import defaultdict

    by_type: dict[int, list] = defaultdict(list)
    for p in pockets:
        pid = int(p["id"])
        by_type[rep_to_type[pid]].append((pid, fps[pid].max_probe_radius_a, p["d1_A"]))
    for ti in sorted(by_type):
        rows = by_type[ti]
        mrs = [r[1] for r in rows]
        d1s = [r[2] for r in rows]
        print(
            f"  type {ti}: n={len(rows)}  max_r [{min(mrs):.3f},{max(mrs):.3f}]  "
            f"d1 [{min(d1s):.2f},{max(d1s):.2f}]  rep={min(r[0] for r in rows)}"
        )

    # large pockets
    print("\nLargest pockets (top max_r):")
    rows = sorted(
        [(fps[int(p["id"])].max_probe_radius_a, int(p["id"]), p["d1_A"], rep_to_type[int(p["id"])]) for p in pockets],
        reverse=True,
    )[:12]
    for mr, pid, d1, ti in rows:
        print(f"  rep {pid:3d}: type {ti}  max_r={mr:.3f}  d1={d1:.2f}")

    # cross-type cosines among top bands if multiple large bands merged
    if len(merge_log):
        for gap, mc, eff, ga, gb in merge_log:
            if mean_max_r(ga, fps) > 5.0 or mean_max_r(gb, fps) > 5.0:
                ra, rb = ga[0], gb[0]
                ch = {}
                a, b = fps[ra], fps[rb]
                if a.psd_vec is not None and b.psd_vec is not None:
                    ch["psd"] = float(_cosine(a.psd_vec, b.psd_vec))
                if len(a.global_sh_power):
                    ch["sh"] = float(_cosine(a.global_sh_power, b.global_sh_power))
                if len(a.global_radial):
                    ch["radial"] = float(_cosine(a.global_radial, b.global_radial))
                if a.conn_sh_power is not None and b.conn_sh_power is not None:
                    ch["conn"] = float(_cosine(a.conn_sh_power, b.conn_sh_power))
                print(f"  large-merge pair rep {ra} vs {rb} channels: {ch}")
                print(f"    direct overall cos = {fingerprint_cosine_similarity(a,b):.4f}")


def main():
    slug = sys.argv[1] if len(sys.argv) > 1 else "cubtc"
    backend = "rtree" if slug == "zif20" else "mic"
    trace(slug.upper(), CIFS[slug], backend)


if __name__ == "__main__":
    main()

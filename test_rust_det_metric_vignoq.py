#!/usr/bin/env python3
"""Validate rust_det_metric full pipeline on VIGNOQ_clean (PLD, boundaries, windows)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from ase.io import read

import mofography as mgr
import test_metric_watershed as tmw
from mofography.compute import calculate_padding_thickness, find_pore_limiting_radius
from mofography.rust_fun import dgrid_from_atoms_rust, do_periodic_segmentation_rust

ROOT = Path(__file__).resolve().parent
VIGNOQ_CIF = ROOT / "examples" / "TestCIFs" / "VIGNOQ_clean.cif"
DEFAULT_HTML = ROOT / "vignoq_rust_det_metric.html"

SPACING = 0.5
MASK_THICKNESS = 0.05
H_PARAM = 0.2
EXTRA_PADDING_ANGSTROM = 1.0
PLD_TARGET = 1.378
PLD_TOL = 0.05
SELF_LOOP_MAX_TARGET = 1.38
SELF_LOOP_MAX_TOL = 0.15


def _is_vignoq_cif(cif: Path) -> bool:
    return cif.name.upper().startswith("VIGNOQ")


def _check_maxima_region_indexing(seg: dict) -> None:
    """maxima_info label L must equal regions[pos] for in-cell seeds (reindex contract)."""
    regions = np.asarray(seg["regions"])
    mi = np.asarray(seg["maxima_info"], dtype=int)
    shape = np.array(regions.shape)
    mism = 0
    for row in mi:
        lab = int(row[0])
        pos = tuple((row[1:4] + row[4:7] * shape).astype(int))
        if any(p < 0 or p >= regions.shape[i] for i, p in enumerate(pos)):
            continue
        if int(regions[pos]) != lab:
            mism += 1
    if mism:
        raise RuntimeError(
            f"{mism} maxima_info labels disagree with regions at seed voxels"
        )


def _boundary_stats_on_regions(
    regions_padded: np.ndarray,
    dgrid: np.ndarray,
    padding: int,
    pocket1: set[int],
    pocket2: set[int],
) -> tuple[dict, dict]:
    padded = np.pad(dgrid, padding, mode="wrap")
    void = np.pad(dgrid > MASK_THICKNESS, padding, mode="wrap")
    offsets_half = tmw.positive_half_offsets(tmw.neighbor_offsets_26())
    self_stats = tmw.summarize(
        tmw.boundary_values_for_label_sets(
            regions_padded, padded, void, offsets_half, pocket1, None
        )
    )
    inter_stats = tmw.summarize(
        tmw.boundary_values_for_label_sets(
            regions_padded, padded, void, offsets_half, pocket1, pocket2
        )
    )
    return self_stats, inter_stats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cif", type=Path, default=VIGNOQ_CIF)
    parser.add_argument(
        "--html",
        type=Path,
        nargs="?",
        const=DEFAULT_HTML,
        default=None,
        help="Write Plotly HTML (default path if flag alone)",
    )
    args = parser.parse_args()

    atoms = read(str(args.cif.expanduser().resolve()))
    radii = mgr.default_vdw_radii_angstrom(atoms)
    dgrid = dgrid_from_atoms_rust(atoms, spacing=SPACING, radii=radii).astype(np.float32)
    padding = int(
        calculate_padding_thickness(dgrid, atoms, SPACING)
        + np.ceil(EXTRA_PADDING_ANGSTROM / SPACING)
    )

    plr = float(find_pore_limiting_radius(dgrid))
    print(f"PLD (dgrid percolation): {plr:.4f} A")

    metric = do_periodic_segmentation_rust(
        dgrid,
        atoms,
        spacing=SPACING,
        mask_thickness=MASK_THICKNESS,
        h=H_PARAM,
        extra_padding=EXTRA_PADDING_ANGSTROM,
        compute_connections=True,
        return_padded_regions=True,
        watershed_backend="rust_det_metric",
    )
    regions = np.asarray(metric["regions"], dtype=np.uint32)
    regions_padded = np.asarray(metric["regions_padded"], dtype=np.uint32)
    pbc_groups = [set(int(x) for x in g) for g in metric["pbc_groups"]]
    pocket1, pocket2 = pbc_groups[0], pbc_groups[1]

    _check_maxima_region_indexing(metric)
    print(
        f"Indexing: regions.max={int(regions.max())}, "
        f"unique_regions={len(metric['unique_regions'])}, "
        f"maxima_info rows={len(metric['maxima_info'])}"
    )

    self_stats, inter_stats = _boundary_stats_on_regions(
        regions_padded, dgrid, padding, pocket1, pocket2
    )
    print(
        "Boundary self-loop (padded grid diagnostic): "
        f"max={self_stats['max']:.4f} A"
    )
    print(
        "Boundary inter-pocket: "
        f"max={inter_stats['max']:.4f} A"
    )

    vignoq = _is_vignoq_cif(args.cif)
    if vignoq and abs(plr - PLD_TARGET) > PLD_TOL:
        raise SystemExit(
            f"FAIL: PLD {plr:.4f} A outside [{PLD_TARGET - PLD_TOL:.3f}, "
            f"{PLD_TARGET + PLD_TOL:.3f}] A"
        )

    conns = metric.get("connections") or []
    if not conns:
        raise SystemExit("FAIL: no connections computed")

    self_mins = []
    all_mins = []
    for c in conns:
        d = np.asarray(c[-1], dtype=float)
        if d.size == 0:
            continue
        b = float(np.min(d))
        all_mins.append(b)
        pair = [int(x) for x in c[0]]
        if pair[0] in pocket1 and pair[1] in pocket1:
            self_mins.append(b)
    print(
        "Connection bottlenecks: "
        f"max={max(all_mins):.4f} A, "
        f"self-loop max={max(self_mins) if self_mins else float('nan'):.4f} A "
        f"(n={len(self_mins)})"
    )
    if vignoq:
        if not self_mins:
            raise SystemExit("FAIL: no self-loop connections in pocket-1")
        if abs(max(self_mins) - SELF_LOOP_MAX_TARGET) > SELF_LOOP_MAX_TOL:
            raise SystemExit(
                f"FAIL: self-loop connection bottleneck {max(self_mins):.4f} A "
                f"outside ~{SELF_LOOP_MAX_TARGET} +/- {SELF_LOOP_MAX_TOL} A"
            )
        print("PASS: rust_det_metric pipeline (PLD, indexing, connections)")
    else:
        print("PASS: rust_det_metric pipeline (indexing, connections)")


    if args.html is not None:
        tmw.write_visualization_html(
            out_path=args.html,
            mat_atoms=atoms,
            dgrid=dgrid,
            seg=metric,
            regions_standard=regions,
            regions_metric=regions,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())

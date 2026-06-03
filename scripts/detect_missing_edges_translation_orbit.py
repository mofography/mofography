#!/usr/bin/env python
"""Translation-orbit missing-edge detector + HTML visualizer.

This version does NOT depend on identifying an asymmetric unit or symmetry ops.
It detects incompletely occupied translation shells for each base-node pair and
adds synthetic internal connections for missing orbit members.

Usage:
  conda run -n mgr_env python scripts/detect_missing_edges_translation_orbit.py \
    examples/TestCIFs/NU-1000.cif /tmp/nu1000_translation_orbit.html
"""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path
from itertools import product

import numpy as np

MOFOGRAPHY_ROOT = Path("/Users/arungopalan/repos/mofography")
TIPMOF_ROOT = Path.home() / "repos/tipmof"
sys.path.insert(0, str(MOFOGRAPHY_ROOT))
sys.path.insert(0, str(TIPMOF_ROOT))

from mofography.cif_io import read_cif_to_ase
import mofography as mgr
from tipmof.tipmof import collect_periodic_edges


def canonical_key(i: int, j: int, offset: tuple[int, int, int]):
    if i < j:
        return (i, j, offset)
    if i > j:
        return (j, i, tuple(-x for x in offset))
    if offset == (0, 0, 0):
        return None
    neg = tuple(-x for x in offset)
    return (i, i, offset) if offset < neg else (i, i, neg)


def cart_dist(fi: np.ndarray, fj: np.ndarray, off: tuple[int, int, int], cell: np.ndarray) -> float:
    d = fj + np.array(off, float) - fi
    return float(np.linalg.norm(d @ cell.T))


def cluster_shell(items, tol: float):
    if not items:
        return []
    s = sorted(items, key=lambda x: x[1])
    shells = []
    cur = [s[0]]
    ref = s[0][1]
    for key, d in s[1:]:
        if ref > 0 and abs(d - ref) / ref < tol:
            cur.append((key, d))
        else:
            shells.append(cur)
            cur = [(key, d)]
            ref = d
    shells.append(cur)
    return shells


def build_tipmof_base_mapping(maxima_info: np.ndarray):
    maxima_pbc = np.vstack([m[1:4] for m in maxima_info if np.all(m[4:7] == 0)])
    maxima_pbc_dict = {tuple(m): i + 1 for i, m in enumerate(maxima_pbc)}
    return np.array([maxima_pbc_dict[tuple(m[1:4])] for m in maxima_info], dtype=int)


def build_base_offset_label_index(maxima_info: np.ndarray, orig_map: np.ndarray):
    """Map base -> list[(offset, label)] plus lookup (base,offset)->labels."""
    base_to_offset_labels = {}
    by_key = {}
    for idx, row in enumerate(maxima_info):
        base = int(orig_map[idx])
        off = tuple(int(x) for x in row[4:7])
        label = int(row[0])
        base_to_offset_labels.setdefault(base, []).append((off, label))
        by_key.setdefault((base, off), []).append(label)
    for b in base_to_offset_labels:
        base_to_offset_labels[b] = sorted(base_to_offset_labels[b], key=lambda x: (x[0], x[1]))
    for k in by_key:
        by_key[k] = sorted(by_key[k])
    return base_to_offset_labels, by_key


def build_close_inside_label_set(maxima_info: np.ndarray, dgrid_shape, close_or_inside_maxima: np.ndarray):
    """Return maxima labels that correspond to close_or_inside maxima positions.

    Compare in *global voxel* space: row[1:4] + row[4:7]*shape.
    """
    shape = np.array(dgrid_shape, int)
    ci_set = {tuple(int(x) for x in p) for p in np.asarray(close_or_inside_maxima, int)}
    labels = set()
    for row in maxima_info:
        vox = row[1:4].astype(int)
        off = row[4:7].astype(int)
        g = tuple(int(x) for x in (vox + off * shape))
        if g in ci_set:
            labels.add(int(row[0]))
    return labels


def choose_label_pair(
    base_i: int,
    base_j: int,
    rel_off: tuple[int, int, int],
    base_to_offset_labels,
    by_key,
    preferred_labels: set[int] | None = None,
):
    """Choose concrete maxima labels (li, lj) satisfying oj - oi = rel_off.

    Preference order:
      1) both labels in preferred set (close/inside labels)
      2) minimal offset magnitude
      3) deterministic by label ids
    """
    cands = []
    pref = preferred_labels or set()

    for oi, li in base_to_offset_labels.get(base_i, []):
        oj = (oi[0] + rel_off[0], oi[1] + rel_off[1], oi[2] + rel_off[2])
        labels_j = by_key.get((base_j, oj))
        if not labels_j:
            continue
        lj = labels_j[0]
        score_off = abs(oi[0]) + abs(oi[1]) + abs(oi[2]) + abs(oj[0]) + abs(oj[1]) + abs(oj[2])
        in_pref = int(li in pref and lj in pref)
        # sort key: prefer in_pref=1 first (negative), then smaller offsets, then labels
        cands.append((-in_pref, score_off, li, lj, oi, oj))

    if not cands:
        return None
    cands.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    _, _, li, lj, oi, oj = cands[0]
    return li, lj, oi, oj


def run(
    cif_path: str,
    output_html: str,
    offset_range: int = 1,
    shell_tol: float = 0.15,
    cutoff_scale: float = 1.05,
):
    atoms = read_cif_to_ase(cif_path)
    atoms = mgr.update_with_covalent_bonds(atoms)
    mgr.compute_aromatic_cycles(atoms)

    dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=0.5)
    seg = mgr.do_periodic_segmentation_rust(
        dgrid,
        atoms,
        spacing=0.5,
        mask_thickness=0.05,
        h=0.2,
        threshold_abs=1.5,
        min_distance=1.5,
        extra_padding=5,
        compute_connections=True,
        return_dgrid=True,
    )

    from ase.geometry import cellpar_to_cell

    cell = cellpar_to_cell(atoms.cell.cellpar())
    maxima_info = np.asarray(seg["maxima_info"])
    orig_map = build_tipmof_base_mapping(maxima_info)

    preferred_labels = build_close_inside_label_set(
        maxima_info, dgrid.shape, np.asarray(seg["close_or_inside_maxima"])
    )

    edge_map, _ = collect_periodic_edges(seg["connections"], maxima_info)
    existing = set(edge_map.keys())

    is_base = np.all(maxima_info[:, 4:7] == 0, axis=1)
    base_nodes = sorted(set(int(orig_map[i]) for i in np.where(is_base)[0]))

    frac_pos = {}
    for b in base_nodes:
        idx = int(np.where(orig_map == b)[0][0])
        row = maxima_info[idx]
        frac_pos[b] = np.round(row[1:4] / np.array(dgrid.shape) + row[4:7], 6)

    pair_edges = {}
    for i, j, off in existing:
        pair_edges.setdefault((i, j), set()).add(off)

    candidates_offsets = list(product(range(-offset_range, offset_range + 1), repeat=3))

    missing = []
    for (i, j), ex_offs in sorted(pair_edges.items()):
        fi, fj = frac_pos[i], frac_pos[j]
        ex_dists = [cart_dist(fi, fj, off, cell) for off in ex_offs]
        cutoff = max(ex_dists) * cutoff_scale

        cands = []
        for off in candidates_offsets:
            k = canonical_key(i, j, tuple(int(x) for x in off))
            if k is None:
                continue
            d = cart_dist(fi, fj, k[2], cell)
            if d <= cutoff:
                cands.append((k, d))

        for shell in cluster_shell(cands, tol=shell_tol):
            skeys = {k for k, _ in shell}
            s_exist = skeys & existing
            s_miss = skeys - existing
            if not s_exist or not s_miss:
                continue
            for k, d in shell:
                if k in s_miss:
                    missing.append(
                        {
                            "edge": k,
                            "dist": d,
                            "shell_ex": len(s_exist),
                            "shell_miss": len(s_miss),
                            "shell_total": len(shell),
                        }
                    )

    seen = set()
    missing_u = []
    for m in missing:
        if m["edge"] in seen:
            continue
        seen.add(m["edge"])
        missing_u.append(m)

    base_to_offset_labels, by_key = build_base_offset_label_index(maxima_info, orig_map)

    mapped = []
    for m in missing_u:
        i, j, off = m["edge"]
        pair = choose_label_pair(i, j, off, base_to_offset_labels, by_key, preferred_labels)
        if pair is None:
            continue
        li, lj, oi, oj = pair
        mm = dict(m)
        mm.update({"label_i": li, "label_j": lj, "off_i": oi, "off_j": oj})
        mapped.append(mm)

    print(f"Base nodes: {len(base_nodes)} | Existing edges: {len(existing)}")
    print(f"Preferred close/inside labels found: {len(preferred_labels)}")
    print(f"Missing translation-orbit members: {len(mapped)}")

    # ---- summary tables ----
    dz0 = sum(1 for m in mapped if m["edge"][2][2] == 0)
    dznz = len(mapped) - dz0
    print("\nSummary: offset class")
    print("  class        count")
    print("  -----------  -----")
    print(f"  Δz = 0       {dz0:>5}")
    print(f"  Δz ≠ 0       {dznz:>5}")

    # pair-level summary
    pair_stats = {}
    for m in mapped:
        i, j, off = m["edge"]
        key = (i, j)
        ps = pair_stats.setdefault(key, {"count": 0, "best_occ": 1.0, "best_sig": None})
        ps["count"] += 1
        occ = m["shell_ex"] / max(1, m["shell_total"])
        if occ < ps["best_occ"]:
            ps["best_occ"] = occ
            ps["best_sig"] = (off, m["shell_ex"], m["shell_total"], m["dist"])

    print("\nSummary: top pairs by missing-edge count")
    print("  pair      miss_count  worst_shell_occ  representative")
    print("  --------  ----------  ---------------  ------------------------------")
    for (i, j), ps in sorted(pair_stats.items(), key=lambda kv: (-kv[1]["count"], kv[0]))[:10]:
        sig = ps["best_sig"]
        if sig is None:
            rep = "-"
            occ_txt = "-"
        else:
            off, exn, tot, d = sig
            rep = f"off={off} d={d:.2f}A"
            occ_txt = f"{exn}/{tot}"
        print(f"  ({i:>2},{j:>2})     {ps['count']:>6}      {occ_txt:>7}         {rep}")

    # lowest-occupancy shells (most suspicious)
    ranked = sorted(
        mapped,
        key=lambda m: (m["shell_ex"] / max(1, m["shell_total"]), m["dist"])
    )
    print("\nSummary: lowest-occupancy candidate shells")
    print("  occ    base_edge                    dist")
    print("  -----  ---------------------------  -------")
    for m in ranked[:12]:
        i, j, off = m["edge"]
        occ = m["shell_ex"] / max(1, m["shell_total"])
        print(f"  {m['shell_ex']}/{m['shell_total']:<4}  ({i:>2},{j:>2},{str(off):<13})   {m['dist']:>6.2f}A")

    print("\nDetailed missing edges")
    for m in sorted(mapped, key=lambda x: x["dist"]):
        i, j, off = m["edge"]
        print(
            f"  base({i},{j},{off}) |d|={m['dist']:.2f}A "
            f"[{m['shell_ex']}ex/{m['shell_total']} shell] -> "
            f"labels {m['label_i']}-{m['label_j']} (offs {m['off_i']}->{m['off_j']})"
        )

    seg_aug = copy.deepcopy(seg)
    conns = list(seg_aug["connections"])

    existing_internal_pairs = {
        tuple(sorted((int(c[0][0]), int(c[0][1]))))
        for c in conns
        if len(c) >= 3 and bool(c[1]) and str(c[2]) == "internal"
    }

    added = 0
    for m in mapped:
        li, lj = int(m["label_i"]), int(m["label_j"])
        pair = tuple(sorted((li, lj)))
        if pair in existing_internal_pairs:
            continue

        idx_i = np.where(maxima_info[:, 0].astype(int) == li)[0]
        idx_j = np.where(maxima_info[:, 0].astype(int) == lj)[0]
        if len(idx_i) == 0 or len(idx_j) == 0:
            continue

        p_i = maxima_info[idx_i[0], 1:4] / (np.array(dgrid.shape) - 1.0)
        p_j = maxima_info[idx_j[0], 1:4] / (np.array(dgrid.shape) - 1.0)
        mid_frac = ((p_i + p_j) / 2.0) % 1.0
        mid_cart = mid_frac @ cell.T

        conns.append(
            np.array(
                [
                    np.array([li, lj], dtype=int),
                    True,
                    "internal",
                    np.array([mid_cart], dtype=float),
                    np.array([2.0], dtype=float),
                ],
                dtype=object,
            )
        )
        existing_internal_pairs.add(pair)
        added += 1

    seg_aug["connections"] = conns
    print(f"Added synthetic internal edges: {added}")

    out = mgr.show_mof_threejs(
        seg_aug,
        atoms,
        plot_framework_representation="scatter",
        plot_surface=True,
        plot_surface_step_size=2,
        plot_connections=True,
        plot_maxima=True,
        plot_windows=False,
        plot_cell=True,
        output_html=output_html,
    )
    print(f"Saved HTML: {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cif", help="Input CIF file")
    ap.add_argument("output", nargs="?", default=None, help="Output HTML path")
    ap.add_argument("--offset-range", type=int, default=1)
    ap.add_argument("--shell-tol", type=float, default=0.15)
    ap.add_argument("--cutoff-scale", type=float, default=1.05)
    args = ap.parse_args()

    cif = str(Path(args.cif).resolve())
    out = args.output or f"/tmp/{Path(cif).stem}_translation_orbit_missing.html"
    run(cif, out, offset_range=args.offset_range, shell_tol=args.shell_tol, cutoff_scale=args.cutoff_scale)


if __name__ == "__main__":
    main()

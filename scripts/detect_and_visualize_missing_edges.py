#!/usr/bin/env python
"""Detect missing pore-graph edges and visualize with synthetic internal links.

Key point: connections[] uses maxima *labels* (from maxima_info[:,0]), not
close_or_inside region ids. This script maps missing base-edge keys to concrete
maxima-label pairs (e.g., 13-16) so they appear in the HTML plot.
"""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

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


def cluster_shell(items: list[tuple[tuple[int, int, tuple[int, int, int]], float]], tol: float = 0.15):
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
    """Reproduce collect_periodic_edges base-node mapping exactly."""
    maxima_pbc = np.vstack([m[1:4] for m in maxima_info if np.all(m[4:7] == 0)])
    maxima_pbc_dict = {tuple(m): i + 1 for i, m in enumerate(maxima_pbc)}
    orig_map = np.array([maxima_pbc_dict[tuple(m[1:4])] for m in maxima_info], dtype=int)
    return orig_map


def build_base_offset_to_label(maxima_info: np.ndarray, orig_map: np.ndarray):
    """Map (base_label, offset) -> maxima label used by connections[]."""
    out: dict[tuple[int, tuple[int, int, int]], int] = {}
    for idx, row in enumerate(maxima_info):
        base = int(orig_map[idx])
        off = tuple(int(x) for x in row[4:7])
        label = int(row[0])
        key = (base, off)
        if key not in out or label < out[key]:
            out[key] = label
    return out


def map_close_inside_to_base_offsets(
    maxima_info: np.ndarray,
    orig_map: np.ndarray,
    dgrid_shape: tuple[int, int, int],
    close_or_inside_maxima: np.ndarray,
):
    """Map close_or_inside maxima coords to (base, offset), then to maxima labels.

    close_or_inside_maxima may include coords outside [0, shape), so we solve:
      ci_pos ~= base_voxel + offset * dgrid_shape
    using base rows (offset=0) as anchors.
    """
    base_rows = []
    for idx, row in enumerate(maxima_info):
        if np.all(row[4:7] == 0):
            base_rows.append((int(orig_map[idx]), row[1:4].astype(int)))

    shape = np.array(dgrid_shape, int)
    mapped: dict[tuple[int, tuple[int, int, int]], tuple[int, int, int]] = {}

    for ci_pos in close_or_inside_maxima:
        ci = np.array(ci_pos, dtype=float)
        best = None
        for base, vox in base_rows:
            off_f = (ci - vox) / shape
            off = np.round(off_f).astype(int)
            recon = vox + off * shape
            err = float(np.linalg.norm(recon - ci))
            if err < 1e-6:
                score = int(np.abs(off).sum())
                cand = (score, base, tuple(int(x) for x in off), tuple(int(x) for x in vox))
                if best is None or cand < best:
                    best = cand
        if best is not None:
            _, base, off, vox = best
            mapped[(base, off)] = vox
    return mapped


def choose_close_inside_label_pair(
    base_i: int,
    base_j: int,
    rel_off: tuple[int, int, int],
    ci_base_offsets: set[tuple[int, tuple[int, int, int]]],
    base_offset_to_label: dict[tuple[int, tuple[int, int, int]], int],
):
    """Pick maxima-label endpoints for a missing base-edge using close/inside offsets.

    We search all close/inside offsets oi for base_i and require oj = oi + rel_off
    to also be a close/inside offset for base_j.
    """
    cands = []
    for b, oi in ci_base_offsets:
        if b != base_i:
            continue
        oj = (oi[0] + rel_off[0], oi[1] + rel_off[1], oi[2] + rel_off[2])
        if (base_j, oj) not in ci_base_offsets:
            continue
        li = base_offset_to_label.get((base_i, oi))
        lj = base_offset_to_label.get((base_j, oj))
        if li is None or lj is None:
            continue
        score = abs(oi[0]) + abs(oi[1]) + abs(oi[2]) + abs(oj[0]) + abs(oj[1]) + abs(oj[2])
        cands.append((score, li, lj, oi, oj))

    if not cands:
        return None
    cands.sort(key=lambda x: (x[0], x[1], x[2]))
    _, li, lj, oi, oj = cands[0]
    return li, lj, oi, oj


def run(cif_path: str, output_html: str | None):
    cif = str(Path(cif_path).resolve())
    if output_html is None:
        output_html = f"/tmp/{Path(cif).stem}_missing_edges.html"

    atoms = read_cif_to_ase(cif)
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
    base_offset_to_label = build_base_offset_to_label(maxima_info, orig_map)

    ci_base_map = map_close_inside_to_base_offsets(
        maxima_info=maxima_info,
        orig_map=orig_map,
        dgrid_shape=dgrid.shape,
        close_or_inside_maxima=np.asarray(seg["close_or_inside_maxima"]),
    )
    ci_base_offsets = set(ci_base_map.keys())

    edge_map, _ = collect_periodic_edges(seg["connections"], maxima_info)
    existing = set(edge_map.keys())

    is_base = np.all(maxima_info[:, 4:7] == 0, axis=1)
    base_labels = sorted(int(orig_map[i]) for i in np.where(is_base)[0])
    base_labels = sorted(set(base_labels))

    frac_pos = {}
    for node in base_labels:
        idx = int(np.where(orig_map == node)[0][0])
        row = maxima_info[idx]
        frac_pos[node] = np.round(row[1:4] / np.array(dgrid.shape) + row[4:7], 3)

    pair_edges: dict[tuple[int, int], set[tuple[int, int, int]]] = {}
    for i, j, off in existing:
        pair_edges.setdefault((i, j), set()).add(off)

    missing = []
    for (i, j), ex_offs in sorted(pair_edges.items()):
        fi, fj = frac_pos[i], frac_pos[j]
        ex_dists = [cart_dist(fi, fj, off, cell) for off in ex_offs]
        cutoff = max(ex_dists) * 1.05

        cands = []
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                for dk in (-1, 0, 1):
                    k = canonical_key(i, j, (di, dj, dk))
                    if k is None:
                        continue
                    d = cart_dist(fi, fj, k[2], cell)
                    if d <= cutoff:
                        cands.append((k, d))

        shells = cluster_shell(cands, tol=0.15)
        for shell in shells:
            skeys = {k for k, _ in shell}
            shell_exist = skeys & existing
            shell_miss = skeys - existing
            if not shell_exist or not shell_miss:
                continue
            for k, d in shell:
                if k not in shell_miss:
                    continue
                pair = choose_close_inside_label_pair(
                    k[0], k[1], k[2], ci_base_offsets, base_offset_to_label
                )
                if pair is None:
                    continue
                li, lj, oi, oj = pair
                missing.append(
                    {
                        "edge": k,
                        "dist": d,
                        "label_i": li,
                        "label_j": lj,
                        "off_i": oi,
                        "off_j": oj,
                        "shell_ex": len(shell_exist),
                        "shell_miss": len(shell_miss),
                    }
                )

    # dedupe by base-edge key + label pair
    seen = set()
    missing_u = []
    for m in missing:
        key = (m["edge"], tuple(sorted((m["label_i"], m["label_j"]))))
        if key in seen:
            continue
        seen.add(key)
        missing_u.append(m)

    print(f"Missing close/inside edges found: {len(missing_u)}")
    for m in sorted(missing_u, key=lambda x: x["dist"]):
        i, j, off = m["edge"]
        print(
            f"  base({i},{j},{off}) |d|={m['dist']:.2f}A -> labels {m['label_i']}-{m['label_j']} "
            f"(offs {m['off_i']}->{m['off_j']}) [{m['shell_ex']}ex+{m['shell_miss']}mi]"
        )

    seg_aug = copy.deepcopy(seg)
    conns = list(seg_aug["connections"])

    # avoid duplicate exact internal connection labels
    existing_internal_pairs = {
        tuple(sorted((int(c[0][0]), int(c[0][1]))))
        for c in conns
        if len(c) >= 3 and bool(c[1]) and str(c[2]) == "internal"
    }

    added = 0
    for m in missing_u:
        li, lj = int(m["label_i"]), int(m["label_j"])
        pair = tuple(sorted((li, lj)))
        if pair in existing_internal_pairs:
            continue

        # midpoint placeholder
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
                    "internal",  # required for plotting
                    np.array([mid_cart], dtype=float),
                    np.array([2.0], dtype=float),
                ],
                dtype=object,
            )
        )
        existing_internal_pairs.add(pair)
        added += 1

    seg_aug["connections"] = conns
    print(f"Added synthetic internal connections: {added}")

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
    p = argparse.ArgumentParser()
    p.add_argument("cif", help="Input CIF")
    p.add_argument("output", nargs="?", default=None, help="Output HTML path")
    args = p.parse_args()
    run(args.cif, args.output)


if __name__ == "__main__":
    main()

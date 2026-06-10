#!/usr/bin/env python3
"""Prototype window-type report for LTA-96 (zeolite).

Clusters periodic graph windows by diameter + shell fingerprint, writes a
4-panel PDF (counts, radius bars, element-shell bars, cosine matrix) and a
Plotly HTML with orthographic framework subsets around representative windows.

  python examples/window_types_lta96_report.py
  python examples/window_types_lta96_report.py --cif path/to/LTA96.cif --out-dir examples/viz_outputs

Modeled on mgr_sn/mgr_org_scripts/pocket_types_report.py; uses mofography only.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from ase.data import chemical_symbols
from ase.io import read
from matplotlib.gridspec import GridSpec

import mofography as mgr
from mofography.pocket_typing import DEFAULT_MASK_THICKNESS_A, DEFAULT_SPACING_A
from mofography.winding_percolation import (
    _cosine_similarity_safe,
    annotate_graph_edges_with_window_fingerprints,
    make_directed_graph_from_segmentation_results,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CIF = REPO_ROOT / "examples/TestCIFs/LTA96.cif"
TYPE_COLORS = ["#4C72B0", "#55A868", "#DD8452", "#C44E52", "#8172B2", "#937860", "#DA8BC3"]
SEG_KW = dict(
    mask_thickness=DEFAULT_MASK_THICKNESS_A,
    h=0.2,
    min_distance=1.5,
    threshold_abs=1.5,
    extra_padding=5.0,
    return_dgrid=True,
)
SHELL_RADIUS_A = 3.5
RADIUS_BAND_GAP_A = 0.18
MERGE_COSINE = 0.88
BALL_QUERY_A = 5.0
VIEW_DIR = (0.0, 0.0, 1.0)


@dataclass
class WindowRecord:
    """One canonical graph window instance."""

    wid: int
    u: int
    v: int
    key: int
    window_radius_a: float
    cart_pos: np.ndarray
    fp_vec: np.ndarray
    pocket_radii: Tuple[float, float] = (0.0, 0.0)


@dataclass
class WindowTypingResult:
    ok: bool
    error: str = ""
    windows: List[WindowRecord] = field(default_factory=list)
    rep_to_type: Dict[int, int] = field(default_factory=dict)
    type_summary: List[Dict[str, Any]] = field(default_factory=list)
    pairwise_cosine: Dict[str, float] = field(default_factory=dict)


def _pair_key(ti: int, tj: int) -> str:
    a, b = min(ti, tj), max(ti, tj)
    return f"{a}_vs_{b}"


def _window_cart_from_edge(data: dict, shape: Tuple[int, int, int], atoms) -> Optional[np.ndarray]:
    positions = data.get("window_positions") or ()
    if not positions:
        return None
    radii = data.get("window_radii") or ()
    pick = 0
    if radii and len(radii) == len(positions):
        pick = int(np.argmax(np.asarray(radii, dtype=float)))
    try:
        idx = np.asarray(positions[pick], dtype=float)
    except Exception:
        return None
    if idx.shape != (3,):
        return None
    return np.asarray(
        mgr.find_coord_from_indices(idx, shape=shape, ase_atoms=atoms),
        dtype=float,
    ).reshape(3)


def collect_window_records(results: dict, atoms, shell_radius: float = SHELL_RADIUS_A) -> List[WindowRecord]:
    graph = make_directed_graph_from_segmentation_results(results)
    annotate_graph_edges_with_window_fingerprints(
        graph, results, atoms=atoms, shell_radius=float(shell_radius)
    )
    shape = tuple(int(x) for x in np.asarray(results["dgrid"]).shape)
    seen: set = set()
    records: List[WindowRecord] = []
    wid = 0
    for u, v, key, data in graph.edges(keys=True, data=True):
        if not data.get("window_fp_ready", False):
            continue
        rel = tuple(data.get("relative_offset", (0, 0, 0)))
        sig = (min(int(u), int(v)), max(int(u), int(v)), rel, round(float(data.get("window_radius", 0.0)), 3))
        if sig in seen:
            continue
        seen.add(sig)
        cart = _window_cart_from_edge(data, shape, atoms)
        if cart is None:
            continue
        fp = np.asarray(data.get("window_fp_vec"), dtype=float)
        pr_u = float(graph.nodes[u].get("pocket_radius", 0.0))
        pr_v = float(graph.nodes[v].get("pocket_radius", 0.0))
        records.append(
            WindowRecord(
                wid=wid,
                u=int(u),
                v=int(v),
                key=int(key),
                window_radius_a=float(data.get("window_radius", 0.0)),
                cart_pos=cart,
                fp_vec=fp,
                pocket_radii=(pr_u, pr_v),
            )
        )
        wid += 1
    return records


def assign_window_types(
    windows: Sequence[WindowRecord],
    *,
    radius_band_gap: float = RADIUS_BAND_GAP_A,
    merge_cosine: float = MERGE_COSINE,
) -> WindowTypingResult:
    if not windows:
        return WindowTypingResult(ok=False, error="no fingerprinted windows")

    order = sorted(range(len(windows)), key=lambda i: windows[i].window_radius_a)
    parent = list(range(len(windows)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    for ai, i in enumerate(order):
        wi = windows[i]
        for j in order[ai + 1 :]:
            wj = windows[j]
            if wj.window_radius_a - wi.window_radius_a > radius_band_gap:
                break
            if _cosine_similarity_safe(wi.fp_vec, wj.fp_vec) >= merge_cosine:
                union(i, j)

    clusters: Dict[int, List[int]] = defaultdict(list)
    for i in range(len(windows)):
        clusters[find(i)].append(i)

    rep_to_type: Dict[int, int] = {}
    type_summary: List[Dict[str, Any]] = []
    reps_by_type: Dict[int, int] = {}
    sorted_clusters = sorted(clusters.items(), key=lambda kv: np.mean([windows[i].window_radius_a for i in kv[1]]))

    for tid, (_root, members) in enumerate(sorted_clusters):
        radii = [windows[i].window_radius_a for i in members]
        rep = min(members, key=lambda i: abs(windows[i].window_radius_a - float(np.median(radii))))
        reps_by_type[tid] = rep
        for i in members:
            rep_to_type[windows[i].wid] = tid
        wrep = windows[rep]
        type_summary.append(
            {
                "type_id": tid,
                "n_windows": len(members),
                "mean_window_radius_A": round(float(np.mean(radii)), 3),
                "std_window_radius_A": round(float(np.std(radii)), 3),
                "representative_wid": wrep.wid,
                "pocket_pair": f"{wrep.u}-{wrep.v}",
                "pocket_radii_A": [round(x, 2) for x in wrep.pocket_radii],
            }
        )

    types = sorted(reps_by_type)
    pairwise: Dict[str, float] = {}
    for i, ti in enumerate(types):
        for tj in types[i:]:
            ri, rj = reps_by_type[ti], reps_by_type[tj]
            pairwise[_pair_key(ti, tj)] = _cosine_similarity_safe(
                windows[ri].fp_vec, windows[rj].fp_vec
            )

    return WindowTypingResult(
        ok=True,
        windows=list(windows),
        rep_to_type=rep_to_type,
        type_summary=type_summary,
        pairwise_cosine=pairwise,
    )


def _top_element_channels(fp: np.ndarray, top_n: int = 6) -> List[Tuple[str, float]]:
    idx = np.argsort(fp)[::-1]
    out: List[Tuple[str, float]] = []
    for j in idx:
        if fp[j] <= 1e-6:
            break
        sym = chemical_symbols[int(j) + 1] if int(j) + 1 < len(chemical_symbols) else f"Z{j+1}"
        out.append((sym, float(fp[j])))
        if len(out) >= top_n:
            break
    return out


def plot_window_report_pdf(typed: WindowTypingResult, out_pdf: Path, material_name: str) -> None:
    windows = typed.windows
    types = sorted({typed.rep_to_type[w.wid] for w in windows})
    reps = {tid: next(w.wid for w in windows if typed.rep_to_type[w.wid] == tid) for tid in types}
    n_types = len(types)

    fig = plt.figure(figsize=(11, 9))
    gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.28)

    ax_a = fig.add_subplot(gs[0, 0])
    counts = [ts["n_windows"] for ts in typed.type_summary]
    labels = [f"T{ts['type_id']}\n{ts['mean_window_radius_A']:.2f} Å" for ts in typed.type_summary]
    ax_a.bar(range(n_types), counts, color=[TYPE_COLORS[i % len(TYPE_COLORS)] for i in range(n_types)])
    ax_a.set_xticks(range(n_types))
    ax_a.set_xticklabels(labels, fontsize=8)
    ax_a.set_ylabel("Window count")
    ax_a.set_title("(a) Window types by count", fontweight="bold")

    ax_b = fig.add_subplot(gs[0, 1])
    for i, ts in enumerate(typed.type_summary):
        tid = int(ts["type_id"])
        radii = [w.window_radius_a for w in windows if typed.rep_to_type[w.wid] == tid]
        ax_b.hist(
            radii,
            bins=min(12, max(4, len(radii))),
            alpha=0.55,
            color=TYPE_COLORS[i % len(TYPE_COLORS)],
            label=f"Type {tid}",
        )
    ax_b.set_xlabel("Window radius (Å)")
    ax_b.set_ylabel("Count")
    ax_b.legend(fontsize=8)
    ax_b.set_title("(b) Window radius distribution", fontweight="bold")

    ax_c = fig.add_subplot(gs[1, 0])
    bar_w = 0.8 / max(n_types, 1)
    elem_order: List[str] = []
    for tid in types:
        for sym, _ in _top_element_channels(windows[reps[tid]].fp_vec, top_n=5):
            if sym not in elem_order:
                elem_order.append(sym)
    x_pos = np.arange(len(elem_order))
    for i, tid in enumerate(types):
        fp = windows[reps[tid]].fp_vec
        weights = {sym: val for sym, val in _top_element_channels(fp, top_n=118)}
        heights = [weights.get(sym, 0.0) for sym in elem_order]
        off = (i - (n_types - 1) / 2) * bar_w
        ax_c.bar(
            x_pos + off,
            heights,
            bar_w * 0.9,
            color=TYPE_COLORS[i % len(TYPE_COLORS)],
            label=f"Type {tid}",
        )
    ax_c.set_ylabel("Norm. shell weight")
    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels(elem_order, fontsize=8)
    ax_c.legend(fontsize=7, loc="upper right")
    ax_c.set_title("(c) Framework shell (top elements)", fontweight="bold")

    ax_d = fig.add_subplot(gs[1, 1])
    mat = np.eye(n_types)
    for i, ti in enumerate(types):
        for j, tj in enumerate(types):
            if i != j:
                mat[i, j] = typed.pairwise_cosine.get(_pair_key(ti, tj), 0.0)
    im = ax_d.imshow(mat, vmin=0, vmax=1, cmap="RdYlGn")
    ax_d.set_xticks(range(n_types))
    ax_d.set_yticks(range(n_types))
    ax_d.set_xticklabels([f"T{t}" for t in types])
    ax_d.set_yticklabels([f"T{t}" for t in types])
    for i in range(n_types):
        for j in range(n_types):
            ax_d.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=9)
    plt.colorbar(im, ax=ax_d, fraction=0.046)
    ax_d.set_title("(d) Shell fingerprint cosine", fontweight="bold")

    fig.suptitle(f"Window types — {material_name}", fontweight="bold", y=0.98)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _framework_subset_trace(mat_updated, subset: Sequence[int], *, name: str, color_idx: int):
    """Bonded framework mesh for atom indices (avoids plot_subset_of_atoms API drift)."""
    if not subset:
        return None
    sg = mgr.sgraph_from_atoms(mat_updated).subgraph(subset)
    sg = nx.convert_node_labels_to_integers(sg, first_label=0)
    subatoms = mgr.get_subset_atoms(mat_updated, subset=list(subset))
    edges = list(sg.edges())
    if edges:
        subatoms.info["covalent_bonds"] = np.vstack(edges)
        subatoms.info["covalent_bond_orders"] = np.ones(len(edges), dtype=float)
    mesh = mgr.mesh_for_atoms(subatoms)
    trace = mgr.plot_mesh([mesh], stylize=True, name=name, opacity=0.92)
    trace.update(
        legendgroup=name,
        legendgrouptitle_text=name,
        lighting=dict(ambient=0.2, diffuse=0.9, specular=0.4),
    )
    return trace


def _atoms_near_points(atoms, points: np.ndarray, r_max: float) -> List[List[int]]:
    """Ball query: atom indices within r_max of each query point (freud AABB)."""
    aabb = mgr.AABB_on_atoms(atoms)
    out: List[List[int]] = []
    for pt in np.asarray(points, dtype=float).reshape(-1, 3):
        nl = aabb.query(
            pt.reshape(1, 3),
            query_args=dict(mode="ball", r_max=float(r_max)),
        ).toNeighborList(sort_by_distance=True)
        out.append(sorted(set(int(x) for x in nl.point_indices)))
    return out


def build_window_types_html(
    atoms,
    mat_updated,
    typed: WindowTypingResult,
    out_html: Path,
    ball_radius: float = BALL_QUERY_A,
) -> None:
    windows = typed.windows
    types = sorted({typed.rep_to_type[w.wid] for w in windows})
    reps = {}
    for tid in types:
        rep_wid = next(w.wid for w in windows if typed.rep_to_type[w.wid] == tid)
        reps[tid] = windows[rep_wid]

    traces: List[Any] = [mgr.plot_cell(atoms)]
    for i, tid in enumerate(types):
        w = reps[tid]
        color = TYPE_COLORS[i % len(TYPE_COLORS)]
        atom_lists = _atoms_near_points(mat_updated, w.cart_pos.reshape(1, 3), ball_radius)
        subset = atom_lists[0] if atom_lists else []
        fw_trace = _framework_subset_trace(
            mat_updated,
            subset,
            name=f"T{tid} framework",
            color_idx=i,
        )
        if fw_trace is not None:
            traces.append(fw_trace)
        traces.append(
            mgr.scatter3d(
                w.cart_pos.reshape(1, 3),
                name=f"T{tid} window ({w.window_radius_a:.2f} Å)",
                marker_color=color,
                marker_size=14,
                legendgroup=f"type{tid}",
            )
        )

    fig = go.Figure(data=traces, layout=mgr.dark_layout)
    mgr.make_fig_orthographic(fig)
    fig.update_layout(
        title=f"LTA window types — orthographic framework subsets (r={ball_radius:.1f} Å)",
        scene_camera=dict(eye=dict(x=VIEW_DIR[0], y=VIEW_DIR[1], z=VIEW_DIR[2])),
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(out_html))


def run_pipeline(cif_path: Path, spacing: float = DEFAULT_SPACING_A) -> Tuple[Any, Any, dict, WindowTypingResult]:
    atoms = read(str(cif_path))
    mat_updated = mgr.update_with_covalent_bonds(atoms)
    dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=spacing).astype(np.float32)
    results = mgr.do_periodic_segmentation_rust(dgrid, atoms, spacing=spacing, **SEG_KW)
    results["mat_atoms"] = atoms
    windows = collect_window_records(results, atoms)
    typed = assign_window_types(windows)
    if not typed.ok:
        raise RuntimeError(typed.error)
    return atoms, mat_updated, results, typed


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="LTA-96 window-type prototype report")
    parser.add_argument("--cif", type=Path, default=DEFAULT_CIF, help="CIF path (default: LTA96)")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "examples/viz_outputs")
    parser.add_argument("--spacing", type=float, default=DEFAULT_SPACING_A)
    args = parser.parse_args(argv)

    cif = args.cif.expanduser().resolve()
    if not cif.is_file():
        print(f"CIF not found: {cif}", file=sys.stderr)
        return 1

    stem = cif.stem.replace(" ", "_")
    out_dir = args.out_dir.expanduser().resolve()
    pdf_path = out_dir / f"{stem}_window_types_report.pdf"
    html_path = out_dir / f"{stem}_window_types_panel.html"

    print(f"Segmenting {cif.name} …")
    atoms, mat_updated, _results, typed = run_pipeline(cif, spacing=args.spacing)

    print(f"Detected {len(typed.type_summary)} window types from {len(typed.windows)} canonical windows:")
    for ts in typed.type_summary:
        print(
            f"  Type {ts['type_id']}: n={ts['n_windows']}, "
            f"r_win={ts['mean_window_radius_A']:.3f} Å, pockets {ts['pocket_pair']}"
        )

    plot_window_report_pdf(typed, pdf_path, material_name=stem)
    build_window_types_html(atoms, mat_updated, typed, html_path)
    print(f"Wrote {pdf_path}")
    print(f"Wrote {html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

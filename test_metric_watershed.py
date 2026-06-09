#!/usr/bin/env python3
"""
Prototype metric-aware bottleneck watershed on VIGNOQ_clean.

Compares:
1) Standard bottleneck priority flood (all neighbor weights = 1)
2) Metric-aware bottleneck priority flood (priorities scaled by metric edge lengths)

Writes optional Plotly HTML for visual inspection.
"""

from __future__ import annotations

import argparse
import heapq
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from ase.io import read

import mofography as mgr
from mofography.compute import calculate_padding_thickness, find_pore_limiting_radius
from mofography.rust_fun import dgrid_from_atoms_rust, do_periodic_segmentation_rust

Offset = Tuple[int, int, int]

ROOT = Path(__file__).resolve().parent
SPACING = 0.5
MASK_THICKNESS = 0.05
H_PARAM = 0.2
EXTRA_PADDING_ANGSTROM = 1.0
W_MIN = 0.5
W_MAX = 2.0

VIGNOQ_CIF = ROOT / "examples" / "TestCIFs" / "VIGNOQ_clean.cif"
DEFAULT_HTML = ROOT / "vignoq_metric_aware_segmentation.html"


def neighbor_offsets_26() -> List[Offset]:
    return [
        (dx, dy, dz)
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
        for dz in (-1, 0, 1)
        if not (dx == 0 and dy == 0 and dz == 0)
    ]


def positive_half_offsets(offsets: Sequence[Offset]) -> List[Offset]:
    out: List[Offset] = []
    for dx, dy, dz in offsets:
        if (dx > 0) or (dx == 0 and dy > 0) or (dx == 0 and dy == 0 and dz > 0):
            out.append((dx, dy, dz))
    return out


def flat_index(i: int, j: int, k: int, shape: Tuple[int, int, int]) -> int:
    _, ny, nz = shape
    return int(i * ny * nz + j * nz + k)


def compute_metric_inv_weights(
    cell_matrix: np.ndarray,
    grid_shape: Tuple[int, int, int],
    spacing: float,
    offsets: Sequence[Offset],
    w_min: float = W_MIN,
    w_max: float = W_MAX,
) -> Tuple[Dict[Offset, float], Dict[str, float]]:
    nx, ny, nz = grid_shape
    step = cell_matrix @ np.diag([1.0 / nx, 1.0 / ny, 1.0 / nz])
    metric = step.T @ step

    inv_w: Dict[Offset, float] = {}
    weights: List[float] = []
    for dx, dy, dz in offsets:
        n = np.array([dx, dy, dz], dtype=np.float64)
        l_actual = float(np.sqrt(n @ metric @ n))
        axes_touched = abs(dx) + abs(dy) + abs(dz)
        l0 = spacing * np.sqrt(float(axes_touched))
        w = float(np.clip(l_actual / l0, w_min, w_max))
        weights.append(w)
        inv_w[(dx, dy, dz)] = 1.0 / w

    w_arr = np.asarray(weights, dtype=np.float64)
    stats = {
        "w_min": float(np.min(w_arr)),
        "w_max": float(np.max(w_arr)),
        "w_mean": float(np.mean(w_arr)),
        "anisotropy_ratio": float(np.max(w_arr) / np.min(w_arr)),
    }
    return inv_w, stats


def build_unmerged_markers_from_maxima_info(
    maxima_info: np.ndarray, dgrid_shape: Tuple[int, int, int], padding: int
) -> np.ndarray:
    markers = np.zeros(
        tuple(np.array(dgrid_shape, dtype=int) + 2 * int(padding)), dtype=np.uint32
    )
    base_shape = np.array(dgrid_shape, dtype=int)

    for row in np.asarray(maxima_info, dtype=int):
        label = int(row[0])
        base = row[1:4]
        offs = row[4:7]
        pos = base + padding + offs * base_shape
        if np.any(pos < 0) or np.any(pos >= np.array(markers.shape)):
            continue
        markers[tuple(pos)] = np.uint32(label)
    return markers


def metric_aware_watershed(
    padded_dgrid: np.ndarray,
    markers: np.ndarray,
    inv_w: Dict[Offset, float],
    offsets: Sequence[Offset],
    mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    shape = padded_dgrid.shape
    labels = np.zeros(shape, dtype=np.uint32)

    if mask is None:
        mask = np.ones(shape, dtype=bool)
    else:
        mask = mask.astype(bool, copy=False)

    heap: List[Tuple[float, int, int, int, int, int]] = []
    seeds = np.argwhere((markers > 0) & mask)
    for i, j, k in seeds:
        lab = int(markers[i, j, k])
        labels[i, j, k] = np.uint32(lab)
        prio = float(padded_dgrid[i, j, k])
        heapq.heappush(
            heap,
            (-prio, -flat_index(int(i), int(j), int(k), shape), int(i), int(j), int(k), lab),
        )

    nx, ny, nz = shape
    while heap:
        neg_prio, _neg_flat, i, j, k, lab = heapq.heappop(heap)
        parent_prio = -neg_prio

        for dx, dy, dz in offsets:
            ni, nj, nk = i + dx, j + dy, k + dz
            if ni < 0 or nj < 0 or nk < 0 or ni >= nx or nj >= ny or nk >= nz:
                continue
            if not mask[ni, nj, nk]:
                continue
            if labels[ni, nj, nk] != 0:
                continue

            cand = float(padded_dgrid[ni, nj, nk]) * float(inv_w[(dx, dy, dz)])
            next_prio = min(parent_prio, cand)
            labels[ni, nj, nk] = np.uint32(lab)
            heapq.heappush(
                heap,
                (
                    -next_prio,
                    -flat_index(ni, nj, nk, shape),
                    ni,
                    nj,
                    nk,
                    lab,
                ),
            )

    return labels


def crop_padded_labels(
    labels: np.ndarray, padding: int, core_shape: Tuple[int, int, int]
) -> np.ndarray:
    p = int(padding)
    nx, ny, nz = core_shape
    return labels[p : p + nx, p : p + ny, p : p + nz].astype(np.uint32, copy=False)


def summarize(values: np.ndarray) -> Dict[str, float]:
    if values.size == 0:
        return {"count": 0.0, "min": np.nan, "median": np.nan, "max": np.nan}
    return {
        "count": float(values.size),
        "min": float(np.min(values)),
        "median": float(np.median(values)),
        "max": float(np.max(values)),
    }


def summarize_connection_bottlenecks(
    connections: Sequence,
    pocket1: Optional[set[int]] = None,
) -> Dict[str, float]:
    """Min diameter along each connection path; optional self-loop filter."""
    all_mins: List[float] = []
    self_mins: List[float] = []
    for c in connections:
        d = np.asarray(c[-1], dtype=np.float64)
        if d.size == 0:
            continue
        b = float(np.min(d))
        all_mins.append(b)
        if pocket1 is not None:
            pair = [int(x) for x in c[0]]
            if pair[0] in pocket1 and pair[1] in pocket1:
                self_mins.append(b)
    out: Dict[str, float] = {
        "count": float(len(all_mins)),
        "min": float(np.min(all_mins)) if all_mins else np.nan,
        "median": float(np.median(all_mins)) if all_mins else np.nan,
        "max": float(np.max(all_mins)) if all_mins else np.nan,
    }
    if pocket1 is not None:
        out["self_count"] = float(len(self_mins))
        out["self_max"] = float(np.max(self_mins)) if self_mins else np.nan
    return out

def boundary_values_for_label_sets(
    labels: np.ndarray,
    padded_dgrid: np.ndarray,
    mask: np.ndarray,
    offsets_half: Sequence[Offset],
    set_a: set[int],
    set_b: Optional[set[int]] = None,
) -> np.ndarray:
    nx, ny, nz = labels.shape
    out: List[float] = []

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if not mask[i, j, k]:
                    continue
                la = int(labels[i, j, k])
                if la == 0:
                    continue

                for dx, dy, dz in offsets_half:
                    ni, nj, nk = i + dx, j + dy, k + dz
                    if ni < 0 or nj < 0 or nk < 0 or ni >= nx or nj >= ny or nk >= nz:
                        continue
                    if not mask[ni, nj, nk]:
                        continue

                    lb = int(labels[ni, nj, nk])
                    if lb == 0 or lb == la:
                        continue

                    if set_b is None:
                        cond = (la in set_a) and (lb in set_a)
                    else:
                        cond = (la in set_a and lb in set_b) or (
                            la in set_b and lb in set_a
                        )

                    if cond:
                        out.append(
                            float(
                                min(
                                    padded_dgrid[i, j, k],
                                    padded_dgrid[ni, nj, nk],
                                )
                            )
                        )

    return np.asarray(out, dtype=np.float64)


def print_comparison(
    std_self: Dict[str, float],
    met_self: Dict[str, float],
    std_inter: Dict[str, float],
    met_inter: Dict[str, float],
) -> None:
    header = (
        f"{'Boundary class':<26}"
        f"{'Run':<10}"
        f"{'Count':>10}"
        f"{'Min(A)':>12}"
        f"{'Median(A)':>12}"
        f"{'Max(A)':>12}"
    )
    print("\n=== Metric-aware watershed comparison (VIGNOQ) ===")
    print(header)
    print("-" * len(header))

    rows = [
        ("self-loop (pocket-1)", "standard", std_self),
        ("self-loop (pocket-1)", "metric", met_self),
        ("inter-pocket (1<->2)", "standard", std_inter),
        ("inter-pocket (1<->2)", "metric", met_inter),
    ]
    for cls, run, s in rows:
        print(
            f"{cls:<26}{run:<10}"
            f"{int(s['count']):>10d}"
            f"{s['min']:>12.4f}"
            f"{s['median']:>12.4f}"
            f"{s['max']:>12.4f}"
        )


def trim_connections(connections: Sequence, unique_regions: np.ndarray) -> list:
    return [c for c in connections if np.any(np.isin(c[0], unique_regions))]


def write_visualization_html(
    out_path: Path,
    mat_atoms,
    dgrid: np.ndarray,
    seg: dict,
    regions_standard: np.ndarray,
    regions_metric: np.ndarray,
) -> None:
    mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
    regions_pbc = np.asarray(seg["regions_pbc"], dtype=np.uint32)
    regions_rust = np.asarray(seg["regions"], dtype=np.uint32)
    maxima = np.asarray(seg["all_maxima"])
    maxima_pbc = np.asarray(seg["maxima_pbc"])
    connections = seg.get("connections", [])
    unique_regions = np.asarray(seg["unique_regions"])
    connections_trimmed = trim_connections(connections, unique_regions)

    cscale = px.colors.qualitative.Pastel
    isovals = [1.0, 1.5]

    frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)
    isosurface_pbc = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=isovals[0],
        ase_atoms=mat_atoms,
        region_labels=regions_pbc,
        opacity=1,
        cscale=cscale,
        stylize=True,
        name="Regions PBC (rust)",
    )
    isosurface_rust = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=isovals[0],
        ase_atoms=mat_atoms,
        region_labels=regions_rust,
        opacity=0.85,
        cscale=cscale,
        stylize=True,
        name="Regions (rust)",
    )
    isosurface_standard = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=isovals[0],
        ase_atoms=mat_atoms,
        region_labels=regions_standard,
        opacity=0.85,
        cscale=cscale,
        stylize=True,
        name="Regions (proto standard)",
    )
    isosurface_metric = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=isovals[0],
        ase_atoms=mat_atoms,
        region_labels=regions_metric,
        opacity=0.85,
        cscale=cscale,
        stylize=True,
        name="Regions (proto metric)",
    )
    cell_trace = mgr.plot_cell(mat_atoms, color="orange")

    maxima_trace = mgr.plot_maxima(
        dgrid,
        maxima,
        mat_atoms,
        show_spheres=False,
        marker_size=7.5,
        marker_colorscale=cscale,
        name="Maxima",
    )
    maxima_trace_pbc = mgr.plot_maxima(
        dgrid,
        maxima_pbc,
        mat_atoms,
        show_spheres=False,
        marker_size=7.5,
        marker_colorscale=cscale,
        name="Maxima (pbc)",
    )

    connections_traces = mgr.plot_connections_as_it_is(
        connections_trimmed,
        maxima,
        regions_rust.shape,
        mat_atoms,
        cscale=cscale,
        plot_disconnected=False,
        plot_periodic=True,
    )
    if len(connections_traces) >= 2:
        connections_traces[-2].marker.size = 1

    windows_trace = mgr.plot_windows(
        connections_trimmed,
        name="Windows",
        marker_color="yellow",
        marker_size=10,
    )

    fig = go.Figure(
        data=[
            *frame_trace,
            cell_trace,
            maxima_trace_pbc,
            isosurface_pbc,
            *connections_traces,
            windows_trace,
            isosurface_rust,
            isosurface_standard,
            isosurface_metric,
            maxima_trace,
        ],
        layout=mgr.dark_layout,
    )
    fig.update_layout(scene_aspectmode="data", legend_groupclick="toggleitem")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(out_path))
    print(f"Wrote visualization: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cif",
        type=Path,
        default=VIGNOQ_CIF,
        help="Input CIF (default: VIGNOQ_clean.cif)",
    )
    parser.add_argument(
        "--html",
        type=Path,
        nargs="?",
        const=DEFAULT_HTML,
        default=None,
        help="Write Plotly HTML (default path if flag alone)",
    )
    args = parser.parse_args()

    cif_path = args.cif.expanduser().resolve()
    if not cif_path.exists():
        raise FileNotFoundError(f"Missing CIF: {cif_path}")

    atoms = read(cif_path)
    print(f"Loaded: {cif_path.name} (n_atoms={len(atoms)})")
    print(f"Cell angles (deg): {atoms.cell.angles()}")

    radii = mgr.default_vdw_radii_angstrom(atoms)
    dgrid = dgrid_from_atoms_rust(atoms, spacing=SPACING, radii=radii).astype(np.float32)
    print(f"dgrid shape: {dgrid.shape}, spacing: {SPACING} A")
    plr = float(find_pore_limiting_radius(dgrid))
    print(f"PLD (dgrid percolation): {plr:.4f} A")

    base_padding = int(calculate_padding_thickness(dgrid, atoms, SPACING))
    padding = base_padding + int(np.ceil(EXTRA_PADDING_ANGSTROM / SPACING))
    padded_dgrid = np.pad(dgrid, int(padding), mode="wrap")
    padded_void = np.pad(dgrid > MASK_THICKNESS, int(padding), mode="wrap")
    print(
        f"padding_thickness: base={base_padding}, +extra={int(np.ceil(EXTRA_PADDING_ANGSTROM / SPACING))}, total={padding}"
    )
    print(f"padded grid shape: {padded_dgrid.shape}")

    seg = do_periodic_segmentation_rust(
        dgrid=dgrid,
        mat_atoms=atoms,
        spacing=SPACING,
        mask_thickness=MASK_THICKNESS,
        h=H_PARAM,
        extra_padding=EXTRA_PADDING_ANGSTROM,
        compute_connections=bool(args.html is not None),
        return_padded_regions=False,
        return_dgrid=False,
    )
    maxima_info = np.asarray(seg["maxima_info"], dtype=int)
    pbc_groups = [set(int(x) for x in g) for g in seg["pbc_groups"]]

    if len(pbc_groups) < 2:
        raise RuntimeError(
            f"Need at least 2 pbc_groups for pocket-1/pocket-2 comparison, got {len(pbc_groups)}"
        )

    pocket1 = pbc_groups[0]
    pocket2 = pbc_groups[1]
    print(
        f"Using pockets: pocket-1 size={len(pocket1)} labels, pocket-2 size={len(pocket2)} labels"
    )

    markers = build_unmerged_markers_from_maxima_info(
        maxima_info=maxima_info,
        dgrid_shape=dgrid.shape,
        padding=padding,
    )
    print(f"unmerged marker seeds in padded grid: {int(np.count_nonzero(markers))}")

    offsets = neighbor_offsets_26()
    offsets_half = positive_half_offsets(offsets)

    inv_w_standard = {o: 1.0 for o in offsets}
    cell_matrix = np.asarray(atoms.cell.array, dtype=np.float64)
    inv_w_metric, weight_stats = compute_metric_inv_weights(
        cell_matrix=cell_matrix,
        grid_shape=dgrid.shape,
        spacing=SPACING,
        offsets=offsets,
    )
    print(
        "metric weights: "
        f"min={weight_stats['w_min']:.4f}, max={weight_stats['w_max']:.4f}, "
        f"ratio={weight_stats['anisotropy_ratio']:.4f}"
    )

    labels_standard = metric_aware_watershed(
        padded_dgrid=padded_dgrid,
        markers=markers,
        inv_w=inv_w_standard,
        offsets=offsets,
        mask=padded_void,
    )
    labels_metric = metric_aware_watershed(
        padded_dgrid=padded_dgrid,
        markers=markers,
        inv_w=inv_w_metric,
        offsets=offsets,
        mask=padded_void,
    )

    regions_standard = crop_padded_labels(labels_standard, padding, dgrid.shape)
    regions_metric = crop_padded_labels(labels_metric, padding, dgrid.shape)
    print(
        "cropped region label counts: "
        f"standard={len(np.unique(regions_standard)) - 1}, "
        f"metric={len(np.unique(regions_metric)) - 1}, "
        f"rust={len(np.unique(seg['regions'])) - 1}"
    )

    std_self = summarize(
        boundary_values_for_label_sets(
            labels=labels_standard,
            padded_dgrid=padded_dgrid,
            mask=padded_void,
            offsets_half=offsets_half,
            set_a=pocket1,
            set_b=None,
        )
    )
    met_self = summarize(
        boundary_values_for_label_sets(
            labels=labels_metric,
            padded_dgrid=padded_dgrid,
            mask=padded_void,
            offsets_half=offsets_half,
            set_a=pocket1,
            set_b=None,
        )
    )
    std_inter = summarize(
        boundary_values_for_label_sets(
            labels=labels_standard,
            padded_dgrid=padded_dgrid,
            mask=padded_void,
            offsets_half=offsets_half,
            set_a=pocket1,
            set_b=pocket2,
        )
    )
    met_inter = summarize(
        boundary_values_for_label_sets(
            labels=labels_metric,
            padded_dgrid=padded_dgrid,
            mask=padded_void,
            offsets_half=offsets_half,
            set_a=pocket1,
            set_b=pocket2,
        )
    )

    print_comparison(std_self, met_self, std_inter, met_inter)

    seg_metric = do_periodic_segmentation_rust(
        dgrid=dgrid,
        mat_atoms=atoms,
        spacing=SPACING,
        mask_thickness=MASK_THICKNESS,
        h=H_PARAM,
        extra_padding=EXTRA_PADDING_ANGSTROM,
        compute_connections=False,
        return_padded_regions=True,
        watershed_backend="rust_det_metric",
    )
    rust_metric_regions = np.asarray(seg_metric["regions"], dtype=np.uint32)
    diff = int(np.sum(regions_metric != rust_metric_regions))
    print(
        f"rust_det_metric full pipeline vs standalone proto flood: {diff} voxel diff "
        "(expected: merged metric flood + cleanup differs from proto-only flood)"
    )
    rust_self = summarize(
        boundary_values_for_label_sets(
            labels=np.asarray(seg_metric["regions_padded"]),
            padded_dgrid=padded_dgrid,
            mask=padded_void,
            offsets_half=offsets_half,
            set_a=pocket1,
            set_b=None,
        )
    )
    if abs(rust_self["max"] - met_self["max"]) > 0.05:
        raise RuntimeError(
            f"rust_det_metric self-loop boundary max {rust_self['max']:.4f} A "
            f"differs from prototype {met_self['max']:.4f} A"
        )
    print("rust_det_metric: PASS (functional parity with prototype metric flood)")

    print("\nExpected targets:")
    print("- metric self-loop max near ~1.38 A (not ~4.14 A)")
    print("- inter-pocket boundaries near ~1.48 A in both runs")

    seg_metric_viz = do_periodic_segmentation_rust(
        dgrid=dgrid,
        mat_atoms=atoms,
        spacing=SPACING,
        mask_thickness=MASK_THICKNESS,
        h=H_PARAM,
        extra_padding=EXTRA_PADDING_ANGSTROM,
        compute_connections=True,
        watershed_backend="rust_det_metric",
    )
    conn_stats = summarize_connection_bottlenecks(
        seg_metric_viz.get("connections", []), pocket1=pocket1
    )
    print(
        "connection bottlenecks (metric seg): "
        f"max={conn_stats['max']:.4f} A, "
        f"self-loop max={conn_stats.get('self_max', float('nan')):.4f} A"
    )

    if args.html is not None:
        write_visualization_html(
            out_path=args.html,
            mat_atoms=atoms,
            dgrid=dgrid,
            seg=seg_metric_viz,
            regions_standard=regions_standard,
            regions_metric=regions_metric,
        )


if __name__ == "__main__":
    main()

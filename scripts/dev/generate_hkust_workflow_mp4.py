#!/usr/bin/env python3
"""
Generate a staged MP4 workflow animation for HKUST-1 (Cu-BTC).

Stages:
1) Framework + unit cell (camera orbits around cell center)
2) Local maxima
3) Segmentation surfaces: raw regions (opaque) -> PBC merge (opaque) -> detected pockets / clusters (opaque)
4) Transparent cavity surface + windows
5) Pore graph on faint transparent cavity surface (framework + cell for context)
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from ase.geometry import complete_cell
from ase.io import read
from skimage.measure import regionprops


def _repo_root() -> Path:
    return Path(__file__).resolve().parent


def _resolve_cif(repo_root: Path) -> Path:
    candidates = [
        repo_root / "mofography" / "examples" / "TestCIFs" / "P1-Cu-BTC.cif",
        repo_root / "examples" / "TestCIFs" / "P1-Cu-BTC.cif",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError("Could not find P1-Cu-BTC.cif in expected example paths.")


def _cell_center_cart(ase_atoms) -> np.ndarray:
    """Cartesian coordinates of fractional (0.5, 0.5, 0.5) — unit-cell body center."""
    cell = complete_cell(ase_atoms.get_cell()).T  # columns = lattice vectors (Å)
    frac = np.array([0.5, 0.5, 0.5], dtype=float)
    return np.dot(cell, frac)


def _camera_from_orbit(
    theta: float,
    center_cart: np.ndarray,
    *,
    x_bias: float,
    orbit_r: float,
) -> dict:
    """View roughly along +x toward ``center_cart``; orbit in the yz plane; z is up."""
    ex = float(center_cart[0] + x_bias)
    ey = float(center_cart[1] + orbit_r * np.cos(theta))
    ez = float(center_cart[2] + orbit_r * np.sin(theta))
    return dict(
        eye=dict(x=ex, y=ey, z=ez),
        center=dict(
            x=float(center_cart[0]),
            y=float(center_cart[1]),
            z=float(center_cart[2]),
        ),
        up=dict(x=0, y=0, z=1),
    )


def prepare_hkust_workflow_export(repo_root: Path | None = None) -> dict:
    """Build the full HKUST workflow figure and export metadata (used by MP4 main + benchmarks)."""
    if repo_root is None:
        repo_root = _repo_root()
    sys.path.insert(0, str(repo_root))
    import mofography as mgr

    cif_path = _resolve_cif(repo_root)
    output_dir = repo_root / "outputs" / "hkust_workflow_animation"
    frames_dir = output_dir / "frames"
    output_dir.mkdir(parents=True, exist_ok=True)
    frames_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/5] Loading structure: {cif_path}")
    mat_atoms = read(str(cif_path))
    mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

    print("[2/5] Running periodic segmentation")
    spacing = 0.5
    dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=spacing)
    results = mgr.do_periodic_segmentation_rust(
        dgrid,
        mat_atoms,
        spacing=spacing,
        mask_thickness=1.0,
        h=0.2,
        threshold_abs=1.0,
        min_distance=1.5,
        return_dgrid=True,
    )

    dgrid = results["dgrid"]
    regions = results["regions"]
    regions_pbc = results["regions_pbc"]
    maxima_info = results["maxima_info"]
    pbc_groups = results["pbc_groups"]
    connections = results["connections"]
    all_maxima = results["all_maxima"]
    unique_regions = results["unique_regions"]

    maxima = all_maxima[: regions.max()]
    all_maxima_radii = dgrid[tuple(maxima_info[:, 1:4].T)]
    maxima_radii = all_maxima_radii[: regions.max()]

    connections_trimmed = [
        c for c in connections if np.all(np.isin(c[0], unique_regions))
    ]
    if len(connections_trimmed) == 0:
        connections_trimmed = connections

    print("[3/5] Building clustering and workflow traces")
    rag = mgr.rag_from_connections(connections_trimmed, maxima)
    rag = mgr.add_maxima_to_rag(rag, maxima, maxima_radii, dgrid.shape, mat_atoms)
    rag = mgr.add_pbc_to_rag(rag, pbc_groups)
    props = regionprops(regions)
    bins = np.linspace(0, 50, 30)
    rag = mgr.add_psd_to_rag_with_pbc(rag, props, pbc_groups, bins, mat_atoms, dgrid.shape)
    rag = mgr.add_clustering_to_rag_geometry_js(rag, props, pbc_groups, mat_atoms, dgrid.shape)
    regions_cluster = mgr.clustering_grid_from_regions(
        regions, rag, cluster_key="cluster_geom_pbc"
    )

    center_cart = _cell_center_cart(mat_updated)
    pos = mat_updated.get_positions()
    span = float(np.max(pos.max(axis=0) - pos.min(axis=0)))
    span = max(span, 1.0)
    # Zoomed-in camera: modest offset along +x, tight yz orbit (same Å units as structure).
    x_bias = 0.52 * span
    orbit_r = 0.36 * span

    cscale = px.colors.qualitative.Pastel
    isoval = 2.0

    frame_traces = mgr.plot_framework_as_mesh(
        mat_updated, show_rings=False, stylize=True
    )
    cell_trace = mgr.plot_cell(mat_updated, color="black", name="unit cell")
    maxima_trace = mgr.plot_maxima(
        dgrid,
        maxima,
        mat_updated,
        name="local maxima",
        marker_size=9,
        marker_colorscale="Viridis",
    )

    # Opaque segmentation surfaces while we introduce each labeling (pockets read clearly).
    surface_regions = mgr.plot_surface_dgrid(
        dgrid,
        isoval,
        mat_updated,
        region_labels=regions,
        select="all",
        cscale=cscale,
        stylize=True,
        opacity=1.0,
        name="segmented regions",
        hover_label="Region : ",
    )
    surface_pbc = mgr.plot_surface_dgrid(
        dgrid,
        isoval,
        mat_updated,
        region_labels=regions_pbc - 1,
        select="all",
        cscale=cscale,
        stylize=True,
        opacity=1.0,
        name="PBC-merged regions",
        hover_label="PBC group : ",
    )
    surface_cluster = mgr.plot_surface_dgrid(
        dgrid,
        isoval,
        mat_updated,
        region_labels=regions_cluster - 1,
        select="all",
        cscale=cscale,
        stylize=True,
        opacity=1.0,
        name="detected pockets (clusters)",
        hover_label="Cavity type : ",
    )

    surface_cluster_dim_windows = copy.deepcopy(surface_cluster)
    surface_cluster_dim_windows.update(opacity=0.28, name="cavity surface (for windows)")

    surface_cluster_dim_graph = copy.deepcopy(surface_cluster)
    surface_cluster_dim_graph.update(opacity=0.14, name="cavity surface (for pore graph)")

    windows_trace = mgr.plot_windows(
        connections_trimmed,
        name="windows",
        marker_color="purple",
        marker_size=8,
        opacity=0.95,
    )

    connections_traces = mgr.plot_connections_as_it_is(
        connections_trimmed,
        all_maxima,
        dgrid.shape,
        mat_updated,
        plot_periodic=True,
        plot_disconnected=False,
        line_color="crimson",
        marker_size=7,
        dot_size=1.0,
    )

    all_traces = [
        *frame_traces,
        cell_trace,
        maxima_trace,
        surface_regions,
        surface_pbc,
        surface_cluster,
        surface_cluster_dim_windows,
        surface_cluster_dim_graph,
        windows_trace,
        *connections_traces,
    ]

    trace_ids: dict[str, list[int] | int] = {}
    idx = 0
    trace_ids["frame"] = list(range(idx, idx + len(frame_traces)))
    idx += len(frame_traces)
    trace_ids["cell"] = idx
    idx += 1
    trace_ids["maxima"] = idx
    idx += 1
    trace_ids["surface_regions"] = idx
    idx += 1
    trace_ids["surface_pbc"] = idx
    idx += 1
    trace_ids["surface_cluster"] = idx
    idx += 1
    trace_ids["surface_cluster_dim_windows"] = idx
    idx += 1
    trace_ids["surface_cluster_dim_graph"] = idx
    idx += 1
    trace_ids["windows"] = idx
    idx += 1
    trace_ids["connections"] = list(range(idx, len(all_traces)))

    def vis_for(keys: list[str]) -> list[bool]:
        visibility = [False] * len(all_traces)
        for key in keys:
            value = trace_ids[key]
            if isinstance(value, list):
                for j in value:
                    visibility[j] = True
            else:
                visibility[value] = True
        return visibility

    # Preview: 720p, low fps, fewer PNGs so Kaleido finishes quickly. Set True for a fast draft.
    preview = False
    if preview:
        render_w, render_h = 1280, 720
        fps = 4
        title_font_size = 22
        n_fw = 30
        n_mx = 30
        n_surf = 30
        n_win = 24
        n_graph = 26
    else:
        # 16:9 from width 1080 — much faster than 1920×1080, still usable for review/talks.
        render_w, render_h = 1080, int(round(1080 * 9 / 16))
        fps = 6
        title_font_size = 22
        n_fw = 52
        n_mx = 52
        n_surf = 52
        n_win = 40
        n_graph = 44

    stages: list[tuple[str, list[bool], int]] = [
        ("Stage 1: Framework + unit cell", vis_for(["frame", "cell"]), n_fw),
        ("Stage 2: Local maxima (cage centers)", vis_for(["frame", "cell", "maxima"]), n_mx),
        (
            "Stage 3a: Segmented void surface (raw regions)",
            vis_for(["frame", "cell", "maxima", "surface_regions"]),
            n_surf,
        ),
        (
            "Stage 3b: Same surface with periodic boundaries merged",
            vis_for(["frame", "cell", "maxima", "surface_pbc"]),
            n_surf,
        ),
        (
            "Stage 3c: Detected pockets (geometry + PBC clusters)",
            vis_for(["frame", "cell", "maxima", "surface_cluster"]),
            n_surf,
        ),
        (
            "Stage 4: Windows on semi-transparent cavity surface",
            vis_for(
                [
                    "frame",
                    "cell",
                    "maxima",
                    "surface_cluster_dim_windows",
                    "windows",
                ]
            ),
            n_win,
        ),
        (
            "Stage 5: Pore graph over faint cavity surface",
            vis_for(
                [
                    "frame",
                    "cell",
                    "maxima",
                    "surface_cluster_dim_graph",
                    "connections",
                ]
            ),
            n_graph,
        ),
    ]

    fig = go.Figure(data=all_traces, layout=mgr.white_layout)
    fig.update_layout(
        width=render_w,
        height=render_h,
        scene=dict(aspectmode="data"),
        legend=dict(
            x=0.99,
            y=0.99,
            xanchor="right",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.75)",
        ),
    )

    total_frames = int(sum(s[2] for s in stages))

    def frame_hook_per_stage(
        f: go.Figure, stage_index: int, frame_in_stage: int, n_frames: int
    ) -> None:
        # Stages 0–4: two full rotations; later stages: 1.5 rotations (still readable).
        turns = 2.0 if stage_index <= 4 else 1.5
        denom = max(n_frames - 1, 1)
        frac = frame_in_stage / denom
        theta = frac * turns * 2.0 * np.pi
        f.update_layout(
            scene_camera=_camera_from_orbit(
                theta, center_cart, x_bias=x_bias, orbit_r=orbit_r
            )
        )

    export_stages = [
        (f"HKUST-1 (Cu-BTC) workflow — {title}", vis, n)
        for title, vis, n in stages
    ]

    output_mp4 = output_dir / "hkust_workflow_staged.mp4"
    return {
        "fig": fig,
        "export_stages": export_stages,
        "fps": fps,
        "render_w": render_w,
        "render_h": render_h,
        "title_font_size": title_font_size,
        "total_frames": total_frames,
        "frames_dir": frames_dir,
        "output_mp4": output_mp4,
        "frame_hook_per_stage": frame_hook_per_stage,
        "center_cart": center_cart,
        "x_bias": x_bias,
        "orbit_r": orbit_r,
    }


def main() -> None:
    repo_root = _repo_root()
    ctx = prepare_hkust_workflow_export(repo_root)
    fig = ctx["fig"]
    total_frames = ctx["total_frames"]
    frames_dir = ctx["frames_dir"]
    output_mp4 = ctx["output_mp4"]

    from mofography.export_mp4 import export_staged_figure_to_mp4

    print(
        f"[4/5] Rendering {total_frames} Kaleido frames to {frames_dir} (tqdm progress bar)…"
    )
    export_staged_figure_to_mp4(
        fig,
        output_mp4,
        stages=ctx["export_stages"],
        fps=ctx["fps"],
        width=ctx["render_w"],
        height=ctx["render_h"],
        scale=1,
        frame_hook_per_stage=ctx["frame_hook_per_stage"],
        title_font_size=ctx["title_font_size"],
        frames_dir=frames_dir,
        cleanup_frames=False,
    )

    print(f"[5/5] Done. MP4 written to: {output_mp4}")


if __name__ == "__main__":
    main()

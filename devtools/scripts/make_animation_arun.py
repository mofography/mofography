#!/usr/bin/env python3
"""
Build the Plotly figure from ``devtools/scripts/make_animation_arun.org`` and export
PNG frames (Kaleido) + H.264 MP4 (ffmpeg).

Trace order in the figure (stable indices for visibility)::

    [*frame_trace, cell, maxima_distinct, regions, regions_pbc, clusters,
     *connections_traces, windows]

Stages are **explicit** (not a simple prefix cumulative list): framework + cell
from the start; then maxima; then region surfaces in order; a beat with regions
+ windows (PBC/cluster surfaces off); then internal pore-graph traces; then
periodic pore-graph traces. Camera orbit runs at **one quarter** the previous
angular speed per stage (see ``_orbit_hook_factory``).

Requirements
------------
- ``kaleido`` (Plotly static image export)
- ``ffmpeg`` on PATH (to stitch PNGs into MP4)

Kaleido v1+ persistent rendering: ``mofography.export_mp4.export_staged_figure_to_mp4``
wraps all ``write_image`` calls in ``kaleido.start_sync_server()`` / ``stop_sync_server()``
by default (Plotly/Kaleido issue #400). Disable with ``--no-kaleido-sync`` or tune workers
with ``--kaleido-workers N`` when supported.

Examples
--------
  python devtools/scripts/make_animation_arun.py \\
      --cif examples/TestCIFs/P1-Cu-BTC.cif \\
      --out outputs/tutorial_reveal.mp4

  python devtools/scripts/make_animation_arun.py --keep-frames --frames-dir outputs/frames_tutorial
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from ase.io import read
from skimage.measure import regionprops

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import mofography as mgr
from mofography.export_mp4 import export_staged_figure_to_mp4
from mofography.visualize import make_bg_image


def build_tutorial_figure(cif_path: Path) -> tuple[go.Figure, list[tuple[str, list[bool]]]]:
    """Return figure and list of (stage_title, visibility_per_trace)."""
    mat_atoms = read(str(cif_path))

    spacing = 0.25
    h = 0.2
    mask_thickness = 0.5
    threshold_abs = 1.5
    extra_padding = 5
    min_distance = 1.5

    dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=20000).astype(
        np.float32
    )
    results = mgr.do_periodic_segmentation_rust(
        dgrid,
        mat_atoms,
        mask_thickness=mask_thickness,
        h=h,
        spacing=spacing,
        return_dgrid=True,
        threshold_abs=threshold_abs,
        min_distance=min_distance,
        return_padded_regions=True,
        extra_padding=extra_padding,
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
    maxima_radii = dgrid[tuple(maxima_info[:, 1:4].T)][: regions.max()]
    connections_trimmed = [
        c for c in connections if np.any(np.isin(c[0], unique_regions))
    ]

    rag = mgr.rag_from_connections(connections_trimmed, maxima)
    rag = mgr.add_maxima_to_rag(
        rag, maxima, maxima_radii, dgrid.shape, mat_atoms
    )
    rag = mgr.add_pbc_to_rag(rag, pbc_groups)
    props = regionprops(regions)
    bins = np.linspace(0, 50, 30)
    rag = mgr.add_psd_to_rag_with_pbc(
        rag, props, pbc_groups, bins, mat_atoms, dgrid.shape
    )
    rag = mgr.add_clustering_to_rag_geometry_js(
        rag, props, pbc_groups, mat_atoms, dgrid.shape
    )
    regions_cluster = mgr.clustering_grid_from_regions(
        regions, rag, cluster_key="cluster_geom_pbc"
    )

    cscale = px.colors.qualitative.Pastel
    isovals = [1.0, 1.5]
    mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

    isosurface_segmentation = mgr.plot_surface_dgrid(
        dgrid,
        stylize=True,
        isovalue=isovals[0],
        region_labels=regions,
        ase_atoms=mat_atoms,
        name="Regions",
        opacity=1,
        cscale=cscale,
    )
    isosurface_pbc = mgr.plot_surface_dgrid(
        dgrid,
        stylize=True,
        isovalue=isovals[0],
        region_labels=regions_pbc,
        ase_atoms=mat_atoms,
        name="Regions PBC",
        opacity=1,
        cscale=cscale,
    )
    cell_trace = mgr.plot_cell(mat_atoms, color="orange")
    maxima_trace = mgr.plot_maxima(
        dgrid,
        maxima,
        mat_atoms,
        marker_size=10,
        show_spheres=False,
        marker_colorscale=cscale,
        name="Maxima (distinct)",
    )
    frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)
    isosurface_clusters = mgr.plot_surface_dgrid(
        dgrid,
        stylize=True,
        isovalue=isovals[0],
        region_labels=regions_cluster,
        name="Pocket Types by geometry",
        ase_atoms=mat_atoms,
        opacity=1,
        cscale=cscale,
    )
    connections_traces = mgr.plot_connections_as_it_is(
        connections_trimmed,
        all_maxima,
        regions.shape,
        mat_atoms,
        cscale=cscale,
        plot_disconnected=False,
        plot_periodic=True,
    )
    # Windows only from internal connections (exclude periodic); full graph still uses connections_trimmed.
    connections_for_windows = [c for c in connections_trimmed if c[2] == "internal"]
    windows_trace = mgr.plot_windows(
        connections_for_windows,
        name="Windows",
        marker_color="yellow",
        marker_size=10,
        marker_opacity=0.8,
    )

    # Trace order chosen for animation: framework + cell + maxima first, then
    # surfaces, then pore graph slices, then windows.
    traces = [
        *frame_trace,
        cell_trace,
        maxima_trace,
        isosurface_segmentation,
        isosurface_pbc,
        isosurface_clusters,
        *connections_traces,
        windows_trace,
    ]
    fig = go.Figure(traces, layout=mgr.white_layout)
    # Legend grows with trace count and can shrink the 3D viewport; reserve width and
    # use a smaller legend so the scene (and unit cell framing) stays more stable.
    fig.update_layout(
        margin=dict(t=20, b=20, l=20, r=220),
        legend=dict(
            x=1.02,
            y=0.98,
            xanchor="left",
            yanchor="top",
            font=dict(size=9),
            groupclick="togglegroup",
            itemsizing="constant",
        ),
    )

    n = len(fig.data)
    n_frame = len(frame_trace)
    n_conn = len(connections_traces)

    def _all_false() -> list[bool]:
        return [False] * n

    def _idx_framework() -> list[int]:
        return list(range(n_frame))

    i_cell = n_frame
    i_max = n_frame + 1
    i_regions = n_frame + 2
    i_pbc = n_frame + 3
    i_cluster = n_frame + 4
    i_conn0 = n_frame + 5
    i_windows = n_frame + 5 + n_conn

    if i_windows >= n:
        raise RuntimeError("Figure trace index mismatch (windows).")

    def _conn_idx(*substrings: str) -> list[int]:
        out: list[int] = []
        for j in range(n_conn):
            tr = fig.data[i_conn0 + j]
            nm = getattr(tr, "name", "") or ""
            if any(s in nm for s in substrings):
                out.append(i_conn0 + j)
        return out

    idx_internal_nodes = _conn_idx("Internal Maxima")
    idx_internal_edges = _conn_idx("Internal Connections")
    idx_periodic = _conn_idx("Periodic Connections", "Periodic Maxima")

    stages: list[tuple[str, list[bool]]] = []

    def stage(title: str, true_idx: set[int]) -> None:
        vis = [i in true_idx for i in range(n)]
        stages.append((title, vis))

    fw = set(_idx_framework())

    # 0 — Framework + unit cell (video start)
    s0 = fw | {i_cell}
    stage("Framework + unit cell", s0)

    # 1 — + distinct maxima (pocket centers)
    stage("+ Maxima (distinct)", s0 | {i_max})

    # 2 — + regions surface
    stage("+ Regions surface", s0 | {i_max, i_regions})

    # 3 — + regions PBC surface
    stage("+ Regions PBC surface", s0 | {i_max, i_regions, i_pbc})

    # 4 — + pocket types (cluster) surface
    stage("+ Pocket types (cluster) surface", s0 | {i_max, i_regions, i_pbc, i_cluster})

    # 5 — Back to regions + windows (turn PBC + cluster surfaces off)
    stage("Regions + windows", s0 | {i_max, i_regions, i_windows})

    # 6 — + internal pore graph (nodes + edges)
    s6 = set(s0) | {i_max, i_regions, i_windows}
    for j in idx_internal_nodes + idx_internal_edges:
        s6.add(j)
    stage("+ Pore graph (internal)", s6)

    # 7 — + periodic pore graph (only if those traces exist)
    if idx_periodic:
        s7 = set(s6)
        for j in idx_periodic:
            s7.add(j)
        stage("+ Pore graph (periodic)", s7)

    return fig, stages


def _camera_dict(
    eye: tuple[float, float, float],
    center: tuple[float, float, float],
    up: tuple[float, float, float],
) -> dict:
    return {
        "eye": {"x": eye[0], "y": eye[1], "z": eye[2]},
        "center": {"x": center[0], "y": center[1], "z": center[2]},
        "up": {"x": up[0], "y": up[1], "z": up[2]},
    }


def _orbit_hook_factory(
    base_eye: tuple[float, float, float],
    center: tuple[float, float, float],
    up: tuple[float, float, float],
):
    """Eye-orbit in scene space; **quarter** angular speed vs original (π/2 per stage sweep)."""

    def hook(
        fig: go.Figure,
        stage_index: int,
        frame_in_stage: int,
        n_frames_in_stage: int,
    ) -> None:
        t = frame_in_stage / max(n_frames_in_stage, 1)
        # Previously one full turn per stage (2π·t). Slow to ¼: (π/2)·t per stage.
        theta0 = math.atan2(base_eye[1], base_eye[0])
        theta = theta0 + (math.pi / 2.0) * t
        r = math.hypot(base_eye[0], base_eye[1])
        fig.update_layout(
            scene_camera=_camera_dict(
                eye=(r * math.cos(theta), r * math.sin(theta), base_eye[2]),
                center=center,
                up=up,
            )
        )

    return hook


def main() -> None:
    p = argparse.ArgumentParser(description="Tutorial Plotly figure → PNG frames → MP4.")
    p.add_argument(
        "--cif",
        type=Path,
        default=Path("examples/TestCIFs/P1-Cu-BTC.cif"),
        help="Input CIF (same default as org tutorial).",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("outputs/tutorial_reveal_plotly.mp4"),
        help="Output MP4 path.",
    )
    p.add_argument("--fps", type=int, default=12)
    p.add_argument("--width", type=int, default=1280)
    p.add_argument("--height", type=int, default=720)
    p.add_argument(
        "--frames-per-stage",
        type=int,
        default=24,
        help="PNG count per cumulative stage (orbit samples this interval).",
    )
    p.add_argument("--crf", type=int, default=20, help="libx264 CRF (lower = better).")
    p.add_argument(
        "--keep-frames",
        action="store_true",
        help="Keep PNG directory (implies --frames-dir if not set).",
    )
    p.add_argument(
        "--frames-dir",
        type=Path,
        default=None,
        help="Directory for frame_%04d.png (default: temp unless --keep-frames).",
    )
    p.add_argument(
        "--no-orbit",
        action="store_true",
        help="Disable per-frame camera motion.",
    )
    p.add_argument(
        "--camera-eye",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=(2.0, 0.0, 1.35),
        help="Camera eye vector (default: 2.0 0.0 1.35).",
    )
    p.add_argument(
        "--camera-center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=(0.0, 0.0, 0.0),
        help="Camera center vector (default: 0.0 0.0 0.0).",
    )
    p.add_argument(
        "--camera-up",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=(0.0, 0.0, 1.0),
        help="Camera up vector (default: 0.0 0.0 1.0).",
    )
    p.add_argument(
        "--bg-type",
        type=str,
        default=None,
        choices=("radial", "left_to_right", "bottom_to_top"),
        help="Optional background gradient type (uses mofography.visualize.make_bg_image).",
    )
    p.add_argument(
        "--bg-cmap",
        type=str,
        default="plasma_r",
        help="Colormap name for --bg-type (default: plasma_r).",
    )
    p.add_argument(
        "--no-kaleido-sync",
        action="store_true",
        help="Do not use kaleido.start_sync_server (slower: new engine per frame on Kaleido 1.x).",
    )
    p.add_argument(
        "--kaleido-workers",
        type=int,
        default=None,
        metavar="N",
        help="Optional parallel workers for kaleido.start_sync_server(n=N) when supported.",
    )
    args = p.parse_args()

    cif = args.cif.expanduser().resolve()
    if not cif.is_file():
        raise SystemExit(f"CIF not found: {cif}")

    fig, stage_vis_specs = build_tutorial_figure(cif)
    camera_eye = tuple(args.camera_eye)
    camera_center = tuple(args.camera_center)
    camera_up = tuple(args.camera_up)
    fig.update_layout(scene_camera=_camera_dict(camera_eye, camera_center, camera_up))
    if args.bg_type is not None:
        fig.update_layout(images=[make_bg_image(type=args.bg_type, cmap=args.bg_cmap)])
    stages = [
        (title, vis, args.frames_per_stage) for title, vis in stage_vis_specs
    ]

    frames_dir = args.frames_dir
    if args.keep_frames and frames_dir is None:
        frames_dir = args.out.parent / (args.out.stem + "_frames")

    hook = (
        None
        if args.no_orbit
        else _orbit_hook_factory(
            base_eye=camera_eye,
            center=camera_center,
            up=camera_up,
        )
    )

    out = export_staged_figure_to_mp4(
        fig,
        args.out,
        stages=stages,
        fps=args.fps,
        width=args.width,
        height=args.height,
        frame_hook_per_stage=hook,
        crf=args.crf,
        frames_dir=frames_dir,
        cleanup_frames=not args.keep_frames,
        show_progress=True,
        kaleido_sync_server=not args.no_kaleido_sync,
        kaleido_sync_workers=args.kaleido_workers,
    )
    print(f"Wrote {out.resolve()}")
    if frames_dir is not None:
        print(f"Frames directory: {frames_dir.resolve()}")


if __name__ == "__main__":
    main()

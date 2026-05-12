#!/usr/bin/env python3
"""
Workflow animation HTML for mofography - FIXED VERSION.

Stages:
1) Framework + unit cell (legacy centered z-spin camera)
2) Add local maxima
3a) Opaque segmented surface (raw regions)
3b) Opaque segmented surface (regions_pbc)
3c) Opaque segmented surface (regions_cluster)
4) Transparent surface + windows
5) Transparent surface + pore graph (internal)
6) Transparent surface + pore graph (periodic emphasized)
7) Final complete picture (everything opaque)

Key fixes:
- Camera uses legacy centered z-spin: eye rotates in xy plane, center=(0,0,0), z fixed
- Surfaces use plot_surface_dgrid directly (correct 3D Mesh3d)
- All traces are 3D (Mesh3d, Scatter3d) - no 2D traces
"""

from __future__ import annotations

import argparse
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
    return Path(__file__).resolve().parent.parent.parent


def _camera_legacy_z_spin(theta: float, radius: float = 2.0, z: float = 1.2) -> dict:
    """Legacy centered camera style: rotate eye in xy plane, keep z fixed/up.
    
    This is the "original centered look" from cif_camera_animation_html.py:
    - Eye rotates in xy plane, fixed z
    - center = (0, 0, 0)
    - up = (0, 0, 1)
    - No "planet orbit" around off-origin center
    """
    return dict(
        eye=dict(
            x=float(radius * np.cos(theta)),
            y=float(radius * np.sin(theta)),
            z=float(z),
        ),
        center=dict(x=0.0, y=0.0, z=0.0),
        up=dict(x=0, y=0, z=1),
    )


def build_workflow_figure(
    atoms,
    *,
    spacing: float = 0.25,
    frames: int = 120,
    fps: float = 12.0,
    frames_per_rotation: float = 120.0,
    camera_eye: tuple[float, float, float] = (2.0, 0.0, 1.2),
    camera_up: tuple[float, float, float] = (0.0, 0.0, 1.0),
    camera_center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    redraw: bool = False,
) -> go.Figure:
    import mofography as mgr

    atoms_updated = mgr.update_with_covalent_bonds(atoms)

    # Trace recipe follows devtools/scripts/make_animation_arun.org
    # ── Segmentation ────────────────────────────────────────────────────────────
    dgrid = mgr.dgrid_from_atoms_cpu(atoms, spacing=spacing, block_size=20000).astype(np.float32)
    seg = mgr.do_periodic_segmentation_rust(
        dgrid,
        atoms,
        spacing=spacing,
        mask_thickness=0.5,
        h=0.2,
        threshold_abs=1.5,
        min_distance=1.5,
        return_dgrid=True,
        return_padded_regions=True,
        extra_padding=5,
    )

    dgrid = seg["dgrid"]
    regions = seg["regions"]
    regions_pbc = seg["regions_pbc"]
    maxima_info = seg["maxima_info"]
    pbc_groups = seg["pbc_groups"]
    connections = seg["connections"]
    all_maxima = seg["all_maxima"]
    unique_regions = seg["unique_regions"]
    maxima_pbc = seg["maxima_pbc"]

    maxima = all_maxima[: regions.max()]
    all_maxima_radii = dgrid[tuple(maxima_info[:, 1:4].T)]
    maxima_radii = all_maxima_radii[: regions.max()]

    connections_trimmed = [
        c for c in connections if np.any(np.isin(c[0], unique_regions))
    ]
    if not connections_trimmed:
        connections_trimmed = connections

    rag = mgr.rag_from_connections(connections_trimmed, maxima)
    rag = mgr.add_maxima_to_rag(rag, maxima, maxima_radii, dgrid.shape, atoms)
    rag = mgr.add_pbc_to_rag(rag, pbc_groups)
    props = regionprops(regions)
    bins = np.linspace(0, 50, 30)
    rag = mgr.add_psd_to_rag_with_pbc(rag, props, pbc_groups, bins, atoms, dgrid.shape)
    rag = mgr.add_clustering_to_rag_geometry_js(rag, props, pbc_groups, atoms, dgrid.shape)
    regions_cluster = mgr.clustering_grid_from_regions(
        regions, rag, cluster_key="cluster_geom_pbc"
    )

    cscale = px.colors.qualitative.Pastel
    isoval = 1.0

    # ── Traces (based on make_animation_arun.org) ──────────────────────────────

    frame_traces = mgr.plot_framework_as_mesh(
        atoms_updated, show_rings=False, stylize=True
    )
    cell_trace = mgr.plot_cell(atoms, color="orange", name="unit cell")
    maxima_trace = mgr.plot_maxima(
        dgrid,
        maxima,
        atoms,
        marker_size=10,
        show_spheres=False,
        marker_colorscale=cscale,
        name="Maxima (distinct)",
    )
    surf_regions_trace = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=isoval,
        ase_atoms=atoms,
        region_labels=regions,
        select="all",
        opacity=1.0,
        cscale=cscale,
        stylize=True,
        name="Regions",
    )
    surf_pbc_trace = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=isoval,
        ase_atoms=atoms,
        region_labels=regions_pbc,
        select="all",
        opacity=1.0,
        cscale=cscale,
        stylize=True,
        name="Regions PBC",
    )
    surf_cluster_trace = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=isoval,
        ase_atoms=atoms,
        region_labels=regions_cluster,
        select="all",
        opacity=1.0,
        cscale=cscale,
        stylize=True,
        name="Pocket Types by geometry",
    )

    # dimmed version for manual window/graph overlay
    surf_cluster_trans = copy.deepcopy(surf_cluster_trace)
    surf_cluster_trans.update(opacity=0.22, name="Pocket Types (transparent)")

    windows_trace = mgr.plot_windows(
        connections_trimmed,
        name="Windows",
        marker_color="yellow",
        marker_size=10,
        opacity=0.8,
    )

    graph_periodic = mgr.plot_connections_as_it_is(
        connections_trimmed, all_maxima, dgrid.shape, atoms_updated,
        cscale=cscale, plot_disconnected=False, plot_periodic=True,
    )

    # Optional pbc maxima trace included but initially hidden
    maxima_trace_pbc = mgr.plot_maxima(
        dgrid,
        maxima_pbc,
        atoms,
        show_spheres=False,
        marker_size=7.5,
        marker_colorscale=cscale,
        name="Maxima (pbc)",
    )

    all_traces = (
        list(frame_traces)
        + [
            cell_trace,
            maxima_trace,
            surf_regions_trace,
            surf_pbc_trace,
            surf_cluster_trace,
            surf_cluster_trans,
            windows_trace,
            maxima_trace_pbc,
        ]
        + graph_periodic
    )

    # Ensure every trace is present and toggleable in legend.
    # Keep framework/cell on by default; park others in legend-only mode.
    initial_visible_names = {"Framework", "unit cell"}
    for i, t in enumerate(all_traces):
        n = getattr(t, "name", "")
        if not isinstance(n, str) or not n.strip():
            t.name = f"trace_{i}"
            n = t.name
        t.showlegend = True
        if (n in initial_visible_names) or ("framework" in n.lower()):
            t.visible = True
        else:
            t.visible = "legendonly"

    # Build camera-only frames so user can freely toggle legend traces while spinning.
    frames_anim: list[go.Frame] = []
    total_frames = max(int(frames), 1)
    fpr = max(float(frames_per_rotation), 1.0)
    ex0, ey0, ez0 = map(float, camera_eye)
    upx, upy, upz = map(float, camera_up)
    cx, cy, cz = map(float, camera_center)
    rxy = float(np.hypot(ex0, ey0))
    phi0 = float(np.arctan2(ey0, ex0))
    for i in range(total_frames):
        theta = (i / fpr) * 2.0 * np.pi
        phi = phi0 + theta
        cam = dict(
            eye=dict(x=float(rxy * np.cos(phi)), y=float(rxy * np.sin(phi)), z=ez0),
            up=dict(x=upx, y=upy, z=upz),
            center=dict(x=cx, y=cy, z=cz),
        )
        frames_anim.append(
            go.Frame(
                name=str(i),
                layout=dict(scene_camera=cam),
            )
        )
    init_cam = dict(
        eye=dict(x=ex0, y=ey0, z=ez0),
        up=dict(x=upx, y=upy, z=upz),
        center=dict(x=cx, y=cy, z=cz),
    )
    frame_duration_ms = int(round(1000.0 / max(float(fps), 0.001)))

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = go.Figure(
        data=all_traces,
        frames=frames_anim,
        layout=go.Layout(
            template="plotly_white",
            title=dict(
                text="Workflow traces loaded. Toggle traces in legend, then press Play.",
                x=0.5,
                xanchor="center",
                font=dict(size=16, color="black"),
            ),
            scene=dict(
                aspectmode="data",
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                camera=init_cam,
            ),
            margin=dict(l=0, r=0, t=50, b=0),
            legend=dict(
                x=0.99, y=0.99, xanchor="right", yanchor="top",
                bgcolor="rgba(255,255,255,0.75)",
                groupclick="toggleitem",
            ),
            updatemenus=[
                dict(
                    type="buttons", showactive=False,
                    y=0, x=0, xanchor="left", yanchor="bottom", pad=dict(t=40),
                    buttons=[
                        dict(
                            label="Play", method="animate",
                            args=[None, dict(frame=dict(duration=frame_duration_ms, redraw=redraw),
                                            transition=dict(duration=0), fromcurrent=True)],
                        ),
                        dict(
                            label="Pause", method="animate",
                            args=[[None], dict(frame=dict(duration=0, redraw=False),
                                              mode="immediate", transition=dict(duration=0))],
                        ),
                    ],
                )
            ],
            sliders=[
                dict(
                    active=0,
                    steps=[
                        dict(
                            args=[[fr.name],
                                 dict(frame=dict(duration=0, redraw=redraw),
                                      mode="immediate", transition=dict(duration=0))],
                            label=str(i), method="animate",
                        )
                        for i, fr in enumerate(frames_anim)
                    ],
                    x=0.05, y=0, len=0.9, xanchor="left", yanchor="top",
                    pad=dict(b=0, t=40),
                    currentvalue=dict(prefix="Frame: ", visible=True),
                )
            ],
        ),
    )
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CIF → Plotly HTML workflow animation (browser playback + screen record)."
    )
    parser.add_argument("cif", type=Path, help="Path to CIF/ASE-readable structure file.")
    parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Output HTML path (default: <stem>_workflow_animation.html next to CIF).",
    )
    parser.add_argument(
        "--spacing", type=float, default=0.25,
        help="Distance-grid spacing in Å (default: 0.25).",
    )
    parser.add_argument("--frames", type=int, default=120, help="Total animation frames (default: 120).")
    parser.add_argument("--fps", type=float, default=12.0, help="Playback frames per second (default: 12).")
    parser.add_argument(
        "--frames-per-rotation",
        type=float,
        default=120.0,
        help="Number of frames per full 360-degree rotation (default: 120).",
    )
    parser.add_argument(
        "--eye",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=(2.0, 0.0, 1.2),
        help="Initial camera eye coordinates (default: 2.0 0.0 1.2).",
    )
    parser.add_argument(
        "--up",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=(0.0, 0.0, 1.0),
        help="Camera up vector (default: 0.0 0.0 1.0).",
    )
    parser.add_argument(
        "--center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=(0.0, 0.0, 0.0),
        help="Camera center/look-at point (default: 0.0 0.0 0.0).",
    )
    parser.add_argument(
        "--redraw",
        action="store_true",
        help="Force full redraw on each frame (usually slower). Default is off for smoother camera-only playback.",
    )
    parser.add_argument(
        "--cdn", action="store_true",
        help="Load plotly.js from CDN (smaller HTML file).",
    )
    args = parser.parse_args()

    cif_path = args.cif.resolve()
    if not cif_path.is_file():
        sys.exit(f"Not a file: {cif_path}")

    out = args.output.resolve() if args.output else cif_path.with_name(
        f"{cif_path.stem}_workflow_animation.html"
    )

    repo_root = _repo_root()
    sys.path.insert(0, str(repo_root))

    print(f"Reading: {cif_path}")
    atoms = read(str(cif_path))

    print(
        f"Building animation (spacing={args.spacing}, "
        f"frames={args.frames}, fps={args.fps}, "
        f"frames_per_rotation={args.frames_per_rotation})…"
    )
    fig = build_workflow_figure(
        atoms,
        spacing=args.spacing,
        frames=args.frames,
        fps=args.fps,
        frames_per_rotation=args.frames_per_rotation,
        camera_eye=tuple(args.eye),
        camera_up=tuple(args.up),
        camera_center=tuple(args.center),
        redraw=args.redraw,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    include_js = "cdn" if args.cdn else True
    fig.write_html(str(out), include_plotlyjs=include_js, full_html=True)
    print(f"Wrote: {out}")
    print("Open in browser, zoom as needed, press Play, then screen-record.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Build a self-contained Plotly HTML: framework + unit cell with animated camera orbit.

Open the HTML in a browser, adjust zoom with the mouse, then use **Play** (and optional
slider) and screen-record. No Kaleido / PNG export.

Example::

    python devtools/scripts/cif_camera_animation_html.py \\
        mofography/examples/TestCIFs/P1-Cu-BTC.cif \\
        -o /tmp/hkust_orbit.html
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from ase.io import read


def _repo_root() -> Path:
    # devtools/scripts/this_file.py -> repo root
    return Path(__file__).resolve().parent.parent.parent


def camera_legacy_z_spin(theta: float, radius: float = 2.0, z: float = 1.2) -> dict:
    """Original centered camera style: rotate eye in xy plane, keep z fixed/up."""
    return dict(
        eye=dict(
            x=float(radius * np.cos(theta)),
            y=float(radius * np.sin(theta)),
            z=float(z),
        ),
        center=dict(x=0.0, y=0.0, z=0.0),
        up=dict(x=0, y=0, z=1),
    )


def build_animation_figure(
    atoms,
    *,
    n_frames: int = 120,
    frame_duration_ms: int = 40,
    turns: float = 2.0,
    template: str = "plotly_white",
) -> go.Figure:
    import mofography as mgr

    atoms = mgr.update_with_covalent_bonds(atoms)
    traces = [
        *mgr.plot_framework_as_mesh(atoms, show_rings=False, stylize=True),
        mgr.plot_cell(atoms, color="black", name="unit cell"),
    ]

    # Original "first good" look: centered unit cell, in-place z-up spin.
    radius = 2.0
    z_eye = 1.2

    denom = max(n_frames - 1, 1)
    frames: list[go.Frame] = []
    for i in range(n_frames):
        theta = (i / denom) * turns * 2.0 * np.pi
        cam = camera_legacy_z_spin(theta, radius=radius, z=z_eye)
        frames.append(
            go.Frame(
                name=str(i),
                layout=dict(scene_camera=cam),
            )
        )

    initial_cam = camera_legacy_z_spin(0.0, radius=radius, z=z_eye)
    layout = go.Layout(
        template=template,
        title=dict(
            text="Camera orbit (Play) — zoom/pan with mouse first if you like",
            x=0.5,
            xanchor="center",
        ),
        scene=dict(
            aspectmode="data",
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            camera=initial_cam,
        ),
        margin=dict(l=0, r=0, t=50, b=0),
        updatemenus=[
            dict(
                type="buttons",
                showactive=False,
                y=0,
                x=0,
                xanchor="left",
                yanchor="bottom",
                pad=dict(t=40),
                buttons=[
                    dict(
                        label="Play",
                        method="animate",
                        args=[
                            None,
                            dict(
                                frame=dict(
                                    duration=frame_duration_ms,
                                    redraw=True,
                                ),
                                transition=dict(duration=0),
                                fromcurrent=True,
                            ),
                        ],
                    ),
                    dict(
                        label="Pause",
                        method="animate",
                        args=[
                            [None],
                            dict(
                                frame=dict(duration=0, redraw=False),
                                mode="immediate",
                                transition=dict(duration=0),
                            ),
                        ],
                    ),
                ],
            )
        ],
        sliders=[
            dict(
                active=0,
                steps=[
                    dict(
                        args=[
                            [fr.name],
                            dict(
                                frame=dict(duration=0, redraw=True),
                                mode="immediate",
                                transition=dict(duration=0),
                            ),
                        ],
                        label=str(i),
                        method="animate",
                    )
                    for i, fr in enumerate(frames)
                ],
                x=0.05,
                y=0,
                len=0.9,
                xanchor="left",
                yanchor="top",
                pad=dict(b=0, t=40),
                currentvalue=dict(prefix="Frame: ", visible=True),
            )
        ],
    )

    fig = go.Figure(data=traces, layout=layout, frames=frames)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CIF → Plotly HTML with camera-orbit animation (browser + screen record)."
    )
    parser.add_argument("cif", type=Path, help="Path to a structure file readable by ASE.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output .html path (default: <stem>_orbit.html next to the CIF).",
    )
    parser.add_argument(
        "--n-frames",
        type=int,
        default=120,
        help="Number of camera keyframes (more = smoother / larger HTML).",
    )
    parser.add_argument(
        "--frame-ms",
        type=int,
        default=40,
        help="Milliseconds per frame when using Play.",
    )
    parser.add_argument(
        "--turns",
        type=float,
        default=2.0,
        help="How many full 360° rotations over all frames.",
    )
    parser.add_argument(
        "--cdn",
        action="store_true",
        help="Load plotly.js from CDN (smaller HTML). Default bundles plotly.js.",
    )
    args = parser.parse_args()

    cif_path = args.cif.resolve()
    if not cif_path.is_file():
        sys.exit(f"Not a file: {cif_path}")

    out = args.output
    if out is None:
        out = cif_path.with_name(f"{cif_path.stem}_orbit.html")
    else:
        out = out.resolve()

    repo_root = _repo_root()
    sys.path.insert(0, str(repo_root))

    print(f"Reading {cif_path}")
    atoms = read(str(cif_path))

    print(
        f"Building figure ({args.n_frames} frames, {args.turns} turns, {args.frame_ms} ms/frame)…"
    )
    fig = build_animation_figure(
        atoms,
        n_frames=args.n_frames,
        frame_duration_ms=args.frame_ms,
        turns=args.turns,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    include_js = "cdn" if args.cdn else True
    fig.write_html(str(out), include_plotlyjs=include_js, full_html=True)
    print(f"Wrote {out}")
    print("Open in Chrome, zoom as needed, then Play and screen-record.")


if __name__ == "__main__":
    main()

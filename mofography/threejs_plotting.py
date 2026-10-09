"""Three.js HTML builders for MOFography (parity with mgr_cli_py threejs_plotting)."""
from __future__ import annotations

import os
import tempfile
from typing import Any, Dict, Optional

import numpy as np

from .movies import (
    MovieFigure,
    _make_trace,
    make_percolation_clusters_trace,
    make_percolation_windows_trace,
)


def show_percolation_threejs(
    analysis: Dict[str, Any],
    atoms,
    probe_radius: float,
    options: Optional[Dict[str, Any]] = None,
) -> str:
    """Build self-contained Three.js HTML for PLD / percolation tab parity."""
    fig = MovieFigure()
    opts = options or {}

    fig.add_cell(atoms)

    maxima = analysis.get("maxima")
    connections = analysis.get("connections", [])
    seg = {
        "dgrid": analysis["dgrid"],
        "regions": analysis["regions"],
        "all_maxima": maxima if maxima is not None else np.zeros((0, 3)),
        "connections": connections,
    }

    if connections:
        try:
            fig.add_connections(seg, atoms, hide_periodic=True)
        except Exception:
            pass

    if maxima is not None and len(maxima):
        try:
            fig.add_maxima(seg, atoms, hide_periodic=True)
        except Exception:
            pass

    pld_locs = analysis.get("pld_locations")
    if pld_locs is not None and len(pld_locs):
        cell = atoms.get_cell().array.T
        shape = np.asarray(analysis["dgrid"].shape, dtype=float)
        frac = np.asarray(pld_locs, dtype=float) / shape
        pos = np.dot(cell, frac.T).T
        fig.add_trace(
            _make_trace(
                trace_id="pld_locations",
                trace_type="points",
                name="PLD locations",
                data={
                    "positions": pos.reshape(-1).tolist(),
                    "colors": [[0.0, 1.0, 0.0]] * len(pos),
                    "sizes": [5.0] * len(pos),
                    "n_points": len(pos),
                },
                style="pld_locations",
                hover=[
                    f"PLD location<br>grid: [{int(g[0])}, {int(g[1])}, {int(g[2])}]"
                    for g in pld_locs
                ],
                visible=True,
            )
        )

    percolation_history = analysis.get("percolation_history", [])
    if percolation_history:
        cluster_traces = make_percolation_clusters_trace(
            percolation_history,
            analysis["dgrid"],
            analysis["regions"],
            atoms,
            probe_radius=probe_radius,
            isovalue=analysis.get("isovalue", 1.6),
            step_size=2,
            cluster_opacity=float(opts.get("plot_surface_opacity", 1.0)),
            taubin_smooth=bool(opts.get("plot_surface_taubin_smooth", False)),
        )
        fig.add_traces(cluster_traces)
        for trace in cluster_traces:
            tid = trace["id"]
            fig.groups[tid] = {"traces": [tid], "label": trace["name"]}

    if connections and percolation_history:
        window_traces = make_percolation_windows_trace(
            analysis,
            atoms,
            probe_radius=probe_radius,
            pld_tolerance=float(opts.get("pld_tolerance", 0.25)),
        )
        if window_traces.get("windows", {}).get("data", {}).get("n_points", 0):
            fig.add_trace(window_traces["windows"])
        if window_traces.get("pld_critical", {}).get("data", {}).get("n_points", 0):
            fig.add_trace(window_traces["pld_critical"])

    # Scene units are Cartesian Å: aim at the cell center and back off in proportion to the cell size.
    cell = np.asarray(atoms.get_cell().array, dtype=float)
    center = cell.T @ np.array([0.5, 0.5, 0.5])
    extent = float(np.linalg.norm(cell.sum(axis=0))) or 1.0
    view_dir = np.array([2.0, 0.0, 1.35])
    eye = center + view_dir / np.linalg.norm(view_dir) * 1.8 * extent
    fig.set_camera(
        eye=tuple(float(x) for x in opts.get("camera_eye", eye)),
        center=tuple(float(x) for x in opts.get("camera_center", center)),
        up=(0.0, 0.0, 1.0),
        projection="perspective",
    )

    layout = opts.get("layout", "dark")
    if layout == "white":
        fig.set_appearance(background="#ffffff", scene_bg=0xFFFFFF)
    else:
        fig.set_appearance(
            background="radial-gradient(circle at center, #7a7a7a 0%, #111111 100%)",
            scene_bg=0x111111,
        )

    tmp = tempfile.NamedTemporaryFile(suffix=".html", delete=False)
    try:
        fig.to_html(tmp.name)
        with open(tmp.name, "r", encoding="utf-8") as f:
            return f.read()
    finally:
        os.unlink(tmp.name)

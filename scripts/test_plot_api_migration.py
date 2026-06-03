#!/usr/bin/env python3
"""Smoke test: Plotly vs Three.js API migration (_3js suffixes, Figure3js)."""
from __future__ import annotations

import sys
from pathlib import Path

import plotly.graph_objects as go
from ase import Atoms
from ase.io import read

import mofography as mgr
from mofography import Figure3js, plot_cell_3js, plot_framework_as_mesh_3js
from mofography.viz_backend import is_threejs_trace, validate_threejs_traces


def _fail(msg: str) -> None:
    print(f"FAIL: {msg}")
    sys.exit(1)


def _ok(msg: str) -> None:
    print(f"OK: {msg}")


def check_plotly_on_mgr() -> None:
    if mgr.plot_cell.__module__ != "mofography.visualize":
        _fail(f"mgr.plot_cell should be Plotly (visualize), got {mgr.plot_cell.__module__}")
    _ok("mgr.plot_cell is Plotly (visualize)")


def check_3js_not_shadowing() -> None:
    if not hasattr(mgr, "plot_cell_3js"):
        _fail("mgr.plot_cell_3js missing")
    if mgr.plot_cell_3js is mgr.plot_cell:
        _fail("plot_cell_3js must differ from plot_cell")
    _ok("plot_cell_3js exported separately")


def check_wrong_trace_on_figure() -> None:
    fig3 = Figure3js()
    try:
        fig3.add_trace(go.Scatter3d(x=[0], y=[0], z=[0]))
        _fail("Plotly trace on Figure3js should raise TypeError")
    except TypeError as exc:
        if "Plotly" not in str(exc):
            _fail(f"unexpected error: {exc}")
    _ok("Plotly trace rejected on Figure3js")

    figp = go.Figure()
    threejs_trace = {
        "id": "t",
        "type": "mesh",
        "name": "x",
        "visible": True,
        "style": {},
        "data": {"positions": [], "indices": [], "n_vertices": 0, "n_triangles": 0},
        "_mofography_backend": "threejs",
    }
    try:
        mgr.add_traces_plotly(figp, threejs_trace)
        _fail("Three.js trace on go.Figure should raise TypeError")
    except TypeError as exc:
        if "Three.js" not in str(exc):
            _fail(f"unexpected error: {exc}")
    _ok("Three.js trace rejected on go.Figure (add_traces_plotly)")


def check_build_scene_minimal() -> None:
    cif_candidates = list(Path(".").glob("*.cif"))
    if cif_candidates:
        atoms = read(str(cif_candidates[0]))
    else:
        atoms = Atoms("C2", positions=[[0, 0, 0], [1.5, 0, 0]], cell=[5, 5, 5], pbc=True)

    mat_updated = mgr.update_with_covalent_bonds(atoms)

    cell_bundle = plot_cell_3js(atoms, color="orange", name="cell")
    if not is_threejs_trace(cell_bundle):
        _fail("plot_cell_3js should return a tagged Three.js trace dict")

    fw_bundle = plot_framework_as_mesh_3js(mat_updated, stylize=True, show_rings=False)
    validate_threejs_traces(fw_bundle)

    fig = (
        Figure3js(layout="white", stylize=True)
        .plot_framework_as_mesh_3js(mat_updated, stylize=True, show_rings=False)
        .plot_cell_3js(atoms, color="orange")
    )
    html = fig.to_html()
    if "<html" not in html.lower():
        _fail("Figure3js.to_html() did not return HTML")
    _ok(f"Figure3js scene built ({len(html)} bytes HTML)")


def check_plotly_figure_still_works() -> None:
    atoms = Atoms("C2", positions=[[0, 0, 0], [1.5, 0, 0]], cell=[5, 5, 5], pbc=True)
    mat_updated = mgr.update_with_covalent_bonds(atoms)
    frame = mgr.plot_framework_as_mesh(mat_updated, stylize=True, show_rings=False)
    cell = mgr.plot_cell(atoms, color="orange")
    fig = go.Figure(data=[*frame, cell], layout=mgr.white_layout)
    if len(fig.data) < 2:
        _fail("Plotly figure should have framework + cell traces")
    _ok("Plotly go.Figure with mgr.plot_* traces")


def main() -> None:
    print("=== mofography plot API migration tests ===\n")
    check_plotly_on_mgr()
    check_3js_not_shadowing()
    check_wrong_trace_on_figure()
    check_build_scene_minimal()
    check_plotly_figure_still_works()
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()

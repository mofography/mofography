#!/usr/bin/env python3
"""Cu-BTC label inspection HTML (Plotly + Three.js) per internal_docs/plotting_guide.org.

Layers: framework mesh, cell, regions surface, regions_pbc surface,
all_maxima[:regions.max()], maxima_pbc (representatives).
Same Pastel cscale on surfaces + maxima so rep label colors should align.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import mofography as mgr
from mofography import Figure3js

CIF = ROOT / "examples/TestCIFs/P1-Cu-BTC.cif"
OUT_DIR = ROOT / "examples/viz_outputs"
OUT_PLOTLY = OUT_DIR / "P1-Cu-BTC_label_inspection_plotly.html"
OUT_THREEJS = OUT_DIR / "P1-Cu-BTC_label_inspection_threejs.html"

SPACING = 0.5
MASK_THICKNESS = 0.05
H_VALUE = 0.2
MIN_DISTANCE = 1.5
THRESHOLD_ABS = 1.5
EXTRA_PADDING = 5.0
ISOVALUE = 1.0

_TAUBIN = dict(
    taubin_smooth=True,
    taubin_iterations=3,
    taubin_lambda=0.2175,
    taubin_nu=-0.2175,
    stepsize=1,
)


def _rep_maxima_coords(results) -> tuple[np.ndarray, np.ndarray]:
    """Unit-cell grid indices for rep labels 1..N from maxima_info vs maxima_pbc."""
    mi = np.asarray(results["maxima_info"], dtype=int)
    maxima_pbc = np.asarray(results["maxima_pbc"], dtype=int)
    rep_from_info = np.array(
        [mi[np.where(mi[:, 0] == lab)[0][0], 1:4] for lab in range(1, len(maxima_pbc) + 1)],
        dtype=int,
    )
    return rep_from_info, maxima_pbc


def _print_label_check(results) -> None:
    rep_info, maxima_pbc = _rep_maxima_coords(results)
    n = len(maxima_pbc)
    exact = int(np.sum(np.all(rep_info == maxima_pbc, axis=1)))
    print(f"Representative maxima: {n} reps")
    print(f"  maxima_info[label, 1:4] == maxima_pbc: {exact}/{n}")
    for i in range(min(n, 20)):
        ok = np.array_equal(rep_info[i], maxima_pbc[i])
        if not ok or i < 8:
            print(
                f"  rep {i + 1}: info={rep_info[i]} pbc={maxima_pbc[i]} match={ok}"
            )
    regions = np.asarray(results["regions"])
    regions_pbc = np.asarray(results["regions_pbc"])
    pg = results["pbc_groups"]
    n_dom = sum(
        np.array_equal(
            np.isin(regions, sorted(int(x) for x in g)),
            regions_pbc == min(g),
        )
        for g in pg
    )
    print(f"  domain union match regions_pbc[rep]: {n_dom}/{len(pg)}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cscale = px.colors.qualitative.Pastel

    print(f"Loading {CIF.name}...")
    mat_atoms = read(str(CIF))
    mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

    print("Distance grid + segmentation...")
    t0 = time.time()
    dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=SPACING).astype(np.float32)
    results = mgr.do_periodic_segmentation_rust(
        dgrid,
        mat_atoms,
        spacing=SPACING,
        mask_thickness=MASK_THICKNESS,
        h=H_VALUE,
        min_distance=MIN_DISTANCE,
        threshold_abs=THRESHOLD_ABS,
        extra_padding=EXTRA_PADDING,
        return_dgrid=True,
        compute_connections=True,
    )
    print(f"  segmentation done in {time.time() - t0:.1f}s")

    dgrid = results["dgrid"]
    regions = results["regions"]
    regions_pbc = results["regions_pbc"]
    all_maxima = results["all_maxima"]
    maxima_pbc = results["maxima_pbc"]
    n_cap = int(regions.max()) if regions.size else 0
    maxima_local = all_maxima[:n_cap] if n_cap > 0 else all_maxima[:0]

    _print_label_check(results)

    n_pbc = len(maxima_pbc)
    mi = np.asarray(results["maxima_info"], dtype=int)
    local_labels = mi[: len(maxima_local), 0] if len(maxima_local) else np.array([], dtype=int)
    pbc_colors = list(range(1, n_pbc + 1))

    print("Building Plotly figure...")
    frame_traces = mgr.plot_framework_as_mesh(
        mat_updated, stylize=True, show_rings=True
    )
    cell_trace = mgr.plot_cell(mat_atoms, color="orange", name="cell")
    surface_regions = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=ISOVALUE,
        ase_atoms=mat_atoms,
        region_labels=regions,
        opacity=0.55,
        cscale=cscale,
        stylize=True,
        name="Regions Surface",
        **_TAUBIN,
    )
    surface_regions_pbc = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=ISOVALUE,
        ase_atoms=mat_atoms,
        region_labels=regions_pbc,
        opacity=0.55,
        cscale=cscale,
        stylize=True,
        name="Regions Surface (PBC)",
        **_TAUBIN,
    )
    maxima_local_trace = mgr.plot_maxima(
        dgrid,
        maxima_local,
        mat_atoms,
        show_spheres=False,
        marker_size=9,
        marker_color=local_labels,
        marker_colorscale=cscale,
        name=f"All Maxima [:regions.max()={n_cap}]",
    )
    maxima_pbc_trace = mgr.plot_maxima(
        dgrid,
        maxima_pbc,
        mat_atoms,
        show_spheres=False,
        marker_size=12,
        marker_color=pbc_colors,
        marker_colorscale=cscale,
        name="Maxima (pbc representatives)",
    )

    fig_plotly = go.Figure(
        data=[
            *frame_traces,
            cell_trace,
            surface_regions,
            surface_regions_pbc,
            maxima_local_trace,
            maxima_pbc_trace,
        ],
        layout=mgr.white_layout,
    )
    fig_plotly.update_layout(
        title=(
            "Cu-BTC label inspection — rep maxima (large) should sit on same-color "
            "PBC pocket; local maxima [:regions.max()] mark rep grid indices"
        )
    )
    fig_plotly.write_html(str(OUT_PLOTLY))
    print(f"  Plotly: {OUT_PLOTLY} ({OUT_PLOTLY.stat().st_size / 1024 / 1024:.1f} MB)")

    print("Building Three.js figure...")
    fig3 = (
        Figure3js(layout="white", stylize=True)
        .plot_framework_as_mesh_3js(mat_updated, stylize=True, show_rings=True)
        .plot_cell_3js(mat_atoms, color="orange", name="cell")
        .plot_surface_dgrid_3js(
            dgrid,
            ISOVALUE,
            mat_atoms,
            region_labels=regions,
            cscale=cscale,
            opacity=0.55,
            name="Regions Surface",
            stylize=True,
            **_TAUBIN,
        )
        .plot_surface_dgrid_3js(
            dgrid,
            ISOVALUE,
            mat_atoms,
            region_labels=regions_pbc,
            cscale=cscale,
            opacity=0.55,
            name="Regions Surface (PBC)",
            stylize=True,
            **_TAUBIN,
        )
    )
    fig3.plot_maxima_from_grid_3js(
        dgrid,
        maxima_local,
        mat_atoms,
        regions=regions,
        cscale=cscale,
        name=f"All Maxima [:regions.max()={n_cap}]",
        maxima_info=results.get("maxima_info"),
        hide_periodic=True,
    )
    fig3.plot_maxima_from_grid_3js(
        dgrid,
        maxima_pbc,
        mat_atoms,
        regions=regions_pbc,
        cscale=cscale,
        name="Maxima (pbc representatives)",
        maxima_pbc=maxima_pbc,
        maxima_info=results.get("maxima_info"),
        maxima_pbc_dict=results.get("maxima_pbc_dict", {}),
        hide_periodic=False,
    )
    fig3.to_html(str(OUT_THREEJS))
    print(f"  Three.js: {OUT_THREEJS} ({OUT_THREEJS.stat().st_size / 1024 / 1024:.1f} MB)")
    print("Done. Toggle layers in legend; hover maxima for grid indices + label.")


if __name__ == "__main__":
    main()

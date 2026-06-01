#!/usr/bin/env python3
"""Cu-BTC HTML via Figure3JS — granular plotly-style API demo.

Usage::

    cd ~/repos/mofography
    python examples/gen_cubtc_figure3js.py

Produces ``cubtc_figure3js.html`` in the repo root.

The snippet at the bottom of this file is the user-facing analogue of the
plotly snippet the user provided, but using the Figure3JS backend:

    import plotly.express as px
    cscale = px.colors.qualitative.Pastel
    fig = mgr.Figure3JS()
    fig.plot_surface_dgrid(dgrid, isovals[0], mat_atoms,
                           region_labels=regions, name='Regions', opacity=1,
                           cscale=cscale)
    fig.plot_surface_dgrid(dgrid, isovals[0], mat_atoms,
                           region_labels=regions_pbc, name='Regions PBC',
                           opacity=0.35, cscale=cscale)
    fig.plot_cell(mat_atoms, color='orange')
    fig.plot_maxima(dgrid, maxima, mat_atoms, name='Maxima (distinct)',
                    cscale=cscale)
    fig.draw_framework(mat_atoms)
    fig.plot_connections_as_it_is(connections_trimmed, all_maxima,
                                   dgrid.shape, mat_atoms, cscale=cscale,
                                   plot_disconnected=False, plot_periodic=True)
    fig.plot_windows(connections_trimmed, name='Windows',
                     marker_color='yellow', marker_size=10, marker_opacity=0.8,
                     grid_shape=dgrid.shape, maxima=all_maxima)
    fig.set_camera(eye=(2.0, 0.0, 1.35))
    fig.to_html('cubtc_figure3js.html')
"""

import sys, time
from pathlib import Path

import numpy as np
import plotly.express as px

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import mofography as mgr
from mofography import threejs_viz as tjs

CIF  = ROOT / "mofography/examples/TestCIFs/P1-Cu-BTC.cif"
OUT  = ROOT / "cubtc_figure3js.html"
SPACING = 0.5
ISOVALS = [1.0]

# ---------------------------------------------------------------------------
# Colorscale (same as plotly snippet)
# ---------------------------------------------------------------------------
cscale = tjs.qualitative_colorscale("Pastel")

# ---------------------------------------------------------------------------
# Load and compute covalent bonds
# ---------------------------------------------------------------------------
print("Loading Cu-BTC …")
from ase.io import read
atoms_raw = read(str(CIF))
atoms = mgr.update_with_covalent_bonds(atoms_raw)
print(f"  {len(atoms)} atoms, {len(atoms.info.get('covalent_bonds', []))} bonds")

# ---------------------------------------------------------------------------
# Distance grid
# ---------------------------------------------------------------------------
print("Distance grid …")
t0 = time.time()
dgrid = mgr.dgrid_from_atoms_cpu(atoms, spacing=SPACING, block_size=20000).astype(np.float32)
print(f"  shape={dgrid.shape}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Segmentation
# ---------------------------------------------------------------------------
print("Segmentation …")
t0 = time.time()
seg = mgr.do_periodic_segmentation_rust(
    dgrid, atoms, spacing=SPACING,
    mask_thickness=0.5, h=0.2,
    threshold_abs=1.5, min_distance=1.5,
    compute_connections=True, return_dgrid=False, extra_padding=5,
)
print(f"  done in {time.time()-t0:.1f}s")

regions      = seg["regions"]
regions_pbc  = seg["regions_pbc"]
connections  = seg["connections"]
all_maxima   = seg["all_maxima"]
unique_regions = seg["unique_regions"]

connections_trimmed = np.array(
    [c for c in connections if np.any(np.isin(c[0], unique_regions))],
    dtype=object,
)

# ---------------------------------------------------------------------------
# Build Figure3JS — plotly-style granular API
# ---------------------------------------------------------------------------
print("Building Figure3JS …")

fig = tjs.Figure3JS(layout="dark", stylize=True)

fig.set_camera(eye=(2.0, 0.0, 1.35), center=(0, 0, 0), up=(0, 0, 1))

# Unit-cell wireframe
fig.plot_cell(atoms, color="orange")

# Framework (ball-and-stick mesh)
fig.draw_framework(atoms)

# Segmentation surfaces
fig.plot_surface_dgrid(
    dgrid, ISOVALS[0], atoms,
    region_labels=regions,
    name="Regions", opacity=1.0,
    cscale=cscale, stepsize=2, stylize=True,
)
fig.plot_surface_dgrid(
    dgrid, ISOVALS[0], atoms,
    region_labels=regions_pbc,
    name="Regions PBC", opacity=0.35,
    cscale=cscale, stepsize=2, stylize=True,
    trace_id="surface_regions_pbc",
)

# Maxima (distinct pore centres)
fig.plot_maxima(
    dgrid, all_maxima, atoms,
    name="Maxima (distinct)", cscale=cscale,
)

# Channel connections
fig.plot_connections_as_it_is(
    connections_trimmed, all_maxima, dgrid.shape, atoms,
    cscale=cscale, plot_disconnected=False, plot_periodic=True,
)

# Windows
fig.plot_windows(
    connections_trimmed,
    name="Windows",
    marker_color="yellow",
    marker_size=10,
    marker_opacity=0.8,
    grid_shape=dgrid.shape,
    maxima=all_maxima,
)

# Export
fig.to_html(str(OUT))
print(f"\nSaved {OUT}  ({OUT.stat().st_size / 1024 / 1024:.1f} MB)")
print("\n--- Equivalent plotly-style snippet ---")
print("""
import mofography as mgr
from mofography import threejs_viz as tjs
import plotly.express as px

cscale = tjs.qualitative_colorscale("Pastel")   # or pass px.colors.qualitative.Pastel

fig = tjs.Figure3JS()
fig.plot_cell(mat_atoms, color='orange')
fig.draw_framework(mat_atoms)                    # ball-and-stick
fig.plot_surface_dgrid(dgrid, isovals[0], mat_atoms,
    region_labels=regions, name='Regions', opacity=1, cscale=cscale)
fig.plot_surface_dgrid(dgrid, isovals[0], mat_atoms,
    region_labels=regions_pbc, name='Regions PBC', opacity=0.35, cscale=cscale,
    trace_id='surface_regions_pbc')
fig.plot_maxima(dgrid, maxima, mat_atoms,
    name='Maxima (distinct)', cscale=cscale)
fig.plot_connections_as_it_is(connections_trimmed, all_maxima, dgrid.shape, mat_atoms,
    cscale=cscale, plot_disconnected=False, plot_periodic=True)
fig.plot_windows(connections_trimmed, name='Windows', marker_color='yellow',
    marker_size=10, marker_opacity=0.8,
    grid_shape=dgrid.shape, maxima=all_maxima)
fig.set_camera(eye=(2.0, 0.0, 1.35))
fig.to_html('output.html')
""")

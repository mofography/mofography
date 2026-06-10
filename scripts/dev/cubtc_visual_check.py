#!/usr/bin/env python3
"""Generate a standalone HTML of CuBTC using the updated batched mesh functions."""

import sys
sys.path.insert(0, '/Users/arungopalan/repos/mofography')

import numpy as np
import plotly.graph_objects as go
import mofography as mgr
from ase.io import read
from ase.data import covalent_radii

# Load Cu-BTC
filename = '/Users/arungopalan/repos/mofography/mofography/examples/TestCIFs/P1-Cu-BTC.cif'
mat_atoms = read(filename)
print(f"Loaded {mat_atoms.get_global_number_of_atoms()} atoms")

# Compute bonds
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

# Framework mesh trace (batched path via mesh_for_atoms)
framework_traces = mgr.plot_framework_as_mesh(mat_updated, show_rings=True, stylize=False)

# Compute grid + regions + maxima + connections
spacing = 0.25
radii = covalent_radii[mat_updated.get_atomic_numbers()]
grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing=spacing, block_size=20000)
regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=1.0, threshold_abs=1)
connections = mgr.connections_from_regions_and_dgrid(
    region_labels=regions, maxima=maxima, dgrid=grid.compute(), ase_atoms=mat_updated
)

cage_radii = grid.compute()[tuple(maxima.T)]
cage_labels = [f"{c:.2f} Å" for c in cage_radii]

windows = np.vstack([c[3] for c in connections if c[1]])
window_labels1 = np.array(['-'.join((str(c[0][0]), str(c[0][1]))) for c in connections if c[1] for w in c[3] if c[1]])
window_labels2 = np.array([c[2] for c in connections if c[1] for w in c[3] if c[1]])
window_labels3 = np.array([str(c[4]) for c in connections if c[1] for w in c[3] if c[1]])
window_labels = [w1 + "<br>" + w2 + "<br>" + w3 for w1, w2, w3 in zip(window_labels1, window_labels2, window_labels3)]

# Build figure
fig = go.Figure(data=[
    mgr.plot_surface_dgrid(
        grid.compute(), isovalue=2.00, opacity=1, stepsize=1,
        region_labels=regions, select='all', ase_atoms=mat_updated,
        stylize=True, cscale=None, random_seed=123
    ),
    mgr.plot_cell(mat_updated),
    mgr.scatter3d(
        mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=mat_updated),
        name='local maxima',
        marker_color=grid.compute()[tuple(maxima.T)],
        marker_size=15,
        hovertext=[f"{i+1} {j}" for i, j in enumerate(cage_labels)]
    ),
    mgr.scatter3d(windows, name='windows', marker_color='yellowgreen', marker_size=10, opacity=0.5, hovertext=window_labels),
    *framework_traces,
], layout=mgr.dark_layout)

output_path = '/Users/arungopalan/repos/mofography/cubtc_batched_check.html'
fig.write_html(output_path)
print(f"HTML written to: {output_path}")

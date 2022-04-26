#%%
import sys
import os
import numpy as np 

import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default ='notebook_connected'
path_to_mgr = os.path.abspath('c:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
sys.path.append(path_to_mgr)

#%%
import mofography as mgr
from ase.io import read
mat_atoms = read('./LTA96.cif')
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
mat_graph = mgr.sgraph_from_atoms(mat_updated)



# * Test the aromatic detection
#%%
spacing=0.25
import mofography as mgr
grid = mgr.dgrid_from_atoms(mat_atoms, spacing =spacing, block_size=20000)
#%%
regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=0.5, exclude_border=False, min_distance=int(1/spacing), threshold_abs=2, p_norm=2)
#%%

#%%
connections = mgr.connections_from_regions(region_labels=regions, ase_atoms=mat_updated)
pgraph = mgr.graph_from_connections(connections)
#%%
bonds = mat_updated.info['covalent_bonds']
bond_orders = mat_updated.info['covalent_bond_orders']
go.Figure(data= [
                mgr.plot_surface_dgrid(grid.compute(), isovalue=1, region_labels = regions, ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                mgr.plot_cell(mat_updated),
                mgr.scatter3d(mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=mat_updated), name='local maxima',marker_color = grid.compute()[tuple(maxima.T)], hovertext=list(map(str, np.round(grid.compute()[tuple(maxima.T)],decimals=2))) ),
                *mgr.draw_framework2(mat_updated, bonds=bonds, bond_orders=bond_orders, show_aromatics=False, bond_color='grey', bond_opacity=0.5)
                ], layout= mgr.dark_layout)
# %%
hist, bins = np.histogram(grid.compute()[tuple(maxima.T)], bins=100)
binC = (bins[1:] + bins[:-1])/2
go.Figure(data = [go.Scatter(x=binC, y=hist, mode='lines', name='dhist')], layout=mgr.white_layout)
# %%

#%%
import sys
import os
import numpy as np 

import plotly.graph_objects as go
import plotly.io as pio
from regex import D
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
from ase.data import covalent_radii 

radii = covalent_radii[mat_updated.get_atomic_numbers()] 
grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing =spacing, block_size=20000)
#%%
regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=0.1, exclude_border=False, min_distance=int(0.5/spacing), threshold_abs=0.5, p_norm=2)
#%%

#%%
connections = mgr.connections_from_regions(region_labels=regions, ase_atoms=mat_updated, distance_threshold=1)
pgraph = mgr.graph_from_connections(connections)
#%%
bonds = mat_updated.info['covalent_bonds']
bond_orders = mat_updated.info['covalent_bond_orders']
#%%
windows = np.vstack([c[2] for c in connections if c[1]])
cage_radii = grid.compute()[tuple(maxima.T)]
import pandas as pd
cage_labels = pd.cut(cage_radii, bins=[0,2,4,7], labels = [ 'tiny','small', 'big'])
connect_cage_labels = [cage_labels[c[0]-1] for c in connections if c[1]]
#%%
cations = mgr.read_raspa_pdb('./Movie_LTA96_1.1.1_298.150000_0.000000_component_K_0.pdb')
#%%
interpolator = mgr.interpolate_labels(regions)
for config in cations['coords']:
    interpolator(mgr.get_fractional_coordinates(config, mat_atoms))
# %%
def get_cage_for_ions(coords, cage_labels, interpolator, ase_atoms):
    import numpy as np
    output_labels = []
    from rich.progress import track
    for coord in track(coords):
        # print(coords.shape)
        # print((interpolator(get_fractional_coordinates(coords, ase_atoms=nax)).astype(int)-1).shape)
        try:
            output_labels.append(np.array(cage_labels)[interpolator(mgr.get_fractional_coordinates(coord, ase_atoms=ase_atoms)).astype(int)-1])
        except:
            print(np.max(coord))
            print(np.min(coord))
            pass
    return np.vstack(output_labels)
#%%
fig = go.Figure(data= [
                mgr.plot_surface_dgrid(grid.compute(), isovalue=1.0, stepsize=1, region_labels = regions, ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                mgr.plot_cell(mat_updated),
                mgr.scatter3d(mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=mat_updated), name='local maxima',marker_color = grid.compute()[tuple(maxima.T)], hovertext= [i+' '+j for i,j in zip(cage_labels, list(map(str, np.round(grid.compute()[tuple(maxima.T)],decimals=2))))] ),
                mgr.scatter3d(windows, name='windows',marker_color = 'yellow', marker_size=15, opacity=0.5, hovertext=['-'.join((c[0],c[1])) for c in connect_cage_labels]),
                mgr.scatter3d(cations['coords'][500], name='ions', marker_color =  interpolator(mgr.get_fractional_coordinates(cations['coords'][500], ase_atoms=mat_atoms)).astype(int), marker_size=15, opacity=0.5, hovertext= [str(x) for x in interpolator(mgr.get_fractional_coordinates(cations['coords'][500], ase_atoms=mat_atoms)).astype(int)]),
                *mgr.draw_framework2(mat_updated, bonds=bonds, bond_orders=bond_orders, show_aromatics=False, bond_color='grey', bond_opacity=0.5)
                ], layout= mgr.dark_layout)
fig.write_html('lta.html')
# %%
# hist, bins = np.histogram(grid.compute()[tuple(maxima.T)], bins=100)
# binC = (bins[1:] + bins[:-1])/2
# go.Figure(data = [go.Scatter(x=binC, y=hist, mode='lines', name='dhist')], layout=mgr.white_layout)
# %%


# %%
[str(x) for x in interpolator(mgr.get_fractional_coordinates(cations['coords'][500], ase_atoms=mat_atoms)).astype(int)]
# %%

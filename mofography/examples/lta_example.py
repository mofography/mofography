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
connections = mgr.connections_from_regions(region_labels=regions, maxima=maxima, ase_atoms=mat_updated, distance_threshold=1)
pgraph = mgr.graph_from_connections(connections)
#%%
bonds = mat_updated.info['covalent_bonds']
bond_orders = mat_updated.info['covalent_bond_orders']
#%%
windows = np.vstack([c[3] for c in connections if c[1]]) # pick out the window coordinates from connections
cage_radii = grid.compute()[tuple(maxima.T)]
import pandas as pd
cage_labels = pd.cut(cage_radii, bins=[0,2,4,7], labels = [ 'tiny','small', 'big'])
connected_cage_labels = [cage_labels[c[0]-1] for c in connections if c[1]]
window_labels2 = np.array(['-'.join((c[0],c[1])) for c in connected_cage_labels])
window_labels1 = np.array(['-'.join((str(c[0][0]),str(c[0][1]))) for c in connections if c[1]])
window_labels3 = np.array([c[2] for c in connections if c[1]])
window_labels = [w1+"<br>"+w2+'<br>'+w3 for w1,w2,w3 in zip(window_labels1, window_labels2, window_labels3)]
#%%
cations = mgr.read_raspa_pdb('./Movie_LTA96_1.1.1_298.150000_0.000000_component_K_0.pdb')
#%%
interpolator = mgr.interpolate_labels(regions)
# for config in cations['coords']:
#     interpolator(mgr.get_fractional_coordinates(config, mat_atoms))
# %%
def get_cage_for_ions(coords, cage_labels, interpolator, ase_atoms):
    """
    > For each ion in the `coords` array, find the cage label at the fractional coordinates of the ion
    
    :param coords: a list of coordinates (in cartesian) of the ions you want to find the cages for
    :param cage_labels: a list of strings, each string is a label for a cage
    :param interpolator: This is the interpolator object that you created in the previous step
    :param ase_atoms: The ASE atoms object for the structure you're interested in
    :return: The cage labels for each ion.
    """
    import numpy as np
    output_labels = []
    for coord in coords:
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
index = 500
ion_region_labels = [str(x) for x in interpolator(mgr.get_fractional_coordinates(cations['coords'][index], ase_atoms=mat_atoms)).astype(int)]
import sparse
# markers_sparse = sparse.COO(maxima.T, 1, shape=grid.shape).todense()
fig = go.Figure(data= [

                mgr.plot_surface_dgrid(grid.compute(), isovalue=1.0, opacity=1,stepsize=1, region_labels = regions, select='all',ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                # mgr.plot_surface_dgrid(grid.compute(), isovalue=1.0, opacity=0.5,stepsize=1, region_labels = mgr.apply_pbc3(regions,maxima), select=[22,23,25,26,34],ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                # mgr.plot_surface_dgrid(grid.compute(), isovalue=1.0, opacity=1,stepsize=1, region_labels = mgr.apply_pbc3(regions, maxima), select=np.delete(np.unique(regions),[22,23,25,26,34]),ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                mgr.plot_cell(mat_updated),
                mgr.scatter3d(mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=mat_updated), name='local maxima',marker_color = grid.compute()[tuple(maxima.T)], hovertext= [i+' '+j for i,j in zip(cage_labels, list(map(str, np.round(grid.compute()[tuple(maxima.T)],decimals=2))))] ),
                mgr.scatter3d(windows, name='windows',marker_color = 'yellowgreen', marker_size=20, opacity=0.5, hovertext=window_labels),
                # mgr.scatter3d(wall_window_coordinates, name='wall_windows',marker_color = 'cyan', marker_size=24, opacity=0.75),
                mgr.scatter3d(cations['coords'][index], name='ions', marker_color =  interpolator(mgr.get_fractional_coordinates(cations['coords'][index], ase_atoms=mat_atoms)).astype(int), marker_size=15, opacity=1, hovertext= ion_region_labels),
                *mgr.draw_framework2(mat_updated, bonds=bonds, bond_orders=bond_orders, show_aromatics=False, bond_color='grey', bond_opacity=0.5)
                ], layout= mgr.dark_layout)
fig.write_html('lta.html')
    # %%
# hist, bins = np.histogram(grid.compute()[tuple(maxima.T)], bins=100)
# binC = (bins[1:] + bins[:-1])/2
# go.Figure(data = [go.Scatter(x=binC, y=hist, mode='lines', name='dhist')], layout=mgr.white_layout)
# %%
def nearest_window_and_distance(coords, ase_atoms, windows, window_labels=None, cut_off=1.0, n_neighbors=4):
    """
    > For each ion, find the nearest window and assign the ion to that window if the distance is less
    than the cut-off
    
    :param coords: the coordinates of the ions
    :param ase_atoms: the atoms object from ASE
    :param windows: a list of points in 3D space
    :param window_labels: a list of labels for each window. If None, then the indices of the windows are
    used
    :param cut_off: the maximum distance from the window center to the ion
    :param n_neighbors: the number of nearest neighbors to consider, defaults to 4 (optional)
    :return: The window labels of the windows that are within the cutoff distance of the ions.
    """
    ab = mgr.AABB_on_atoms(ase_atoms, points=windows)
    nls = ab.query(coords, query_args=dict(mode='nearest', exclude_ii=True, num_neighbors = n_neighbors )).toNeighborList(sort_by_distance=True)
    point_indices = nls.point_indices.reshape(-1,n_neighbors)
    distances = nls.distances.reshape(-1, n_neighbors)
    if window_labels is None:
        window_labels = list(range(len(windows)))   
    return np.vstack([[window_labels[i], d] for i,d in zip(point_indices[:,0], distances[:,0]) ]).astype('object')

#%%
labeled_ions =nearest_window_and_distance(np.vstack(cations['coords']), ase_atoms=mat_atoms, windows=windows, window_labels=window_labels, cut_off=1.5)
#%%
from ase.io import read, write
write('./mat_atoms.cube',images=mat_atoms, data=grid)
# %% 

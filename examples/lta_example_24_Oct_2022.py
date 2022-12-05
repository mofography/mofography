
#%%
import sys
import os
import numpy as np 

import plotly.graph_objects as go
import plotly.io as pio
from regex import D
pio.renderers.default ='browser'


path_to_mgr = os.path.abspath('c:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
sys.path.append(path_to_mgr)

#%% [markdown]
# Let's do all the framework stuff first
#%%
import mofography as mgr
from ase.io import read
mat_atoms = read('./TestCIFs/LTA96.cif')
# mat_atoms = read('./Framework_0_initial_1_1_1_P1.cif')


# let's update with covalent bonds and add some extra bonds across the 
# periodic boundaries
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

# Detect the bonds in the framework, this maybe useful for plotting
bonds = mat_updated.info['covalent_bonds']
bond_orders = mat_updated.info['covalent_bond_orders']

# make a structure graph, this is useful for plotting a part of of the 
# framework in the future
mat_graph = mgr.sgraph_from_atoms(mat_updated)


# frame traces 
# in the line below rings refer to the aromatic rings only
frame_traces= mgr.plot_framework_as_mesh(mat_updated, show_rings= True, stylize =True)#, 


#%% # Compute the grid
spacing=0.25
import mofography as mgr
from ase.data import covalent_radii 

radii = covalent_radii[mat_updated.get_atomic_numbers()] 
grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing =spacing, block_size=20000)

#%% # Compute the region segmentation
regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=0.5)
# Plot trace for the maxima 
maxima_trace= mgr.plot_maxima(grid.compute(), maxima, ase_atoms= mat_updated)

#%%
#* Find the windows and connected regions including the windows on the walls
# connections = mgr.connections_from_regions(dgrid = grid.compute(), region_labels=regions, maxima=maxima, ase_atoms=mat_updated)
rag = mgr.get_rag_from_regions_and_grid(regions, grid.compute(), mat_updated, maxima)
connections = mgr.get_connections_for_rag(rag, regions, maxima, mat_updated, grid.compute())
# pgraph = mgr.graph_from_connections(connections, maxima)     

# Plot trace for the windows 
window_trace = mgr.plot_windows(connections)
#%% #* Extract the window coordinates from the previus output, create labels for the 
#* cages and windows based on the cage sizes; tiny, small and big (1.44, 3.9, 6.2 )
import pandas as pd
windows = np.vstack([c[3] for c in connections if c[1]]) # pick out the window coordinates from connections
cage_radii = grid.compute()[tuple(maxima.T)]
cage_labels = pd.cut(cage_radii, bins=[0,2,5,7], labels = [ 'tiny','small', 'big'])
connected_cage_labels = [cage_labels[np.array(c[0])-1] for c in connections if c[1]]
window_labels2 = np.array(['-'.join((c[0],c[1])) for c in connected_cage_labels])
window_labels1 = np.array(['-'.join((str(c[0][0]),str(c[0][1]))) for c in connections if c[1]])
window_labels3 = np.array([c[2] for c in connections if c[1]])
window_labels = [w1+"<br>"+w2+'<br>'+w3 for w1,w2,w3 in zip(window_labels1, window_labels2, window_labels3)]



#%% [markdown]
# Let's do the cation stuff now

#%%
#%% #* read the cation coordinates
# cations = mgr.read_raspa_pdb('./Movie_LTA96_1.1.1_298.150000_0.000000_component_K_0.pdb')
cations = mgr.read_raspa_pdb('./Movie_LTA_my_1.1.1_303.000000_0.000000_component_Na_1.pdb')
# cations = mgr.read_raspa_pdb('./Movie_NaX1.0newmodel2_1.1.1_298.150000_100000.000000_component_Na_0.pdb')

# Create an linear interpolator for the regions for assigning arbitrary points (indices) to regions
interpolator = mgr.interpolate_labels(regions)
# for config in cations['coords']:
#     interpolator(mgr.get_fractional_coordinates(config, mat_atoms))
# %%
# * Assign ions to cages
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


#%% #* Assign ions to cages, we will use the cage labels while plotting the ions
index = 500
ion_region_labels = [str(x) for x in interpolator(mgr.get_fractional_coordinates(cations['coords'][index], ase_atoms=mat_atoms)).astype(int)]

# %% #* Assign ion to the nearest window
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
    :return: The window label and the distance of the window closest to each ion as long as it is within the cutoff distance, otherwise says 'not close enough'.
    """
    ab = mgr.AABB_on_atoms(ase_atoms, points=windows)
    nls = ab.query(coords, query_args=dict(mode='nearest', exclude_ii=True, num_neighbors = n_neighbors )).toNeighborList(sort_by_distance=True)
    point_indices = nls.point_indices.reshape(-1,n_neighbors)
    distances = nls.distances.reshape(-1, n_neighbors)
    if window_labels is None:
        window_labels = list(range(len(windows)))   
    return np.vstack([[window_labels[i], d] if d <= cut_off else ['not_close_enough', d] for i,d in zip(point_indices[:,0], distances[:,0]) ]).astype('object')

#%% #* Label the ions with the label of the nearest window (not_close_enough for ions that are too far)
# this is for all the cations in all the snapshots
labeled_ions =nearest_window_and_distance(np.vstack(cations['coords']), ase_atoms=mat_atoms, windows=windows, window_labels=window_labels, cut_off=1.5)


#%%
#* Print out the important ion fractions
print('Fraction too far from windows: {0}'.format(np.mean(np.mean(labeled_ions[:,0]=='not_close_enough'))))
#* Fraction of ions near the big-big windows
print('Fraction near big-big windows: {0}'.format(np.mean([('big-big' in l) for l in labeled_ions[:,0]])))
print('Fraction near small-big windows: {0}'.format(np.mean([np.logical_or('small-big' in l, 'big-small' in l) for l in labeled_ions[:,0]])))
print('Fraction near small-tiny windows: {0}'.format(np.mean([np.logical_or('small-tiny' in l, 'tiny-small' in l) for l in labeled_ions[:,0]])))
print('Fraction near big-tiny windows: {0}'.format(np.mean([np.logical_or('big-tiny' in l, 'tiny-big' in l) for l in labeled_ions[:,0]])))
# np.mean([('big-big' in l) for l in labeled_ions[:,0]])

# %% 
# this is just for snapshot 500.
index = 500
# ion_region_labels = [str(x) for x in interpolator(mgr.get_fractional_coordinates(cations['coords'][index], ase_atoms=mat_atoms)).astype(int)]

labeled_ions1 = nearest_window_and_distance(cations['coords'][index], ase_atoms=mat_atoms, windows=windows, window_labels=window_labels, cut_off=1.5) 
import plotly

frame_traces= mgr.plot_framework_as_mesh(mat_updated, show_rings=True, stylize =True)#, 


#%%
surface_trace = mgr.plot_surface_dgrid(grid.compute(), isovalue=1, opacity=1,stepsize=1, region_labels = regions, select='all',ase_atoms = mat_updated, cscale = plotly.colors.qualitative.Pastel_r ,stylize=True, random_seed=123), 
#%% [markdown]
# let's ake a mesh trace ofteh cations also, and then we can stylize it

#%%

cation_spheres= mgr.spheres_for_positions(cations['coords'][index],radii=0.8*np.ones(len(cations['coords'][index])) )

cations_trace = mgr.plot_mesh(cation_spheres,  stylize=True, color='magenta',opacity=1, name='K+ ions')


#%%
fig = go.Figure(data= [

                mgr.plot_cell(mat_updated, color='blue'),
                *frame_traces,
                *surface_trace,
                # mgr.plot_surface_dgrid(grid.compute(), isovalue=1, opacity=1,stepsize=1, region_labels = regions, select='all',ase_atoms = mat_updated, cscale = plotly.colors.qualitative.Pastel_r ,stylize=True, random_seed=123), 
                # mgr.plot_surface_dgrid(grid.compute(), isovalue=1.0, opacity=0.5,stepsize=1, region_labels = mgr.apply_pbc3(regions,maxima), select=[22,23,25,26,34],ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                # mgr.plot_surface_dgrid(grid.compute(), isovalue=1.0, opacity=1,stepsize=1, region_labels = mgr.apply_pbc3(regions, maxima), select=np.delete(np.unique(regions),[22,23,25,26,34]),ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                # mgr.scatter3d(mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=mat_updated), name='local maxima',marker_color = grid.compute()[tuple(maxima.T)], marker_size=15, hovertext= [i+' '+j for i,j in zip(cage_labels, list(map(str, np.round(grid.compute()[tuple(maxima.T)],decimals=2))))] ),
                # mgr.scatter3d(windows, name='windows',marker_color = 'yellowgreen', marker_size=10, opacity=0.5, hovertext=window_labels),
                # maxima_trace,
                # window_trace,
                # mgr.scatter3d(cations['coords'][index], name='ions', marker_color =  'yellow', marker_size=15, marker_line_color='white', opacity=0.8, hovertext=labeled_ions1[:,0]),
                # *mgr.draw_framework2(mat_updated, bonds=bonds, bond_orders=bond_orders, show_aromatics=False, bond_color='grey', bond_opacity=0.5)
                cations_trace
                ], layout= mgr.white_layout)
# fig.write_html('lta.html')
cc = 0.12
camera = dict( up=dict(x=1, y=0, z=0), eye=dict(zip(['x','y','z'],[ cc,cc,1.25])), center=dict(zip(['x', 'y','z'],[cc,cc,0])))

fig.update(layout=dict(scene=dict(camera=camera))).show()
#%%
def atom_indices_around_point(point, ase_atoms, r_max=10):
    
    import networkx as nx
    
    ab = mgr.AABB_on_atoms(ase_atoms)
    ab_r = ab.query(point, query_args=dict(mode='ball', r_max=r_max, exclude_ii=True)).toNeighborList(sort_by_distance=True)
    nls = ab_r.point_indices # these are the neighbors of the window
    # let;s make a subgraph of the neighbors
    return nls 
# %%

# # selection = [56, 79]
# selection = [56, 57]
# # selection = [34, 79]

selection = [25, 30]
# selection = [25, 45]
# selection = [25,46]
w = np.where(['-'.join(np.array(selection).astype(str)) in l for l in window_labels])[0]
# w = np.where(['56-57' in l for l in window_labels])[0]
# w = np.where(['34-79' in l for l in window_labels])[0]
atoms_near_window = atom_indices_around_point(windows[w[0]], mat_updated, r_max=5)
#%% # * Plot a subset of the atoms 
def plot_subset_of_atoms(ase_atoms, subset, **kwargs):

    import networkx as nx
    from ase import Atoms
    mat_graph = mgr.sgraph_from_atoms(ase_atoms)
    sg = mat_graph.subgraph(subset)
    sg = nx.relabel.convert_node_labels_to_integers(sg, first_label=0, ordering='default')
    cell = ase_atoms.get_cell()
    #make atoms object from the subgraph
    atoms_ = Atoms(symbols=list(nx.get_node_attributes(sg, 'element').values()), positions = list(nx.get_node_attributes(sg, 'position').values()), cell=cell)
    atoms_.info['atom ids'] = list(nx.get_node_attributes(sg, 'element').keys())
    # bonds = [list(e) for e in sg.edges()]
    # atoms.info['covalent_bonds'] = bonds

    # lets put some defualt values in for the kwargs
    stylize = kwargs.get('stylize', True)
    name= kwargs.get('name', 'subset')

    kwargs.pop('stylize', None)
    kwargs.pop('name', None)


    atoms_mesh = mgr.mesh_for_atoms(atoms_, bonds=np.vstack(list(sg.edges())))#, bond_orders=np.ones(len(sg.edges())), show_aromatics=False, bond_color='grey', bond_opacity=0.5)

    atoms_trace = mgr.plot_mesh([atoms_mesh], stylize = stylize, name=name, **kwargs)

    return atoms_trace
#%%
surface_trace2= mgr.plot_surface_dgrid(grid.compute(), isovalue=1, opacity=1,stepsize=1, region_labels = regions, select=selection,ase_atoms = mat_updated, cscale = plotly.colors.qualitative.Pastel_r ,stylize=True, random_seed=123), 
# fig = go.Figure(data= [mgr.plot_mesh([atoms_mesh], stylize=True, name='Part of framework'), mgr.plot_cell(atoms_)], layout= mgr.white_layout).show()
ww = windows[w[0]]
fig = go.Figure(data= [  *surface_trace2,go.Scatter3d(x=[ww[0]], y=[ww[1]], z=[ww[2]],marker_size=25, marker_color='lightgreen'), plot_subset_of_atoms(mat_updated, atoms_near_window.tolist(), name='Part of framework')], layout= mgr.white_layout).update_layout(scene=dict(camera=camera)).show()


# %%

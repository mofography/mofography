#%%
import sys
import os
import numpy as np 
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default ='browser'
# pio.renderers.default ='browser'

# * change this to path on your machine
path_to_mgr = os.path.abspath('c:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
sys.path.append(path_to_mgr)

#%%
import mofography as mgr
from ase.io import read
# filename = './TestCIFs/NU-1000_S1.cif'
# filename = './LTA96.cif'
# filename = './TestCIFs/MIL-47-RASPA.cif'
# filename = './TestCIFs/P1-Cu-BTC.cif'
filename ='./TestCIFs/P1-IRMOF-1.cif'
mat_atoms = read(filename)
#%%
frame_mesh = mgr.mesh_for_atoms(mat_atoms)
frame_mesh_trace = mgr.plot_framework_mesh(frame_mesh, stylize=True)
#%%
#* If you want to make supercell for some reason, uncomment this line
# * but, it will slow down the calculation might be better to replicate at the end
# import ase
# mat_atoms = ase.build.make_supercell(mat_atoms, P=2*np.eye(3))

#%% #* Detect the bonds in the framework for plotting
# * you can also get bonds and bond orders as shown below
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
mat_graph = mgr.sgraph_from_atoms(mat_updated)
bonds = mat_updated.info['covalent_bonds']
bond_orders = mat_updated.info['covalent_bond_orders']



#%% 
#* Compute the grid
spacing=0.25
import mofography as mgr
from ase.data import covalent_radii 

# * You can use whatever radii you want for the atoms: covalent, ionic, vdw
radii = covalent_radii[mat_updated.get_atomic_numbers()] 
grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing =spacing, block_size=20000)
#%% 
#* Compute the region segmentation and the size of the maxima (cages/pockets)
regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=1.0, exclude_border=False, min_distance=int(2/spacing), threshold_abs=1, p_norm=np.inf)
maxima_radii = grid.compute()[tuple(maxima.T)]

#%% 
#* Find the windows and connected regions including the windows on the walls
connections = mgr.connections_from_regions(region_labels=regions, maxima=maxima, dgrid=grid.compute(), ase_atoms=mat_updated,minimum_window_separation=1, wall_windows=True)
#%% 
# #* Create a pore graph (optional, comment corresponding line in fig)
G = mgr.graph_from_connections(connections, maxima)
G = mgr.add_props_to_graph(G, node_props = [mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=mat_updated)], node_prop_names=['node_coords'])
#%% 
# * Let's make some cage labels with the cage radii (same as maxima_radii)
cage_radii = grid.compute()[tuple(maxima.T)]
cage_labels = [str(c)+' \u212b' for c in np.round(cage_radii, decimals =2)]
#%%
#* Let's make some window labels for showing on hover as information
windows = np.vstack([c[3] for c in connections if c[1]])
window_labels1 = np.array(['-'.join((str(c[0][0]),str(c[0][1]))) for c in connections if c[1] for w in c[3] if c[1]])
window_labels2 = np.array([c[2] for c in connections if c[1] for w in c[3] if c[1]])
window_labels3 = np.array([str(c[4]) for c in connections if c[1] for w in c[3] if c[1]])
window_labels = [w1+"<br>"+w2+"<br>"+w3 for w1,w2,w3 in zip(window_labels1, window_labels2, window_labels3)]
#%% 
#* Plot everything
fig = go.Figure(data= [

                mgr.plot_surface_dgrid(grid.compute(), isovalue=2.00, opacity=1,stepsize=1, region_labels = regions, select='all',ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                mgr.plot_cell(mat_updated),
                mgr.scatter3d(mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=mat_updated), name='local maxima',marker_color = grid.compute()[tuple(maxima.T)], marker_size=15, hovertext= [str(i+1)+' '+j for i,j in enumerate(cage_labels)] ),
                mgr.scatter3d(windows, name='windows',marker_color = 'yellowgreen', marker_size=10, opacity=0.5, hovertext=window_labels),
                *mgr.plot_pore_graph(G),
                frame_mesh_trace,
                # *mgr.draw_framework2(mat_updated, bonds=bonds, bond_orders=bond_orders, show_aromatics=False, bond_color='grey', bond_opacity=0.5)
                ], layout= mgr.dark_layout)
# fig.write_html('lta.html')
# fig.show()


# %% 
#* Get the periodic map that would transform the regions into regions with PBC
# pmap = dict(zip(regions.flatten(), apply_pbc(regions, maxima, maxima_radii, wall_thickness=int(0.5/spacing)).flatten()))#mgr.get_periodic_map(regions, maxima, wall_thickness=0.01)
pmap = mgr.get_periodic_map(regions, maxima, maxima_radii, wall_thickness=3)

#%% 
#* new window labels with PBC mapping applied
windows = np.vstack([c[3] for c in connections if c[1]])
wregs = mgr.apply_periodic_map(np.vstack([c[0] for c in connections if c[1] for w in c[3]]), periodic_map=pmap)
wregs_labels  =np.array(['-'.join((str(w[0]),str(w[1]))) for w in wregs]) # regions
window_labels2 = np.array([c[2] for c in connections if c[1] for w in c[3] if c[1]]) # internal or periodic
window_labels3 = np.array([str(c[4])+ ' \u212b' for c in connections if c[1] for w in c[3] if c[1]]) # size
window_labels = [w1+"<br>"+w2+"<br>"+w3 for w1,w2,w3 in zip(wregs_labels, window_labels2, window_labels3)]

# %% #* Plot again with PBC map applied
fig = go.Figure(data= [

                mgr.plot_surface_dgrid(grid.compute(), isovalue=1.75, opacity=1,stepsize=1, region_labels = mgr.apply_periodic_map(regions, pmap), select='all',ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                # mgr.plot_surface_dgrid(grid.compute(), isovalue=1.5, opacity=1,stepsize=1, region_labels =regions_pbc, select='all',ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                # Plotting the regions with periodic boundary conditions.
                # mgr.plot_surface_dgrid(grid.compute(), isovalue=1.5, opacity=1,stepsize=1, region_labels =apply_pbc(regions, maxima, wall_thickness=0.02), select='all',ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                # mgr.plot_surface_dgrid(grid.compute(), isovalue=1.5, opacity=1,stepsize=1, region_labels = mgr.apply_pbc(mgr.apply_pbc(mgr.apply_pbc(regions, maxima, wall_thickness=0.01, reindex=True), maxima=maxima, wall_thickness=0.01, reindex=True), maxima), select='all',ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                # mgr.plot_surface_dgrid(grid.compute(), isovalue=1.5, opacity=0.5,stepsize=1, region_labels = apply_pbc(regions,maxima), select=[22,23,25,26,34],ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                # mgr.plot_surface_dgrid(grid.compute(), isovalue=1.0, opacity=1,stepsize=1, region_labels = mgr.apply_pbc3(regions, maxima), select=np.delete(np.unique(regions),[22,23,25,26,34]),ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                mgr.plot_cell(mat_updated),
                mgr.scatter3d(mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=mat_updated), name='local maxima',marker_color = grid.compute()[tuple(maxima.T)], marker_size=15, hovertext= [str(mgr.apply_periodic_map(i+1, pmap))+' '+j for i,j in enumerate(cage_labels)] ),
                # mgr.scatter3d(mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=mat_updated), name='local maxima',marker_color = grid.compute()[tuple(maxima.T)], marker_size=15, hovertext= [str(i+1)+' '+j for i,j in enumerate(cage_labels)] ),
                mgr.scatter3d(windows, name='windows',marker_color = 'yellowgreen', marker_size=10, opacity=0.5, hovertext=window_labels),
                # mgr.scatter3d(wall_window_coordinates, name='wall_windows',marker_color = 'cyan', marker_size=24, opacity=0.75),
                # mgr.scatter3d(cations['coords'][index], name='ions', marker_color =  interpolator(mgr.get_fractional_coordinates(cations['coords'][index], ase_atoms=mat_atoms)).astype(int), marker_size=20, opacity=1, hovertext= ion_region_labels),
                *mgr.plot_pore_graph(G),
                frame_mesh_trace,
                *mgr.draw_framework2(mat_updated , bonds=bonds, bond_orders=bond_orders, show_aromatics=False,bond_color='white', bond_opacity=0.25)
                ], layout= mgr.dark_layout)

# * Uncomment this line to write to file
# fig.write_html(filename.split('/')[-1].split('.')[0]+'.html')

# fig.update_layout(paper_bgcolor='darkred')
fig.show()

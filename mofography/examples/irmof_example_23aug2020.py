#%%
#region #* make imports
import sys
import os
# from types import _Cell
import numpy as np

import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default ='notebook_connected'

# * Put your path here
path_to_mgr = os.path.abspath('c:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
sys.path.append(path_to_mgr)
#endregion
#%%
#region #* Read cif and compute bonds
import mofography as mgr
from ase.io import read
# Reading the cif file and storing it in a variable called mat_atoms.
# mat_atoms = read('./TestCIFs/LTA96.cif')
# mat_atoms = read('./TestCIFs/NU-1000_S1.cif')
mat_atoms = read('./TestCIFs/NU-1000_S1.cif')
# mat_atoms = read('./TestCIFs/ZIF-8-RASPA.cif')
# mat_atoms = read('./TestCIFs/P1-Cu-BTC.cif')
# mat_atoms = read('./TestCIFs/P1-IRMOF-1.cif')
# mat_atoms = read('./TestCIFs/MIL-47-P1-RASPA.cif')
# mat_atoms = read('./TestCIFs/VUSKEA_clean.cif')
# mat_atoms = read('./TestCIFs/hypotheticalMOF_1003830_0_0_4_10_10_0.cif')
# mat_atoms = read('./TestCIFs/hypotheticalMOF_2000756_0_0_2_32_32_0.cif')
# mat_atoms = read('./TestCIFs/hypotheticalMOF_5026325_3_0_0_27_27_0.cif')
# mat_atoms = read('./TestCIFs/hypotheticalMOF_7002752_0_0_4_28_28_0.cif')
# mat_atoms = read('./TestCIFs/ADARAA_clean.cif')
# mat_atoms = read('./TestCIFs/ZIF-69-RASPA.cif')
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
mat_graph = mgr.sgraph_from_atoms(mat_updated)
#endregion

# %%
#region #* compute the distance grid
spacing=0.25
import mofography as mgr
from ase.data import covalent_radii

radii = covalent_radii[mat_updated.get_atomic_numbers()]
grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing =spacing, block_size=20000)
#endregion
#%%
##region #* Compute region segmentation and local maxima
regions, maxima = mgr.regions_from_dgrid(grid.compute())#, mask_thickness=1.5,exclude_border=False, min_distance=10 , threshold_abs=3, p_norm=2)
#endregion
#%%
#region #* Compute the poregraph
rag2 = mgr.get_rag_from_regions_and_grid(regions, grid.compute(), mat_updated,maxima)
#endregion
#%%
#region #* Useful PBC information in groups and as mapping
pbc_groups = mgr.get_periodic_map(regions=regions, maxima=maxima, maxima_radii=grid.compute()[tuple(maxima.T)], return_conn=True, wall_thickness=4)
pmap= mgr.get_periodic_map(regions=regions, maxima=maxima, maxima_radii=grid.compute()[tuple(maxima.T)], return_conn=False, wall_thickness=4)
#endregion
#%%
#region# * Make traces for surface, pore graph, unit cell
rag_trace=mgr.plot_rag(rag2,plot_with='maxima', name='p-graph') 
# surface_trace=plot_surface_dgrid2(grid.compute(), 2.0, select=list(pbc_groups[0].union(pbc_groups[2].union(pbc_groups[-3]))),ase_atoms=mat_updated,region_labels=regions, stylize=True)
surface_trace=mgr.plot_surface_dgrid(grid.compute(), 2.0, select='all',ase_atoms=mat_updated,region_labels=regions, stylize=True)
cell_trace=mgr.plot_cell(mat_updated)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, show_rings=False, stylize=True)
#endregion
#%%
#region # * Plot the figure with all the traces
# fig = go.Figure(data=[*frame_trace, cell_trace, surface_trace, rag_trace]).update(layout=mgr.dark_layout)
# fig.show()  
#endregion
#%%
#region #* Let's add some properties to the pore graph (rag) and visualize
from skimage.measure import regionprops
props = regionprops(regions)
# * size dist about centroid with and without pbc 
rag2 = mgr.add_psd_to_rag(rag2,props,  bins=np.linspace(0,30,50), ase_atoms=mat_updated, shape = regions.shape)
rag2 = mgr.add_psd_to_rag_with_pbc(rag=rag2,props =props, pbc_groups = pbc_groups, bins=np.linspace(0,30,50), ase_atoms=mat_updated, shape = regions.shape)
fig1=go.Figure() 
for group in pbc_groups:
    node = list(group)[0]
    fig1.add_trace(go.Scatter(x=rag2.nodes[node]['psd_x_pbc'], y=rag2.nodes[node]['psd_y_pbc'], mode='markers+lines', name=str(rag2.nodes[node]['labels'][0])))
pio.renderers.default = "notebook_connected"
fig1.update_layout(xaxis=dict(title='Cavity size (\u212b)'), yaxis=dict(title='Density histogram')).show()
# fig1.show()
#
#endregion
#%%
#region #* Let's find the unique 1D and 2D channels in the MOF
pca_list, paths_list= mgr.find_paths(regions, rag2)

# * separate into 1D and 2D paths
threshold = 0.99 #* This is our strict criterio for something to be 1D
paths1d_to_plot = np.array(paths_list)[np.where([check[0]>threshold for check in pca_list])[0]].tolist()
pathshd_to_plot = np.array(paths_list)[np.where([check[0]<threshold for check in pca_list])[0]].tolist()
#
#endregion

#%%
#region #* # Finding the isomorphic path groups. (paths that are truly different) 
path_groups_1d  = [list(p)[0] for p in mgr.find_isomorphic_path_groups(paths1d_to_plot, rag2)]
path_groups_hd  = [list(p)[0] for p in mgr.find_isomorphic_path_groups(pathshd_to_plot, rag2)]

#%%
#region #* Finally let's visualize everything
# path_groups_hd = find_isomorphic_path_groups(pathshd_to_plot, rag2)
paths_1d_trace = mgr.plot_paths(rag2, path_groups_1d, legendgroup='1', legendgrouptitle_text='1D Channels')
paths_hd_trace= mgr.plot_paths(rag2,path_groups_hd, legendgroup='2', legendgrouptitle_text='2D & 3D Channels')
#%%
fig = go.Figure(data=[*frame_trace, *paths_1d_trace, *paths_hd_trace,cell_trace, surface_trace, rag_trace]).update(layout=mgr.dark_layout)
fig.show()
#
#endregion
# %%

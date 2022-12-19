
#%% [markdown]
# # Cluster pockets example
# A short example of the sequence of steps used to identify (and visualize) the geometrically similar pockets in a nanaoporous material. 

# These simply are one set of choices for the techniques that seem to work well across the board, but feel free to experiment with different choices in clustering algorithms, and workflow parameter values
#%% [markdown]
# ## Import packages
# Import modulees after setting paths if necessary (replace the path to mgr with your own)
#%%

import sys
import os 
import numpy as np 

# #* cif io 
# from ase.io import read 

# #* for atomic radii
# from ase.data import covalent_radii

# # * cif manipulation
# import ase

# # * For manipulating the distance grid efficiently
# import dask.array as da

# #* make grid finer
# from scipy.ndimage import zoom 

# #* measure some properties ofthe regions 
# from skimage.measure import regionprops

# #* Plotting and visualization
# import plotly.graph_objects as go
# import plotly.io as pio 
# pio.renderers.default='notebook_connected'

# #* colorscale for plotting 
# import plotly.express as px 

# # * if you are using presavaed database dataframe from MOFDB
import pandas as pd 
import pickle 
import random 



#* path to mgr
path_to_mgr = os.path.abspath('C:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
# path_to_mgr = os.path.abspath('./mofography')
sys.path.append(path_to_mgr)

import mofography as mgr
#%%

# #* Read a in a cif file  as an ase.Atoms object
# path_to_cif = os.path.abspath('./mofography/examples/TestCIFs/NU-1000_S1.cif')
# # path_to_cif = './TestCIFs/NU-1000_S1.cif'
# mat_atoms = read(path_to_cif)
# mat_updated = mgr.update_with_covalent_bonds(mat_atoms) # for plotting later


# #%% [markdown]
# # ## Use a database instead (Optional)
# # instead of picking a file we can download adatabase from MOFDb and use that 
# os.chdir('C:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography/examples')

# # * Download and save the database 



# * Load the IZA database from the feather file
iza_atoms = pickle.load(open('iza_atoms_list.pkl', 'rb'))
df = pd.read_feather('iza_df.feather')

#%% [markdown]
# ## Set  params for the workflow
# Set some parameters for the workflow (this works for most cases. you might want to change the bin range if your material has a very large pores)
#%%
# #* Set some basic parameters fo the workflow 
# spacing = 0.5 # grid spacing of the distance grid in angstroms
# mask_thickness =1.5 # thickness of the grid mask, anything within this distance from the atom centers is masked out as background: region 0 
# h=0.2
# compactness = 0 # compactness for region segmentation, this should be set to 0 otherwise it gives incorrect segmentation
# zoom_ratio = 2 # zoom ratio for the distance grid, this is used to make the grid finer for better resolution of the distance map
# wall_thickness =1.0 # thickness of the wall region in angstroms, the number of wall voxels is calculated to at least this perpendicular thickness
# bins = np.linspace(0,30, 50) # the bins for the cavity histogram, used as the feature vector for clustering
# block_size = 50000 # partition size for the distance grid, this is used to avoid memory issues when calculating the distance grid
# a, b, c, =1,1,1 # optional, to make supercells of the material along each direction
# distance_threshold = 0.25 # distance threshold for agglomerative clustering, 
# #%% [markdown]

# ## Execute Workflow
# # The workflow is in the next cell, the purpose of each line is explained in the comments

# #%%

# # # * Pick a random material from the database
# chosen = random.choice(np.arange(len(iza_atoms)))
# # chosen = 108
# # chosen = 107
# print(chosen, df['name'][chosen])
# #%%
# mat_atoms = iza_atoms[chosen]
# mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

# # * compute the distance grid, locate the local maxima and make a list of their size (radii)
# radii = covalent_radii[mat_atoms.get_atomic_numbers()] # assign covalent radii
# grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing =spacing, block_size=block_size) # compute the distance grid (dask array) of given spacing and partition size
# mat_atoms = ase.build.make_supercell(mat_atoms, [[a,0,0],[0,b,0],[0,0,c]]) # (optional) make a supercell of the material
# grid = da.from_array(np.tile(grid.compute(), [a,b,c])) #
# grid = da.from_array(zoom(grid.compute(), zoom_ratio,mode='grid-wrap', order=1 ), chunks=(40,40,40)) # (optional) make the distance grid finer by zooming in, here we are aamking it twice as fine
# regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=mask_thickness, compactness=compactness, h=h) # # compute the region segementaion and array indices of the maxima corresponding to each region
# maxima_radii = grid.compute()[tuple(maxima.T)] # get teh distance values (size) of the maxima

# #* get the connectivity of the pockets, each pbc group is a list of regions that are split up by the PBC, that come together to form a complete cavity.
# # * return_conn = True returns the list of pbc groups, return_conn = False returns the a 3D array of regions where region index is replaced by the pbc group index
# pbc_groups = mgr.apply_pbc(regions, maxima, maxima_radii,mat_atoms, wall_thickness=wall_thickness, return_conn=True, minimum_overlap_fraction=0.5)  # this returns a list of pbc groups
# regions_pbc = mgr.apply_pbc(regions, maxima, maxima_radii,mat_atoms, wall_thickness=wall_thickness, return_conn=False, minimum_overlap_fraction=0.5) # this returns a 3D array 
# #%%
# # * construct the region adjacency graph, a.k.a the pore graph form the regions, grid, maxima indices and the ase atoms object
# rag2 = mgr.get_rag_from_regions_and_grid(regions, grid.compute(),mat_atoms,maxima)


# # * in this part we will be adding properties to the region adjacency graph
# rag2 = mgr.add_pbc_to_rag(rag2, pbc_groups)
# props = regionprops(regions) # * get the region properties from skimage.measure.regionprops
# rag2 = mgr.add_psd_to_rag( rag2,region_props=props,  bins = bins, ase_atoms = mat_atoms, shape =regions.shape) # *  add the cavity histogram without pbc to the RAG
# rag2 = mgr.add_psd_to_rag_with_pbc(rag = rag2, props = props, pbc_groups= pbc_groups, bins = bins, ase_atoms = mat_atoms, shape =regions.shape) # * add the cavity histogram with pbc to the RAG
# # rag2= mgr.add_clustering_to_rag(rag2) # * add the clustering labels to the RAG without pbc 

# #* add energy histogram to the RAG
# forcefield_mixed = mgr.mix_lorentz_berthelot(mat_atoms)
# # rag2 = mgr.add_vdw_hist_to_rag_with_pbc(rag2, props, pbc_groups, mat_atoms, forcefield_mixed, regions.shape)


# # rag2 = mgr.add_vdw_hist_to_rag(rag2, props, mat_atoms, forcefield_mixed, regions.shape, density=True)
# rag2, cluster_obj  = mgr.add_clustering_to_rag_pca(rag2, no_pbc=True, kwargs_cluster=dict(distance_threshold=distance_threshold, kernel = 'linear'), return_cluster_obj=True, kwargs_pca=dict(kernel='linear')) # * do clustering of the regions with PCA and default settings (clustering kernel='linear', linkage='single', distance='euclidean', pca kernel ='cosine')

# # * get the 3D array of cluster labels
# regions_cluster= mgr.clustering_grid_from_regions(regions, rag2)

# #%% [markdown]
# # ## Visualize the results
# # In this section we will visualize the results of the workflow using the built in functions 
# # these built in functions return plotly traces, which are then collected into a plotly figure object and plotted
# # either in the notebook (renderer='notebook_connected') or in a browser (renderer='browser')

# # ### We will visualize the following
# # * The distance isosurface at 1.75 A (needs to be greater than the mask thickness)
# #    * The distance isosurface colored by the region labels before applying pbc
# #    * The distance isosurface colored by the region labels after applying pbc
# #    * The distance isosurface colored by the cluster labels
# #    * The framework atoms and bonds as a mesh 
# #    * The unit cell as a box 
# #    * The local maxima as points colored by their size (radius)


# #%%
# #* set colorscale for plotting, make things consistent across the two figures
# cscale = px.colors.qualitative.Pastel

# # * Plot the distance isosurface at 1.75 A colored in different ways
# surface_trace=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions, cscale = cscale , stylize=True, name = 'Raw segmentation', hover_label='Region :')
# surface_trace_pbc=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions_pbc-1, cscale = cscale , stylize=True, name='PBC groups', hover_label='PBC Group :')
# surface_trace_cluster=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions_cluster, cscale = cscale , stylize=True, name ='Cavity types', hover_label='Cavity type :')

# # * plot the framework atoms and bonds as a mesh
# frame_trace = mgr.plot_framework_as_mesh(mat_updated, show_rings=True, stylize=True)
# # surface_trace=mgr.plot_surface_dgrid(grid.compute(), 1.75, select=list(pbc_groups[0]),ase_atoms=mat_updated,region_labels=regions, stylize=True)

# # * plot the unit cell as a box
# cell_trace=mgr.plot_cell(mat_updated)

# # * plot the local maxima as points colored by their size (radius)
# maxima_trace = mgr.plot_maxima(maxima = maxima, dgrid = grid.compute(),ase_atoms =mat_updated)
# #%%
# # * collect the traces into a figure object
# fig = go.Figure(data=[surface_trace, surface_trace_pbc, surface_trace_cluster, *frame_trace, cell_trace, maxima_trace])
# fig.update(layout=mgr.white_layout)
# fig.update_layout(scene=dict(aspectmode='data'))
# fig.show(renderer='notebook_connected')

# #%% [markdown]
# # ## Visualize the cavity histogram feature for each pbc group ( for each cavity)
# # It is ogten useful to visualize the cavity histogram for each cavity, this can be done using the built in function `plot_psds' 
# # the main purpose is ti see if the cavity histograms do a good job of represeneting the cavities for clustering, 
# # this is usefule for debugging the workflow when the clusering does not work as expected

# # %%
# # * plot the cavity histogram for each pbc group
# csd_traces = mgr.plot_psds(rag2, pbc_groups, cscale=cscale)
# fig_csd = go.Figure(data=csd_traces)
# fig_csd.update_layout(xaxis_title='Distance to centroid (\u212b)', yaxis_title='Normalized histogram' )
# fig_csd.show(renderer='notebook_connected')
# #%%
# pbc_groups
# # %%
# energy_traces= mgr.plot_vdw_hists(rag2, pbc_groups, cscale=cscale)
# fig_ehist= go.Figure(data=energy_traces)
# fig_ehist.update_layout(xaxis_title='Energy histogram (\u212b)', yaxis_title='Normalized histogram' )
# fig_ehist.show(renderer='notebook_connected')
# %%
# rag2 = mgr.add_vdw_hist_to_rag(rag2, props, mat_atoms, forcefield_mixed, regions.shape)









# %%
# def quick_cluster(mat_atoms,
# spacing = 0.5,
# mask_thickness =1.5,
# h=0.2,
# compactness = 1E-4,
# zoom_ratio = 2,
# wall_thickness =1.0,
# bins = np.linspace(0,30, 50),
# block_size = 50000 ,a=1, b=1, c=1,distance_threshold = 0.25, minimum_overlap_fraction=0.5, show_figures=True, renderer='browser', list_return_strings=[]): # distance threshold for agglomerative clustering 
#     # # # * Pick a random material from the database
#     # if chosen is None:
#     #     chosen = random.choice(np.arange(len(iza_atoms)))
#     # # chosen = 108
#     # # chosen = 107
#     # print(chosen, df['name'][chosen])
#     # mat_atoms = iza_atoms[chosen]
#     # from ase.io import read
#     # mat_atoms = read('../examples/TestCIFs/NU-1000_S1.cif')

#     mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

#     # * compute the distance grid, locate the local maxima and make a list of their size (radii)
#     radii = covalent_radii[mat_atoms.get_atomic_numbers()] # assign covalent radii
#     grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing =spacing, block_size=block_size) # compute the distance grid (dask array) of given spacing and partition size
#     mat_atoms = ase.build.make_supercell(mat_atoms, [[a,0,0],[0,b,0],[0,0,c]]) # (optional) make a supercell of the material
#     grid = da.from_array(np.tile(grid.compute(), [a,b,c])) #
#     grid = da.from_array(zoom(grid.compute(), zoom_ratio,mode='grid-wrap', order=1, grid_mode=False), chunks=(40,40,40)) # (optional) make the distance grid finer by zooming in, here we are aamking it twice as fine
#     regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=mask_thickness, compactness=compactness, h=h) # # compute the region segementaion and array indices of the maxima corresponding to each region
#     maxima_radii = grid.compute()[tuple(maxima.T)] # get teh distance values (size) of the maxima

#     #* get the connectivity of the pockets, each pbc group is a list of regions that are split up by the PBC, that come together to form a complete cavity.
#     # * return_conn = True returns the list of pbc groups, return_conn = False returns the a 3D array of regions where region index is replaced by the pbc group index
#     pbc_groups = mgr.apply_pbc(regions, maxima, maxima_radii,mat_atoms, wall_thickness=wall_thickness, return_conn=True, minimum_overlap_fraction=minimum_overlap_fraction)  # this returns a list of pbc groups
#     regions_pbc = mgr.apply_pbc(regions, maxima, maxima_radii,mat_atoms, wall_thickness=wall_thickness, return_conn=False, minimum_overlap_fraction=minimum_overlap_fraction) # this returns a 3D array 
#     # * construct the region adjacency graph, a.k.a the pore graph form the regions, grid, maxima indices and the ase atoms object
#     rag2 = mgr.get_rag_from_regions_and_grid(regions, grid.compute(),mat_atoms,maxima)


#     # * in this part we will be adding properties to the region adjacency graph
#     rag2 = mgr.add_pbc_to_rag(rag2, pbc_groups)
#     props = regionprops(regions) # * get the region properties from skimage.measure.regionprops
#     rag2 = mgr.add_psd_to_rag( rag2,region_props=props,  bins = bins, ase_atoms = mat_atoms, shape =regions.shape) # *  add the cavity histogram without pbc to the RAG
#     rag2 = mgr.add_psd_to_rag_with_pbc(rag = rag2, props = props, pbc_groups= pbc_groups, bins = bins, ase_atoms = mat_atoms, shape =regions.shape) # * add the cavity histogram with pbc to the RAG
#     # rag2= mgr.add_clustering_to_rag(rag2) # * add the clustering labels to the RAG without pbc 

#     #* add energy histogram to the RAG
#     forcefield_mixed = mgr.mix_lorentz_berthelot(mat_atoms)
#     # rag2 = mgr.add_vdw_hist_to_rag_with_pbc(rag2, props, pbc_groups, mat_atoms, forcefield_mixed, regions.shape)


#     rag2, cluster_obj  = mgr.add_clustering_to_rag_pca(rag2, no_pbc=False, kwargs_cluster=dict(distance_threshold=distance_threshold, kernel = 'linear'), return_cluster_obj=True, kwargs_pca=dict(kernel='cosine')) # * do clustering of the regions with PCA and default settings (clustering kernel='linear', linkage='single', distance='euclidean', pca kernel ='cosine')

#     #---------------------------------------------------#
#     # #* add energy histogram to the RAG
#     # props = regionprops(regions)
#     # forcefield_mixed = mgr.mix_lorentz_berthelot(mat_atoms)
#     # rag2 = mgr.add_vdw_hist_to_rag(rag2, props, mat_atoms, forcefield_mixed, regions.shape, use_pbc=True,pbc_groups=pbc_groups,density=False)
    
#     # #* add clustering to rag using energy with pbc
#     # rag2, cluster_obj  = add_clustering_to_rag_pca2(rag2, no_pbc=False, kwargs_cluster=dict(distance_threshold=distance_threshold, kernel = 'linear'), return_cluster_obj=True, kwargs_pca=dict(kernel='cosine')) # * do clustering of the regions with PCA and default settings (clustering kernel='linear', linkage='single', distance='euclidean', pca kernel ='cosine')
#     #-------------------------------------------------#
#     # * get the 3D array of cluster labels
#     regions_cluster= mgr.clustering_grid_from_regions(regions, rag2)

#     if show_figures:
#         #* set colorscale for plotting, make things consistent across the two figures
#         cscale = px.colors.qualitative.Plotly

#         # * Plot the distance isosurface at 1.75 A colored in different ways
#         surface_trace=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions, cscale = cscale , stylize=True, name = 'Raw segmentation', hover_label='Region :')
#         surface_trace_pbc=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions_pbc-1, cscale = cscale , stylize=True, name='PBC groups', hover_label='PBC Group :')
#         surface_trace_cluster=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions_cluster-1, cscale = cscale , stylize=True, name ='Cavity types', hover_label='Cavity type :')

#         # * plot the framework atoms and bonds as a mesh
#         frame_trace = mgr.plot_framework_as_mesh(mat_updated, show_rings=True, stylize=True)
#         # surface_trace=mgr.plot_surface_dgrid(grid.compute(), 1.75, select=list(pbc_groups[0]),ase_atoms=mat_updated,region_labels=regions, stylize=True)

#         # * plot the unit cell as a box
#         cell_trace=mgr.plot_cell(mat_updated)

#         # * plot the local maxima as points colored by their size (radius)
#         maxima_trace = mgr.plot_maxima(maxima = maxima, dgrid = grid.compute(),ase_atoms =mat_updated)
        
#         #* plot the windows of the material 
#         # connections = mgr.get_connections_for_rag(rag2, region_labels=regions, ase_atoms = mat_atoms, maxima=maxima, dgrid=grid.compute())
#         # windows_trace = mgr.plot_windows(connections)
#         # * collect the traces into a figure object
#         fig_full_visual = go.Figure(data=[surface_trace, surface_trace_pbc, surface_trace_cluster, *frame_trace, cell_trace, maxima_trace])
#         # fig_full_visual = go.Figure(data=[surface_trace, surface_trace_pbc, surface_trace_cluster, *frame_trace, cell_trace, maxima_trace, windows_trace])
#         fig_full_visual.update(layout=mgr.white_layout)
#         fig_full_visual.update_layout(scene=dict(aspectmode='data'))
#         fig_full_visual.show(renderer=renderer)


#         # * plot the cavity histogram for each pbc group
#         csd_traces = mgr.plot_psds(rag2, pbc_groups, cscale=cscale)
#         fig_csd = go.Figure(data=csd_traces)
#         fig_csd.update_layout(xaxis_title='Distance to centroid (\u212b)', yaxis_title='Normalized histogram')
#         fig_csd.show(renderer=renderer)


#     # # # * compute and plot the energy histograms for each pbc group 
#     # rag2 = mgr.add_vdw_hist_to_rag(rag2, props, mat_atoms, forcefield_mixed, regions.shape,use_pbc=True,  pbc_groups=pbc_groups, density=True)
#     # ehist_traces = mgr.plot_vdw_hists(rag2, pbc_groups, cscale=cscale)
#     # fig_ehist = go.Figure(data=ehist_traces)
#     # fig_ehist.update_layout(xaxis_title='Energy (K)', yaxis_title='Normalized histogram')
#     # fig_ehist.show(renderer=renderer)
#     if len(list_return_strings)>0:
#         print('Some return objects specified, possible options to pass in list_return_strings are:')
#         list_return_objects =[vars()[name] for name in list_return_strings]

#         return list_return_objects
#     else:
#         print('No return objects specified, possible options to pass in list_return_strings are:')
#         print(vars().keys())    


#%%

chosen = random.choice(np.arange(len(iza_atoms)))
chosen=109
# nice ones are 12, 87, 139, 61
#bad ones are 97, 109,174,199,134,66,73
mat_atoms = iza_atoms[chosen]
print('chosen material index is {0}, name is {1}'.format( chosen,df['name'][chosen]))
use_energy_for_clustering=True
sample_fraction = 1
rag2=mgr.quick_cluster(mat_atoms,spacing=0.5,compactness=0.001, a=1, b=1, c=1, zoom_ratio=2, sample_fraction=sample_fraction, distance_threshold=0.4, use_energy_for_clustering=use_energy_for_clustering, density=True, return_these_vars=['rag2'])[0]
#%%
props = regionprops(regions)
#endregion
rag = mgr.add_vdw_hist_to_rag(rag2, props, mat_atoms, forcefield_mixed, regions.shape,use_pbc=True,  pbc_groups=pbc_groups, density=True)
#%%
# * instead of pulling a random mOF from IZA make a random one using pormake

path_to_poremake= os.path.abspath('C:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/PORMAKE/PORMAKE-master')
sys.path.append(path_to_poremake)
import pormake as pm 
# %%


def make_random_mof(db,seed=None):
    
    import pormake as pm 
    import numpy as np
    if seed is not None:
        np.random.seed(seed)
    builder = pm.Builder()
    locator = pm.Locator()
    topo_list = db.topo_list
    bblist = np.array(db.bb_list)
    edges = np.array([bb for bb in bblist if 'E' in bb])
    nodes = np.array([bb for bb in bblist if 'N' in bb])
    node_conns = [db.get_bb(n).n_connection_points for n in nodes]

    top = db.get_topo(np.random.choice(topo_list))
    print(f'Chosen topology is {top.name}')

    accepted_edge_bbs = [db.get_bb(e) for e in np.random.choice(edges, size=top.n_edge_types, replace=True)]
    print(f'Chosen edges are {[e.name for e in accepted_edge_bbs]}')

    accepted_node_bbs = []
    top_node_options = [nodes[np.where([c==nc for c in node_conns])] for nc in top.unique_cn]
    for i,nopts in enumerate(top_node_options):

        ni = db.get_bb(np.random.choice(nopts))
        rmse = locator.calculate_rmsd(top.unique_local_structures[i], ni)
        # print(f'rmse for {ni.name} is {rmse}')
        while (rmse > 0.3):
            ni = db.get_bb(random.choice(nopts))
            rmse = locator.calculate_rmsd(top.unique_local_structures[i], ni)
            # print(f'rmse for {ni.name} is {rmse}')
        accepted_node_bbs.append(ni)
    print(f'Accepted node {ni.name} with rmse {rmse}')
    edge_types = [tuple(e) for e in top.unique_edge_types]
    edge_bb_dict= dict(zip(edge_types, accepted_edge_bbs))
    node_bb_dict = dict(enumerate(accepted_node_bbs))
    mat_pormake = builder.build_by_type(top, node_bb_dict, edge_bb_dict)
    # mat_atoms = mat_pormake
    return mat_pormake

# %%
# db = pm.Database()
# mat_atoms = make_random_mof(db,seed=np.random.randint(0,100000)).atoms
# mat_atoms
#%%
# rag2=mgr.quick_cluster(mat_atoms,compactness=0,spacing=1,mask_thickness=1.5, isoval=2, zoom_ratio=1, return_these_vars=[])#,spacing =0.5, zoom_ratio=2, distance_threshold=0.20, mop=0.5, wall_thickness=1)
#%%
from ase.data import covalent_radii
import ase
import dask.array as da
from scipy.ndimage import zoom 
from skimage.measure import label, regionprops
#* Set some basic parameters fo the workflow 
spacing = 0.5 # grid spacing of the distance grid in angstroms
mask_thickness =1.5 # thickness of the grid mask, anything within this distance from the atom centers is masked out as background: region 0 
h=0.2
compactness = 1E-3 # compactness for region segmentation, this should be set to 0 otherwise it gives incorrect segmentation
zoom_ratio = 2 # zoom ratio for the distance grid, this is used to make the grid finer for better resolution of the distance map
wall_thickness =1.0 # thickness of the wall region in angstroms, the number of wall voxels is calculated to at least this perpendicular thickness
bins = np.linspace(0,30, 50) # the bins for the cavity histogram, used as the feature vector for clustering
block_size = 50000 # partition size for the distance grid, this is used to avoid memory issues when calculating the distance grid
a, b, c, =1,1,1 # optional, to make supercells of the material along each direction
distance_threshold = 0.25 # distance threshold for agglomerative clustering, 
forcefield_mixed = mgr.mix_lorentz_berthelot(mat_atoms)
use_pbc = True
bins_energy = np.linspace(-3000, 100, 25)
radii = covalent_radii[mat_atoms.get_atomic_numbers()] # assign covalent radii
grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing =spacing, block_size=block_size) # compute the distance grid (dask array) of given spacing and partition size
mat_atoms = ase.build.make_supercell(mat_atoms, [[a,0,0],[0,b,0],[0,0,c]]) # (optional) make a supercell of the material
grid = da.from_array(np.tile(grid.compute(), [a,b,c])) #
grid = da.from_array(zoom(grid.compute(), zoom_ratio,mode='grid-wrap', order=1 ), chunks=(40,40,40)) # (optional) make the distance grid finer by zooming in, here we are aamking it twice as fine
regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=mask_thickness, compactness=compactness, h=h) # # compute the region segementaion and array indices of the maxima corresponding to each region
maxima_radii = grid.compute()[tuple(maxima.T)] # get teh distance values (size) of the maxima

#* get the connectivity of the pockets, each pbc group is a list of regions that are split up by the PBC, that come together to form a complete cavity.
# * return_conn = True returns the list of pbc groups, return_conn = False returns the a 3D array of regions where region index is replaced by the pbc group index
pbc_groups = mgr.apply_pbc(regions, maxima, maxima_radii,mat_atoms, wall_thickness=wall_thickness, return_conn=True, minimum_overlap_fraction=0.5)  # this returns a list of pbc groups
#%% [markdown]
# * Let's get th energies for a couple of the pbc groups and see if we can tell they are similar


#%%



#%%
import dask.bag as db 
pbc_groups_lists = [list(g) for g in pbc_groups]
props = regionprops(regions)
shape = regions.shape
cutoff = 12.8
ase_atoms = mat_atoms
energies_list = [mgr.get_energies_for_group(list(g), props, ase_atoms, forcefield_mixed, shape, cutoff=cutoff,sample_fraction=0.5) for g in pbc_groups]

#%%
histograms_bins= [np.histogram(e, bins=bins_energy, density=False) for e in energies_list]
ehistograms = [hb[0] for hb in histograms_bins]
ebincenters = [(hb[1][:-1]+hb[1][1:])/2 for hb in histograms_bins]

#%% 
# * test the samples using ks2samp
from scipy.stats import ks_2samp
import itertools
import networkx as nx 
ks_graph = nx.Graph()
ks_graph.add_nodes_from(range(len(pbc_groups)))
for i,j in itertools.combinations(range(len(pbc_groups)), 2):
    ks = ks_2samp(energies_list[i], energies_list[j])
    if ks[1] < 0.05:
        ks_graph.add_edge(i,j)
conn = list(nx.connected_components(ks_graph))
energy_groups = [set(itertools.chain.from_iterable(np.array(pbc_groups)[list(c)])) for c in conn]

# for i, c in enumerate(conn):

print(f'the number of types based on energy is {len(energy_groups)}')
print(f'the typs of cavities based on energy are {energy_groups} ')


#     print(f'the ks test for pbc_groups {pbc_groups[i]},{pbc_groups[j]} is {ks}')
#     if ks[1] < 0.05:
#         print(f'{pbc_groups[i]} and {pbc_groups[j]} are from the same dist. 95%')
#     else:
#         print(f'{pbc_groups[i]} and {pbc_groups[j]} are from different dist. 95%')
# # print(f'the ks test for pair 2 is {ks_2samp(energies_list[0], energies_list[2])}')
# print(f'the ks test for pair 3 is {ks_2samp(energies_list[1], energies_list[2])}')


#%%
from dask.diagnostics import ProgressBar
with ProgressBar():
    energies_computed = energies_bag.compute(scheduler='threads')


#%%
import time 
start = time.time()
energies_list = [mgr.get_energies_for_group(list(g), props, ase_atoms, forcefield_mixed, shape, cutoff=cutoff) for g in pbc_groups]
print(f'it took {time.time()-start} seconds')
#%%
if use_pbc:
    if pbc_groups is None:
        raise ValueError('pbc_groups must be provided if use_pbc=True')
    else:
        energies_list = [get_energies_for_group(list(g), props, ase_atoms, forcefield_mixed, shape, cutoff=cutoff) for g in pbc_groups]
        histograms_bins= [np.histogram(e, bins=bins_energy, density=density) for e in energies_list]
        ehistograms = [hb[0] for hb in histograms_bins]
        ebincenters = [(hb[1][:-1]+hb[1][1:])/2 for hb in histograms_bins]
        for i, node in enumerate(rag.nodes()):
            rag.nodes[node]['vdw_hist_y_pbc']=ehistograms[rag.nodes[node]['pbc_group']]
            rag.nodes[node]['vdw_hist_x_pbc']=ebincenters[rag.nodes[node]['pbc_group']]

else:
    energies_list = [get_energies_for_group([node], props, ase_atoms, forcefield_mixed, shape, cutoff=cutoff) for node in rag.nodes()]
    histograms_bins= [np.histogram(e, bins=bins_energy, density=density) for e in energies_list]
    ehistograms = [hb[0] for hb in histograms_bins]
    ebincenters = [(hb[1][:-1]+hb[1][1:])/2 for hb in histograms_bins]


    for i, node in enumerate(rag.nodes()):
        rag.nodes[node]['vdw_hist_y']=ehistograms[i]
        rag.nodes[node]['vdw_hist_x']=ebincenters[i]
# nx.set_node_attributes(rag, dict(zip(rag.nodes,ehists_data)), 'vdw_hist_y_pbc')
# nx.set_node_attributes(rag, dict(zip(rag.nodes,ebins_data)), 'vdw_hist_x_pbc')
if add_energies:
    for i,node in rag.nodes():
        rag.nodes[node]['energies']=energies_list[i]
        
        # energies_nodes = [energies_list[c] for c in cindices]
        # nx.set_node_attributes(rag, dict(zip(rag.nodes, energies_nodes)), 'vdw_energies_pbc')
    
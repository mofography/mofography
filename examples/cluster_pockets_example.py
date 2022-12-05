
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

#* cif io 
from ase.io import read 

#* for atomic radii
from ase.data import covalent_radii

# * cif manipulation
import ase

# * For manipulating the distance grid efficiently
import dask.array as da

#* make grid finer
from scipy.ndimage import zoom 

#* measure some properties ofthe regions 
from skimage.measure import regionprops

#* Plotting and visualization
import plotly.graph_objects as go
import plotly.io as pio 
pio.renderers.default='notebook_connected'

#* colorscale for plotting 
import plotly.express as px 

# * if you are using presavaed database dataframe from MOFDB
import pandas as pd 
import pickle 
import random 



#* path to mgr
path_to_mgr = os.path.abspath('C:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
sys.path.append(path_to_mgr)

import mofography as mgr
#%%

#* Read a in a cif file  as an ase.Atoms object
mat_atoms = read('./TestCIFs/NU-1000_S1.cif')
mat_updated = mgr.update_with_covalent_bonds(mat_atoms) # for plotting later


#%% [markdown]
# ## Use a database instead (Optional)
# instead of picking a file we can download adatabase from MOFDb and use that 
os.chdir('C:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography/examples')

# * Download and save the database 



# * Load the IZA database from the feather file
iza_atoms = pickle.load(open('iza_atoms_list.pkl', 'rb'))
df = pd.read_feather('iza_df.feather')

#%% [markdown]
# ## Set  params for the workflow
# Set some parameters for the workflow (this works for most cases. you might want to change the bin range if your material has a very large pores)
#%%
#* Set some basic parameters fo the workflow 
spacing = 0.5 # grid spacing of the distance grid in angstroms
mask_thickness =1.5 # thickness of the grid mask, anything within this distance from the atom centers is masked out as background: region 0 
h=0.2
compactness = 0 # compactness for region segmentation, this should be set to 0 otherwise it gives incorrect segmentation
zoom_ratio = 2 # zoom ratio for the distance grid, this is used to make the grid finer for better resolution of the distance map
wall_thickness =1.0 # thickness of the wall region in angstroms, the number of wall voxels is calculated to at least this perpendicular thickness
bins = np.linspace(0,30, 50) # the bins for the cavity histogram, used as the feature vector for clustering
block_size = 50000 # partition size for the distance grid, this is used to avoid memory issues when calculating the distance grid
a, b, c, =1,1,1 # optional, to make supercells of the material along each direction
distance_threshold = 0.25 # distance threshold for agglomerative clustering, 
#%% [markdown]

## Execute Workflow
# The workflow is in the next cell, the purpose of each line is explained in the comments

#%%

# * Pick a random material from the database
chosen = random.choice(np.arange(len(iza_atoms)))
chosen = 108
print(chosen, df['name'][chosen])
#%%
mat_atoms = iza_atoms[chosen]
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

# * compute the distance grid, locate the local maxima and make a list of their size (radii)
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
regions_pbc = mgr.apply_pbc(regions, maxima, maxima_radii,mat_atoms, wall_thickness=wall_thickness, return_conn=False, minimum_overlap_fraction=0.5) # this returns a 3D array 
#%%
# * construct the region adjacency graph, a.k.a the pore graph form the regions, grid, maxima indices and the ase atoms object
rag2 = mgr.get_rag_from_regions_and_grid(regions, grid.compute(),mat_atoms,maxima)


# * in this part we will be adding properties to the region adjacency graph
rag2 = mgr.add_pbc_to_rag(rag2, pbc_groups)
props = regionprops(regions) # * get the region properties from skimage.measure.regionprops
rag2 = mgr.add_psd_to_rag( rag2,region_props=props,  bins = bins, ase_atoms = mat_atoms, shape =regions.shape) # *  add the cavity histogram without pbc to the RAG
rag2 = mgr.add_psd_to_rag_with_pbc(rag = rag2, props = props, pbc_groups= pbc_groups, bins = bins, ase_atoms = mat_atoms, shape =regions.shape) # * add the cavity histogram with pbc to the RAG
# rag2= mgr.add_clustering_to_rag(rag2) # * add the clustering labels to the RAG without pbc 

#* add energy histogram to the RAG
forcefield_mixed = mgr.mix_lorentz_berthelot(mat_atoms)
# rag2 = mgr.add_vdw_hist_to_rag_with_pbc(rag2, props, pbc_groups, mat_atoms, forcefield_mixed, regions.shape)


# rag2 = mgr.add_vdw_hist_to_rag(rag2, props, mat_atoms, forcefield_mixed, regions.shape, density=True)
rag2, cluster_obj  = mgr.add_clustering_to_rag_pca(rag2, no_pbc=True, kwargs_cluster=dict(distance_threshold=distance_threshold, kernel = 'linear'), return_cluster_obj=True, kwargs_pca=dict(kernel='linear')) # * do clustering of the regions with PCA and default settings (clustering kernel='linear', linkage='single', distance='euclidean', pca kernel ='cosine')

# * get the 3D array of cluster labels
regions_cluster= mgr.clustering_grid_from_regions(regions, rag2)

#%% [markdown]
# ## Visualize the results
# In this section we will visualize the results of the workflow using the built in functions 
# these built in functions return plotly traces, which are then collected into a plotly figure object and plotted
# either in the notebook (renderer='notebook_connected') or in a browser (renderer='browser')

# ### We will visualize the following
# * The distance isosurface at 1.75 A (needs to be greater than the mask thickness)
#    * The distance isosurface colored by the region labels before applying pbc
#    * The distance isosurface colored by the region labels after applying pbc
#    * The distance isosurface colored by the cluster labels
#    * The framework atoms and bonds as a mesh 
#    * The unit cell as a box 
#    * The local maxima as points colored by their size (radius)


#%%
#* set colorscale for plotting, make things consistent across the two figures
cscale = px.colors.qualitative.Pastel

# * Plot the distance isosurface at 1.75 A colored in different ways
surface_trace=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions, cscale = cscale , stylize=True, name = 'Raw segmentation', hover_label='Region :')
surface_trace_pbc=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions_pbc-1, cscale = cscale , stylize=True, name='PBC groups', hover_label='PBC Group :')
surface_trace_cluster=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions_cluster, cscale = cscale , stylize=True, name ='Cavity types', hover_label='Cavity type :')

# * plot the framework atoms and bonds as a mesh
frame_trace = mgr.plot_framework_as_mesh(mat_updated, show_rings=True, stylize=True)
# surface_trace=mgr.plot_surface_dgrid(grid.compute(), 1.75, select=list(pbc_groups[0]),ase_atoms=mat_updated,region_labels=regions, stylize=True)

# * plot the unit cell as a box
cell_trace=mgr.plot_cell(mat_updated)

# * plot the local maxima as points colored by their size (radius)
maxima_trace = mgr.plot_maxima(maxima = maxima, dgrid = grid.compute(),ase_atoms =mat_updated)
#%%
# * collect the traces into a figure object
fig = go.Figure(data=[surface_trace, surface_trace_pbc, surface_trace_cluster, *frame_trace, cell_trace, maxima_trace])
fig.update(layout=mgr.white_layout)
fig.update_layout(scene=dict(aspectmode='data'))
fig.show(renderer='browser')

#%% [markdown]
# ## Visualize the cavity histogram feature for each pbc group ( for each cavity)
# It is ogten useful to visualize the cavity histogram for each cavity, this can be done using the built in function `plot_psds' 
# the main purpose is ti see if the cavity histograms do a good job of represeneting the cavities for clustering, 
# this is usefule for debugging the workflow when the clusering does not work as expected

# %%
# * plot the cavity histogram for each pbc group
csd_traces = mgr.plot_psds(rag2, pbc_groups, cscale=cscale)
fig_csd = go.Figure(data=csd_traces)
fig_csd.update_layout(xaxis_title='Distance to centroid (\u212b)', yaxis_title='Normalized histogram' )
fig_csd.show(renderer='browser')
#%%
pbc_groups
# # %%
# energy_traces= mgr.plot_vdw_hists(rag2, pbc_groups, cscale=cscale)
# fig_ehist= go.Figure(data=energy_traces)
# fig_ehist.update_layout(xaxis_title='Energy histogram (\u212b)', yaxis_title='Normalized histogram' )
# fig_ehist.show(renderer='browser')
# %%
# rag2 = mgr.add_vdw_hist_to_rag(rag2, props, mat_atoms, forcefield_mixed, regions.shape)
# %%
import networkx as nx 
nx.get_node_attributes(rag2, 'cluster')
# %%
pbc_groups
# %%

rag2 = mgr.add_vdw_hist_to_rag(rag2, props, mat_atoms, forcefield_mixed, regions.shape, use_pbc=True, pbc_groups=pbc_groups,density=False)
energy_traces= mgr.plot_vdw_hists(rag2, pbc_groups, no_pbc=False, cscale=cscale)
fig_ehist= go.Figure(data=energy_traces)
fig_ehist.update_layout(xaxis_title='Energy histogram (\u212b)', yaxis_title='Normalized histogram' )
fig_ehist.show(renderer='browser')
# %%
np.vstack(list(nx.get_node_attributes(rag2, 'vdw_hist_x_pbc').values()))
# %%

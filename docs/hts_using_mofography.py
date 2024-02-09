#%% [markdown]
# # High-Throughput Processing using Mofography and Dask
# In thie notebook, we'll load up a dataset of MOFs, such as the IZA database, create ASE atoms objects from it, create a workflow for a single materials and finally extend and parallelize it over the whole dataset using Dask. This notebook uses a 'local' cluster, meaning your personal laptop or desktop with a single CPU, however, the same code can be used to scale up to a cluster with many CPUs or GPUs. Please refer the Dask documentation for more information on how to set up a cluster here: https://docs.dask.org
#%% [markdown]
# ## Imports and Paths

# Let's import the necessary libraries and set some paths
#%%
#* General imports
import os, sys
import numpy as np 
from pathlib import Path 
import pandas as pd
from tqdm import tqdm

#* Add the path to mofography to the system path
path_to_mgr = Path('C:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
sys.path.append(str(path_to_mgr))


import mofography as mgr

#* Needed for cheminformatics if any
from ase.io import read, write
from ase.build import make_supercell 
from ase.geometry import complete_cell


#* for plotting
import plotly.graph_objects as go

#* Needed for the workflow
import dask
from dask import delayed
from dask.distributed import Client, progress, LocalCluster
# %% [markdown]
# ## Load the dataset 

# I have the dataset saved up, if you don't have it, check out hte getting data from MOFDB tutuorial to create and save the dataframe.

# %%
iza_df = pd.read_feather('C:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography/examples/iza_df.feather')
# %% [markdown]
# the dataset looks like this:

# %% 
print(iza_df.columns)
iza_df.head(10)
# %%
## Create ASE atoms objects from the dataframe
# %%
iza_atoms = mgr.atoms_from_mofdb_df(iza_df)
# %%
## Create a workflow for a single material
# Let's pick one random material and define our usual workflow where we create pore graphs from the atoms object. If this looks unfamiliar, please refer to the previous tutorials.
# %%
mat_atoms = iza_atoms[3] # you can try other materials also
# %%
# compute bonds
mat_updated  = mgr.update_with_covalent_bonds(mat_atoms)
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=0.5) #* compute distance grid
shape = dgrid.shape #* we'll need this for later
regions, maxima = mgr.regions_from_dgrid(dgrid, mask_thickness=1.0) #* compute regions
maxima_radii = dgrid[tuple(maxima.T)] #* grab the radii of the maxima from dgrid
regions_pbc  = mgr.apply_pbc(regions, maxima, maxima_radii, mat_atoms, return_conn=False) #* apply PBC and return as 3D array
pbc_groups = mgr.apply_pbc(regions, maxima, maxima_radii, mat_atoms, return_conn=True) #* return as connected groups
rag = mgr.get_rag_from_regions_and_grid(regions, dgrid, mat_atoms, maxima) #* create pore graph
rag = mgr.add_maxima_to_rag(rag, maxima, maxima_radii, shape, mat_atoms) #* add maxima to the graph
rag = mgr.add_pbc_to_rag(rag, pbc_groups) #* add PBC information to the graph
connections = mgr.get_connections_for_rag(rag, regions, maxima, mat_atoms, dgrid) #* compute windows
# %% #* store the information into the rag 
rag.graph['mat_atoms'] = mat_atoms
rag.graph['dgrid'] = dgrid
rag.graph['regions'] = regions
rag.graph['maxima'] = maxima
rag.graph['maxima_radii'] = maxima_radii
rag.graph['regions_pbc'] = regions_pbc
rag.graph['pbc_groups'] = pbc_groups
rag.graph['connections'] = connections

# %% #* Visualize really quickly
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)
cell_trace = mgr.plot_cell(mat_updated)
maxima_trace = mgr.plot_maxima(dgrid, maxima, mat_atoms)
surface_regions = mgr.plot_surface_dgrid(dgrid, isovalue=1.5, ase_atoms = mat_atoms, region_labels=regions, stylize=True)
rag_trace = mgr.plot_rag(rag)

fig = go.Figure(data=[*frame_trace, cell_trace, maxima_trace, surface_regions, rag_trace], layout=mgr.epic_layout)
fig.update_scenes(aspectmode='data')
fig.show()
# %% [markdown]
## Define workflow using Dask delayed
# Dask delayed essentially defines the compute graph for the workflow, without actually runnning it right away. This is very good for two reasons:
# - We can parallelize the workflow very well. Not just over materials, but, over the steps needed for each material also can be parallelized. 
# - If one the calculations fail, it won't affect the rest of the materials
# This in only a quick and dirty use of Dask without any of the nuances explained in detail. 
#
# We can have another detailed intro tutorial on dask in the future, but, the current dask documentation is very good and I highly recommend it.
#  However, there are a couple ofthings we need to keep in mind, dask delayed objects are immutable, hence, we can't work with them like tuples like how we did above. But, this is an easy fix, we just have to wrap the couple of steps into a function like this : 

# %%
def regions_maxima_radii_from_dgrid(dgrid, mask_thickness = 1.0):
    regions, maxima = mgr.regions_from_dgrid(dgrid, mask_thickness = mask_thickness)
    maxima_radii = dgrid[tuple(maxima.T )]
    return regions, maxima, maxima_radii

def make_rag(regions, dgrid, maxima, maxima_radii, mat_atoms):
    rag = mgr.get_rag_from_regions_and_grid(regions, dgrid, mat_atoms)
    shape = dgrid.shape
    rag = mgr.add_maxima_to_rag(rag, maxima, maxima_radii, shape, mat_atoms)
    rag.graph['dgrid'] = dgrid
    rag.graph['regions'] = regions
    rag.graph['maxima'] = maxima
    rag.graph['maxima_radii'] = maxima_radii
    rag.graph['mat_atoms'] = mat_atoms
    return rag

def connections_for_rag(rag):
    connections = mgr.get_connections_for_rag(rag, rag.graph['regions'], rag.graph['maxima'], rag.graph['mat_atoms'], rag.graph['dgrid'])
    rag.graph['connections'] = connections
    return rag
# %% [markdown]
# Now, we can use these functions to define the workflow using dask delayed.

# %%
delayed_rag_futures = [] 

for mat_atoms in iza_atoms:
    dgrid = delayed(mgr.dgrid_from_atoms_cpu)(mat_atoms, spacing=0.5, precompute=True)
    regions_maxima_radii = delayed(regions_maxima_radii_from_dgrid)(dgrid, mask_thickness=0.5)
    regions_pbc = delayed(mgr.apply_pbc)(regions_maxima_radii[0], regions_maxima_radii[1], regions_maxima_radii[2], mat_atoms)
    rag = delayed(make_rag)(regions_maxima_radii[0], dgrid, regions_maxima_radii[1],regions_maxima_radii[2], mat_atoms)
    rag = delayed(connections_for_rag)(rag)
    delayed_rag_futures.append(rag)
# %%
#* try running one to see if it works
delayed_rag_futures[0].compute()
# %% [markdown]
# ## Start a dask cluster locally. 
# Look at the dask documentation for more information on how to set up a cluster on a remote machine or a cloud service or to change the number of workers, processes, threads etc. 
# %%
cluster = LocalCluster()
client = Client(cluster)

print('THe client is running on : ', client.dashboard_link)
# %% [markdown]
# IF you go to the link, you'll see something like this:  https://docs.dask.org/en/stable/dashboard.html. You can see the progress of the workflow and the resources being used.
# %%
#* submit the workflow using the futures API

futures = dask.persist(*delayed_rag_futures) 
# %%
#* quick progress tracking, run this cell again for updates after any errors

progress(futures)
# %% [markdown]
# after the workflow is done, we can grab the results like this
#%%
results = client.compute(futures)
#%% [markdown]
# grab the results for finished futures. If some are still running or errored, let's set them to None
#%%
rags = [r.result() if r.status == 'finished' else None for r in results]
# %% [markdown]
# Check how many of them failed

(np.array(rags) == None).sum()
#%% [markdown] 
# Notice that some of them failed, this is OK, it is just to show how to handle errors (the updated routines don't fail as much :)). Let's just redefine the workflow with updated functions where will compute the edges and windows first and then make a graph from that, instead of using RAG functions from skimage. 
# %% [markdown]
# ## Redefine the workflow with different functions
# %% #* for single material
connections = mgr.connections_from_regions_and_dgrid(regions, dgrid, maxima, mat_atoms)
rag  = mgr.rag_from_connections(connections, maxima)

#* rest of the steps are the same

#%% [markdown]
# Now, let's redefine the workflow using dask delayed like before, notice the change in the order of the calculations 
# before: dgrid -> regions -> maxima -> regions_pbc -> rag -> connections -> rest of the steps
# after: dgrid -> regions -> maxima -> connections -> rag-> rest of the steps

#%%

def make_rag_updated(regions, dgrid, maxima,maxima_radii, mat_atoms):
    connections = mgr.connections_from_regions_and_dgrid(regions, dgrid, maxima, mat_atoms)
    rag = mgr.rag_from_connections(connections, maxima)
    shape = dgrid.shape
    rag = mgr.add_maxima_to_rag(rag, maxima, maxima_radii, shape, mat_atoms)
    rag.graph['dgrid'] = dgrid
    rag.graph['regions'] = regions
    rag.graph['maxima'] = maxima
    rag.graph['ase_atoms'] = mat_atoms
    rag.graph['connections'] = connections

    return rag

delayed_rag_futures = [] 

for mat_atoms in iza_atoms:
    dgrid = delayed(mgr.dgrid_from_atoms_cpu)(mat_atoms, spacing=0.5, precompute=True)
    # dgrid = delayed(mgr.zoom_grid)(dgrid, 2.0) #* let's zoom the grid make it finer using interpolation
    regions_maxima_radii = delayed(regions_maxima_radii_from_dgrid)(dgrid, mask_thickness=0.5)
    regions_pbc = delayed(mgr.apply_pbc)(regions_maxima_radii[0], regions_maxima_radii[1], regions_maxima_radii[2], mat_atoms)
    #* notice the change in the wrapping function here, defined above
    rag = delayed(make_rag_updated)(regions_maxima_radii[0], dgrid, regions_maxima_radii[1], regions_maxima_radii[2], mat_atoms)
    delayed_rag_futures.append(rag)
# %% #* let's try one
delayed_rag_futures[0].compute()
#%%
client = Client(LocalCluster())

print('THe client is running on : ', client.dashboard_link)
#%%
futures= dask.persist(*delayed_rag_futures)

# %%
results = client.compute(futures)
rags = [r.result() if r.status == 'finished' else None for r in results]

#%%

def connections_from_regions_and_dgrid(region_labels, dgrid, maxima, ase_atoms, wall_windows=True, distance_threshold=1.0, wall_thickness=2):

    """
    It takes a 3D array of region labels, a list of maxima, an ASE atoms object, and a distance
    threshold, and returns a list of connections between regions, and the centers of the windows between
    them
    
    :param region_labels: The labels of the regions, as returned by the skimage.measure.label function
    :param maxima: the maxima of the image
    :param ase_atoms: ASE atoms object
    :param wall_windows: If True, then the windows on the wall are also calculated, defaults to True
    (optional)
    :param distance_threshold: The distance threshold for the agglomerative clustering
    :return: A list of connections between regions.
    """

    def check_connected(regs):
            
        import numpy as np
        from icecream import ic
        ic.disable()
        from sklearn.cluster import AgglomerativeClustering
        i = regs[0]
        j = regs[1]
        ic.disable()
        ic(regs)

        check_flag = np.logical_or(np.logical_and(region_labels == i, outer.get(j)), np.logical_and(region_labels == j, outer.get(i)))
        if np.sum(check_flag) > 0:
            # print("Edges found between "+str(i)+' and '+ str(j) + ' sum: '+str(np.sum(check_flag)))
            windices = np.vstack(np.where(check_flag)).T  # * (N,3)
            window_points = np.dot(A_unit, (windices / region_labels.shape).T).T  # * (N,3)
            # ic(window_points)
            # * No clustering
            # window_centers = np.mean(window_points, axis=0)
            
            # * DBScan
            # X_scaled = StandardScaler().fit_transform(window_points)
            # db = DBSCAN(eps=0.5, min_samples=10).fit(X_scaled)
            # window_centers = np.vstack([window_points[db.labels_==i].mean(axis=0) for i in np.unique(db.labels_) if i != -1])
            
            # if len(window_points) >1:
            # * Agglomerative
            ag = AgglomerativeClustering(n_clusters=None, linkage='single',distance_threshold=minimum_window_separation, compute_full_tree=True).fit(window_points)
            #* we are taking the centroid of the clusters here as the window location
            # window_centers = np.vstack([window_points[ag.labels_==i].mean(axis=0) for i in np.unique(ag.labels_) if i != -1])
            
            # * what if we took the index of the point where the distance value is the largest in the window boundary region, that way we don't have to worry bout
            #* regions at the wall
            d_windices = dgrid[tuple(windices.T)] # distance values at the window indices
            r_windices =region_labels[tuple(windices.T)] # region labels at the window indices
            # window_centers = [(window_points[np.argmax(d_windices[np.logical_and(r_windices[ag.labels_ == uwi]==regs[0],ag.labels_==uwi)])] + window_points[np.argmax(d_windices[np.logical_and(r_windices[ag.labels_==uwi]==regs[1],ag.labels_==uwi)])])/2.0 for uwi in np.unique(ag.labels_) if uwi != -1]

            # print('These are window centers', window_centers)
            window_centers = []
            for uwi in np.unique(ag.labels_):
                # for each cluster that is nota noise, do this 
                if uwi != -1:
                    # np.argmax(d_windices[ag.labels_==uwi])



                    uwi_dvals = d_windices[ag.labels_==uwi] # distance values in the cluster with label uwi
                    uwi_points = window_points[ag.labels_==uwi] # points in the cluster with label uwi
                    uwi_r = r_windices[ag.labels_==uwi]
                    
                    # for points in first region

                    # uwi_r1_indices =np.transpose((uwi_r==regs[0])).nonzero()
                    uwi_r1_indices =(np.where)(uwi_r==regs[0])[0]
                    # print('uwi_r1_indices are ', uwi_r1_indices)
                    # print(uwi_r1_indices.shape)
                    uwi_r1_dvals = uwi_dvals[uwi_r1_indices]
                    uwi_r1_points = uwi_points[uwi_r1_indices]
                    uwi_r1_max_indices = np.argmax(uwi_r1_dvals) #np.where(uwi_r1_dvals == np.max(uwi_r1_dvals))[0] # indices of the maximum distance values in the cluster with label uwi
                    # print(uwi_r1_max_indices)
                    uwi_r1_max_points = uwi_r1_points[uwi_r1_max_indices]

                    
                    
                    # for points in second region
                    uwi_r2_indices =(np.where)(uwi_r==regs[1])[0]
                    uwi_r2_dvals = uwi_dvals[uwi_r2_indices]
                    uwi_r2_points = uwi_points[uwi_r2_indices]
                    uwi_r2_max_indices = np.where(uwi_r2_dvals == np.max(uwi_r2_dvals))[0] # indices of the maximum distance values in the cluster with label uwi
                    uwi_r2_max_points = uwi_r2_points[uwi_r2_max_indices]
                    
                    # take the centroid of the two maxima of each region close to the window
                    uwi_center = (uwi_r1_max_points + uwi_r2_max_points)/2

                    # uwi_center = np.mean(uwi_max_points, axis=0) # mean of the points where the distance is at maximum if ther are multiple a.k.a the window is symmetric on each side
                    # add that to the list of windows
                    window_centers.append(uwi_center) 

            # window_centers = np.vstack(window_centers)
        # else:
            # window_centers = window_points
            window_centers = np.vstack(window_centers)
            window_fractional = get_fractional_coordinates(window_centers,ase_atoms)
            window_fractional = np.round(window_fractional, decimals=4)
            window_radii = np.array(dinterp(window_fractional))
            return np.array([regs, True,'internal', window_centers, window_radii], dtype='object') # returns a list of [[r1,r2], connected or not, array of window_centers]
        else:
            return np.array([regs, False,'internal', None, None], dtype='object')

    from skimage.segmentation import find_boundaries
    import numpy as np
    import itertools
    import ase
    # from sklearn.cluster import DBSCAN 
    # from sklearn.preprocessing import StandardScaler
    from icecream import ic
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.decomposition import PCA
    from mofography import interpolate_labels, get_fractional_coordinates
    import dask.bag as db
    A_unit = ase.geometry.complete_cell(ase_atoms.get_cell()).T
    from icecream import ic
    minimum_window_separation = distance_threshold
        
    dinterp = interpolate_dgrid(dgrid)
    # * Find the outer boundary of each individual region
    regions_to_check = np.delete(np.unique(region_labels),0)
    ic(regions_to_check)
    # outer = {reg:find_boundaries(region_labels == reg, mode='outer') for reg in regions_to_check}
    outer = {reg:find_boundaries(region_labels == reg, mode='outer') for reg in regions_to_check}
    ic(len(outer))
    # * Check for connectedness between regions
    check_list = np.vstack(list(itertools.combinations(regions_to_check, 2))) # list of pairs of region labels
    connections_list = [check_connected(c) for c in check_list] # now this step could cause an overflow
    # connections_list = db.from_sequence(check_list).map(check_connected).compute()
    
    # connected_regions = [[c[0], [c[2]] for c in connections_list if c[1] == True] # Make a new list with all the 
    
    #* Windows on the wall are to be calculated
    if wall_windows:
        wwl= get_wall_windows2(regions=region_labels, maxima=maxima, dgrid=dgrid, ase_atoms=ase_atoms, wall_thickness=wall_thickness)
        if len(wwl)>0:
            connections_list.append(wwl)


    # #* we have one final thing to do 
    # # * if both the regions forming an internal window are on the same wall, then we need to set the window center and the radius to 
    # # * that on the wall, since centroid


    # for c in connections_list:
    #     if np.logical_and(c[1], c[2]=='internal'):

    return np.vstack(connections_list)


def interpolate_dgrid(dgrid):
    """Create a scipy.interpolator.RegularGridInterpolator, using 'nearest' method for the region labels (or any 3D numpy grid) 
    over a unit box ([0,1]) in all axes.

    :param regions: 3D numpy array of region labels
    :type regions: numpy.array 
    :return: a scipy RegularGridInterpolator object that works on the fractional coordinates of atoms
    :rtype: scipy.interpolate.RegularGridInterpolator
    """    
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator as RGI
    xx = np.linspace(0,  1,dgrid.shape[0])
    yy = np.linspace(0,  1,dgrid.shape[1])
    zz = np.linspace(0,  1,dgrid.shape[2])
    dl = RGI((xx, yy, zz),dgrid, method='linear')

    return dl

def get_wall_windows2(regions, maxima, dgrid, ase_atoms, wall_thickness=2):
    import numpy as np
    from icecream import ic
    # from .helper import find_coord_from_indices, interpolate_labels
    wall_windows = []
    shape =regions.shape
    dinterp =interpolate_dgrid(dgrid)
    for a in [0,1,2]:
        ic.disable()
        ic(a)
        llr = regions.take(indices=0 ,axis=a) #* left wall
        llr_d = dgrid.take(indices=0, axis=a) #* left wall distances 
        rrr = regions.take(indices=shape[a]-1 ,axis=a)  #* right wall
        rrr_d = dgrid.take(indices=shape[a]-1, axis=a) #* right wall distances 

        map_list = np.unique(list(zip(llr.flatten(), rrr.flatten())), axis=0)
        not_background = [m for m in map_list if np.all(m!=0)]
        # * there should be some non background mappings in the direction
        if len(not_background) > 0:
            map_list = np.vstack([m for m in map_list if np.all(m!=0)])

            # * flags to check whether it is window or a just a pocket split across walls
            flags = np.vstack([[maxima[k - 1][a]>wall_thickness, maxima[v-1][a]<(shape[a]-wall_thickness)] for k,v in map_list])    
            
            # * this line makes sure maxima are both well inside
            flags = np.logical_and(flags[:,0],flags[:,1])
            
            # window_indices = [np.vstack(np.where(llr==m[0])).mean(axis=1) for i,m in enumerate(map_list) if flags[i]]
            wregs_l  = [(np.array([m[0], m[1]]), True,'periodic') for i, m in enumerate(map_list) if flags[i]]
            wregs_r  = [(np.array([m[1], m[0]]), True,'periodic') for i, m in enumerate(map_list) if flags[i]]
            # instead of taking the mean of the points at the left wall, let's find the location of the maximum in the distance
            # values for the points on that plan
            # window_indices = [np.vstack(np.where(llr==m[0])).mean(axis=1).T for i, m in enumerate(map_list) if flags[i]]
            window_indices_list = [np.argwhere(llr==m[0]) for i, m in enumerate(map_list) if flags[i]]
            window_indices = [indices[np.argmax(llr_d[indices[:,0], indices[:,1]])] for indices in window_indices_list]

            # print('Detected wall windows are :{}'.format(window_indices))

            ic(window_indices)
            if np.any(flags):
                window_indices_l = np.insert(np.vstack(window_indices), [a], np.zeros((len(window_indices),1)), axis=1  )
                window_indices_r = np.insert(np.vstack(window_indices), [a], (shape[a] -1)*np.ones((len(window_indices),1)), axis=1 )
                window_radii_l = np.vstack(dinterp(window_indices_l/regions.shape))
                window_radii_r = np.vstack(dinterp(window_indices_r/regions.shape))
                window_coords_l = find_coord_from_indices(window_indices_l, shape=shape, ase_atoms=ase_atoms)
                window_coords_r = find_coord_from_indices(window_indices_r, shape=shape, ase_atoms=ase_atoms)
                ic(window_coords_l)
                ic(window_coords_r)
                window_rows_l = np.vstack([[*w, np.array([ window_coords_l[i] ]), window_radii_l[i]] for i,w in enumerate(wregs_l)])
                window_rows_r = np.vstack([[*w, np.array([ window_coords_r[i] ]), window_radii_r[i]] for i,w in enumerate(wregs_r)])
                wall_windows.append(window_rows_l) 
                wall_windows.append(window_rows_r) 
                
    if len(wall_windows)>0:
        wall_windows =  np.vstack(wall_windows)
    

    return wall_windows
def rag_from_connections(connections, maxima, only_use_internal=True):
    
    """ Create a basic networkx graph (with no properties annotated) from the array of connections
    :type connections: np.array of shape (Nregions-Choose-2 ,1) where Nregions is the number of unique regions. Each row has 
    three elements, the region indices [r1, r2], a boolean flag showing whether they are connected, and window coordinate(s) if the regions are connected, and None is they are not.
    :param connections: numpy array

    :raises:

    :rtype: a graph object of the class networkx.Graph
    """
    import networkx as nx 
    import numpy as np
    G = nx.Graph()
    # number_of_regions = np.max([c[0] for c  in connections])
    number_of_nodes = len(maxima) #np.max([c[0] for c  in connections])
    # print(number_of_nodes)
    # nx.set_node_attributes(G, values={i+1:i+1 for i in range(number_of_nodes)}, name='labels')
    G.add_nodes_from(range(1, number_of_nodes+1)) 
    if only_use_internal:# range stops one short
        G.add_edges_from([c[0] for c in connections if np.logical_and(c[1], c[2] =='internal')])
    else:
        G.add_edges_from([c[0] for c in connections if c[1]])
    # add labels property to rag
    nx.set_node_attributes(G, values={i+1:[i+1] for i in range(len(G.nodes))}, name='labels')
    return G
# %%
connections = mgr.connections_from_regions_and_dgrid(regions, dgrid, maxima, mat_atoms)
# %%
connections
# %%
rag  = mgr.rag_from_connections(connections, maxima)
# %%

rag.nodes(data=True), rag.edges(data=True)
# %%

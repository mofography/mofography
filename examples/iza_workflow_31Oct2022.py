# make the mofography imports 
#%%
import sys
import os

import numpy as np
# * Put your path here
path_to_mgr = os.path.abspath('c:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')

sys.path.append(path_to_mgr)

import mofography as mgr 
#%%
# read the iza database and the list of ase atoms


import pandas as pd
import pickle 


df = pd.read_feather('C:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography/examples/iza_df.feather')
iza_atoms = pickle.load(open('C:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography/examples/iza_atoms_list.pkl','rb'))

# %%
# let's look at the workflow again for just one object
import numpy as np
from ase.data import covalent_radii
from skimage.measure import regionprops 
from scipy.ndimage import zoom
import dask.array as da 

spacing = 0.5
block_size = 20000
mask_thickness =1.5
bins = np.linspace(0,30, 50)

mat_atoms = iza_atoms[54]
# mat_atoms = np.random.choice(iza_atoms)
radii = covalent_radii[mat_atoms.get_atomic_numbers()]
grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing =spacing, block_size=block_size)
grid = da.from_array(zoom(grid.compute(), 2), chunks=(20,20,20))

regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=mask_thickness, h=0.2)#, mask_thickness=1.5,exclude_border=False, min_distance=10 , threshold_abs=3, p_norm=2)
rag2 = mgr.get_rag_from_regions_and_grid(regions, grid.compute(),mat_atoms,maxima)   

#%%
pbc_groups = mgr.get_periodic_map(regions=regions, maxima=maxima, maxima_radii=grid.compute()[tuple(maxima.T)], return_conn=True, wall_thickness=4)
pmap= mgr.get_periodic_map(regions=regions, maxima=maxima, maxima_radii=grid.compute()[tuple(maxima.T)], return_conn=False, wall_thickness=4)

regions_pbc = np.vectorize(pmap.get)(regions)
props = regionprops(regions)
rag2 = mgr.add_psd_to_rag_with_pbc(rag = rag2, props = props, pbc_groups= pbc_groups, bins = bins, ase_atoms = mat_atoms, shape =regions.shape)
rag2 = mgr.add_clustering_to_rag(rag2)


# %%

# let's  plot and visualize everything 

import plotly.graph_objects as go
import plotly.io as pio 
pio.renderers.default = "browser"

mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, show_rings=False, stylize=True)
# rag_trace=mgr.plot_rag(rag2,plot_with='maxima', name='p-graph') 
surface_trace=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions, stylize=True)
# surface_trace=mgr.plot_surface_dgrid(grid.compute(), 1.75, select='all',ase_atoms=mat_updated,region_labels=regions_pbc, stylize=True)
# surface_trace=mgr.plot_surface_dgrid(grid.compute(), 1.75, select=list(pbc_groups[0]),ase_atoms=mat_updated,region_labels=regions, stylize=True)
cell_trace=mgr.plot_cell(mat_updated)
maxima_trace = mgr.plot_maxima(maxima = maxima, dgrid = grid.compute(),ase_atoms =mat_updated)

# fig = go.Figure(data=[*frame_trace,rag_trace,surface_trace,cell_trace,maxima_trace]).update(layout=mgr.white_layout).show()
fig = go.Figure(data=[*frame_trace,surface_trace,cell_trace,maxima_trace]).update(layout=mgr.white_layout).show()
#%%
go.Figure(data=mgr.plot_psds(rag2, pbc_groups)).show()
# %%
# conver these theps into a prefect flow with no subtasks, let's keep things simple for now
# from prefect import task, flow
import numpy as np
import dask.array as da
from ase.data import covalent_radii
from scipy.ndimage import zoom
from skimage.measure import regionprops
def get_rag_from_atoms(mat_atoms, spacing = 0.5, block_size=20000, mask_thickness=1.5 , bins = np.linspace(0,30, 50)):
    try:
        radii = covalent_radii[mat_atoms.get_atomic_numbers()]
        grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing =spacing, block_size=block_size)
        grid = da.from_array(zoom(grid.compute(), 2), chunks=(20,20,20))
        regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=mask_thickness, h=0.25)
        rag2 = mgr.get_rag_from_regions_and_grid(regions, grid.compute(),mat_atoms,maxima)
        rag2.graph['name'] = mat_atoms.info['name']
        pbc_groups = mgr.get_periodic_map(regions=regions, maxima=maxima, maxima_radii=grid.compute()[tuple(maxima.T)], return_conn=True, wall_thickness=4)
        props = regionprops(regions)
        rag2 = mgr.add_psd_to_rag_with_pbc(rag = rag2, props = props, pbc_groups= pbc_groups, bins = bins, ase_atoms = mat_atoms, shape =regions.shape)
        rag2 = mgr.add_clustering_to_rag(rag2)
        return rag2
    except: 
        return 'NaN'
#%%
import numpy as np
import dask.array as da
from ase.data import covalent_radii
from scipy.ndimage import zoom
from skimage.measure import regionprops
def get_rag_from_atoms_test(mat_atoms, spacing = 0.5, block_size=20000, mask_thickness=1.5 , bins = np.linspace(0,30, 50)):
    radii = covalent_radii[mat_atoms.get_atomic_numbers()]
    grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing =spacing, block_size=block_size)
    grid = da.from_array(zoom(grid.compute(), 2), chunks=(20,20,20))
    regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=mask_thickness, h=0.25)
    rag2 = mgr.get_rag_from_regions_and_grid(regions, grid.compute(),mat_atoms,maxima)
    rag2.graph['name'] = mat_atoms.info['name']
    pbc_groups = mgr.get_periodic_map(regions=regions, maxima=maxima, maxima_radii=grid.compute()[tuple(maxima.T)], return_conn=True, wall_thickness=4)
    props = regionprops(regions)
    rag2 = mgr.add_psd_to_rag_with_pbc(rag = rag2, props = props, pbc_groups= pbc_groups, bins = bins, ase_atoms = mat_atoms, shape =regions.shape)
    rag2 = mgr.add_clustering_to_rag(rag2)
    return rag2

# %%
from dask import delayed 

from dask.distributed import progress
from distributed import Client, LocalCluster

# client = Client()

#%%
graphs = []
for a in iza_atoms:
    graphs.append(delayed(get_rag_from_atoms)(a))
#%%

from dask.diagnostics import ProgressBar
import dask
with ProgressBar():
    iza_graphs =dask.compute(graphs)
#%%
futures = client.persist(graphs)
progress(futures)


# %%
# grab the data
iza_graphs = [f.compute() for f in futures]


# %%
# check how many failed 

np.where([i=='NaN' for i in iza_graphs])
# %%
# save graphs to file 

import pickle

with open('iza_graphs.pkl', 'wb') as f:
    pickle.dump(iza_graphs, f)
# %%
get_rag_from_atoms_test(iza_atoms[54])
# %%
cluster = LocalCluster()
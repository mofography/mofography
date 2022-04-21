#%% [markdown]
### Scalable Workflow using Prefect 
# This notebooks shows the creation of a scalable workflow for computing 
# graphs for a subset of MOFDB
#%% [markdown]
#### Make the necessary imports, inclucing mofography
#%%
import sys
import os
import numpy as np 

import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default ='browser'
path_to_mgr = os.path.abspath('c:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
sys.path.append(path_to_mgr)
#%% [markdown]
#### Compile a dataset easily from MOFDB
# Here we are pulling the IZA database asynchronously
#%%
import mofography as mgr 
params = dict(vf_min=0.2, vf_max=1, database='IZA')
df = mgr.get_data_from_mofdb(params=params)
#%% [markdown]
#### EDA (Optional)
# Time to analyze the data interactively to see correlations, make plots etc.
#%%

# import dtale 
# dtale.show(df).open_browser()
#%% [markdown]
# Convert the cifs to ASE atoms objects
# - ASE atom objects are very useful for all of our tasks
# - If some cifs are not readable due to missing or obscure information, flag and remove those

#%%
atoms_list = mgr.atoms_from_mofdb_df(df)
#%%
cannot_read= np.where(atoms_list=='NaN')[0]
df = df.drop(cannot_read, axis=0)
atoms_list = np.delete(atoms_list, cannot_read, axis=0)

#%%
#### Filter some rows for testing
# Pick only materials with large pores, purely for testing the flow
#%%
atoms_list = atoms_list[df['lcd']>10]
df = df[df['lcd']>10]

#%% [markdown]
#### Now we should think about creating a workflow. 
# These are the steps in the workflow
# - Compute distance grids
# - Compute the region segmentationa
# - Compute the pore graph
# - Annotate the pore graph with properties
# - Compute the graph embedding to get a mixed node properties

#%% [markdown]

#These are some of the  the parameters of the workflow (for each step)
# - Grid spacing, block size
# - mask thickness, absolute peak height
# - bins=100 (for the region histograms)


#%%
from prefect import Flow, task, unmapped

def get_node_coords(seg,shape, ase_atoms):
    """
    > The function takes in a segment (seg) and a shape (shape) and returns the coordinates of the nodes
    in the segment (node_coords)
    
    :param seg: a tuple of the form (segment_id, [list of node indices])
    :param shape: the shape of the unit cell
    :param ase_atoms: the ase atoms object of the MOF
    :return: The coordinates of the atoms in the segment.
    """
    import mofography as mgr
    node_coords = mgr.find_coord_from_indices(seg[1], shape=shape, ase_atoms=ase_atoms)
    return [node_coords.tolist() ]

def get_psds(seg,dgrid):
    """
    > For each region in the segmentation, get the histogram of the distance grid values within that
    region
    
    :param seg: a tuple of (regions,labels)
    :param dgrid: the distance grid
    :return: A list of tuples, each tuple containing the histogram data and the bin edges.
    """
    import numpy as np 
    psds = []    
    regions=seg[0]
    for r in range(1,np.unique(regions).max()): # * exclude region 0, which is the background (framework mask)
        psds.append(np.histogram(dgrid[np.where(regions==r)], bins=100, density=True))
    return list(psds)  



g = task(mgr.dgrid_from_atoms, name='grids', checkpoint=False)
c = task(lambda x:x.compute(), name='compute', checkpoint=True)
r = task(mgr.regions_from_dgrid, name='regions', nout=2)
r1= task(lambda x:x[0], name='get_first')
n = task(mgr.connections_from_regions, name='connect')
gr = task(mgr.graph_from_connections, name='graphs')
a = task(mgr.add_props_to_graph, name='add_props')
nc =task(get_node_coords, name='node_coords')
ew = task(lambda conn:[[c[2] for c in conn if c[1]]], name='edge_windows')
sh = task(lambda x:x.shape, name='get_shape')
ps = task(get_psds, name='get_psds')
ls = task(lambda x:list(*x), name='get_lists')


with Flow('test_grid') as flow1:
    grid = g.map(ase_atoms = atoms_list, spacing = unmapped(0.25), block_size=unmapped(10000))    
    grid1 = c.map(grid)
    shapes = sh.map(grid)
    seg = r.map(dgrid=grid1,mask_thickness=unmapped(1), threshold_abs=unmapped(0.5), min_distance=unmapped(4))
    regions1 = r1.map(seg)
    connections = n.map(region_labels = regions1, ase_atoms=atoms_list)
    graphs = gr.map(connections=connections )
    node_coords = nc.map(seg=seg, shape=shapes,ase_atoms=atoms_list)
    edge_windows = ew.map(connections)
    psds = ps.map(seg=seg, dgrid=grid1)
    graphs_updated = a.map(graphs, node_props=node_coords, node_prop_names=unmapped(['dhists']))#, edge_props=edge_windows, edge_prop_names=unmapped(['window_coords']))
    # grids = g.map(ase_atoms=atoms_list[0:10], spacing = unmapped(0.5), block_size=unmapped(10000))

#%% [markdown]
#### Visualize the flow before run
# see the calculation for each atoms object
#%%
flow1.visualize()


#%% [markdown]
# Run the flow

#%%
state = flow1.run()

#%% [markdown]
#### Visualize the state of the flow
# We can see which calculations succeded vs failed,
# And which step failed
# Note that individual failures don't affect other calculations

#%%
flow1.visualize(flow_state=state)

#%% [markdown]
#### Scaling up to an HPC using Dask
# We can swap out the default executor with Dask
# So scaling up to an HPC simply becomes, swapping out the LocalCluster
# with the SGDCluster
#%%
from prefect.engine.executors import DaskExecutor 
from dask.distributed import LocalCluster, Client

with Client(LocalCluster()) as client:
    executor = DaskExecutor(address=client.scheduler.address)
    state = flow1.run(executor=executor)
#%%

flow1.visualize(flow_state=state)
#%%
#%% [markdown]
#### Keeping track of the workflows using Prefect Server
# We cna use the prefect server and UI to see the workflows



#region #* GARBAGE
#%%
def get_data_from_mofdb(url='https://mof.tech.northwestern.edu/mofs/#.json', params={'vf_min': 0.1, 'vf_max': 1.0, 'database':'IZA','page':1} , headers={'loading': 'cm3(STP)/cm3', 'pressure': 'bar'},npages=None, mgr_path='c:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography'):
    """
    It takes in a url, parameters, headers, number of pages and a path to the mofography repository and
    returns a pandas dataframe of the data. 
    
    The function is asynchronous and uses the aiohttp library to make multiple requests to the MOFDB
    API. 
    
    The function is also dependent on the mofography repository. 
    
    
    :param url: The url of the MOF database, defaults to https://mof.tech.northwestern.edu/mofs/#.json
    (optional)
    :param params: a dictionary of parameters to be passed to the url. The default parameters are: {'vf_min': 0.1, 'vf_max': 1.0, 'database':'IZA','page':1}
    :param headers: This is the headers that you want to pass to the MOFDB API. The default is to get
    the loading and pressure data
    :param npages: number of pages to be scraped. If None, then the number of pages is determined by the
    function to return all the pages
    :param mgr_path: The path to the mofography repository. This is needed to import the mofography
    module, defaults to c:/Users/Arun Gopalan/OneDrive - The University of
    Manchester/my_cloned_repos/mofography (optional)
    :return: A dataframe with the results of the query.
    """

    from aiohttp import ClientSession
    import asyncio
    import nest_asyncio
    nest_asyncio.apply()
    import pandas as pd 
    import sys 
    if not (mgr_path in sys.path):
        sys.path.append(os.path.abspath(mgr_path))
    
    if npages==None:
        npages=mgr.get_number_of_pages(url, headers, params)
    
    # This shouldn't be called in a function (similar to initiating a Dask cluster)
    loop = asyncio.get_event_loop()
    # with ClientSession(headers= headers) as session:
    session = ClientSession(headers=headers)
    fut = asyncio.ensure_future(mgr.get_mofdb_data(session, url=url, params=params, headers=headers, npages=npages))
    loop.run_until_complete(fut)
    session.close()
    data_df = pd.DataFrame(list(np.hstack(([f['results'] for f in fut.result()]))))

    return data_df



#%%

def atoms_from_mofdb_df(mofdb_df):

    """
    This function takes a MOFDB dataframe from get_data_from_mofdb routine, which has a column called 'cif' that contains CIF data as strings, and
    returns a list of ASE atoms objects
    
    :param mofdb_df: the dataframe containing the MOFDB data
    :return: A list of atoms objects.
    """

    def read_cif_string(cif):
        """
        It takes a string of CIF data and returns an ASE atoms object.
        
        :param cif: The CIF file as a string
        :return: the ASE atoms object.
        """
        from ase.io import read 
        import io
        try:
            output = read(io.BytesIO(bytes(cif, encoding='utf-8')), format='cif')
        except: 
            output = 'NaN'
        return output

#endregion
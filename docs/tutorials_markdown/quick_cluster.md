# Quick Clustering Workflow with Periodic Segmentation

This notebook demonstrates the recommended explicit workflow for pore detection and clustering in nanoporous materials using MOFography's new **periodic segmentation** approach.

## What's New?

- Uses `do_periodic_segmentation_rust()` - **5-20x faster** with Rust backend
- Properly handles periodic boundary conditions from the start
- More accurate pore segmentation
- Explicit step-by-step workflow (not a black box)

**Note:** This replaces the old `quick_cluster()` "black box" function. For the legacy approach, see `applying_pbc.ipynb`.

## Workflow Overview

The steps used in this workflow are:

1. Load the structure
2. Compute the distance grid
3. Run periodic segmentation (detects maxima, segments regions, handles PBC, computes connections)
4. Extract results (maxima, regions, connections)
5. Create the pore graph (RAG) from connections
6. Add properties and create visualizations

## Tips and Tricks #1

All the possible arguments to the periodic segmentation function can be looked at using the `help` function:

```python
help(mgr.do_periodic_segmentation_rust)
```


## Imports


```python
# add path to mofography
import os, sys
import numpy as np

# Use relative path from docs folder
path_to_mgr = os.path.abspath('..')
import mofography as mgr

from ase.io import read

path_to_cif = os.path.join(path_to_mgr, 'examples', 'TestCIFs', 'P1-IRMOF-1.cif')

```

## Read the CIF file


```python
mat_atoms = read(path_to_cif)

# figures = mgr.quick_cluster(mat_atoms, return_these_vars=[figures])
```

    c:\Users\arung\miniforge3\envs\magic\Lib\site-packages\ase\io\cif.py:408: UserWarning: crystal system 'triclinic' is not interpreted for space group Spacegroup(1, setting=1). This may result in wrong setting!
      warnings.warn(


## Set Parameters

First, let's set up the parameters for the periodic segmentation workflow:


```python
# Grid and segmentation parameters
spacing = 0.5  # Grid spacing in Angstroms

# Segmentation parameters
mask_thickness = 1.0  # Void space threshold (Angstroms)
h = 0.2  # h-maxima persistence threshold
threshold_abs = 1.5  # Absolute threshold for maxima detection

# Periodic segmentation parameters
minimum_window_separation = 1.5  # Minimum distance between windows
wall_thickness = 0.08  # Boundary detection thickness (fraction of cell)
extra_padding = 10  # Extra padding in Angstroms

```

## Compute Distance Grid


```python
# Compute distance grid using CPU backend
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=50000)

# Apply zoom and smoothing
dgrid = mgr.zoom_grid(dgrid, zoom_factor=1).compute()
dgrid = mgr.periodic_gaussian_filter(dgrid, sigma=spacing*2)

print(f"Distance grid shape: {dgrid.shape}")
print(f"Distance grid range: {dgrid.min():.2f} to {dgrid.max():.2f} Angstroms")

```

## Run Periodic Segmentation

Now let's run the periodic segmentation workflow. This single function call:
- Detects local maxima (with PBC handling)
- Performs watershed segmentation
- Identifies connections between regions
- Returns all outputs in a structured dictionary


```python
results = mgr.do_periodic_segmentation_rust(
    dgrid,
    mat_atoms,
    spacing=spacing,
    mask_thickness=mask_thickness,
    h=h,
    threshold_abs=threshold_abs,
    minimum_window_separation=minimum_window_separation,
    wall_thickness=wall_thickness,
    extra_padding=extra_padding,
    return_padded_regions=True,
    compute_connections=True
)

print("Periodic segmentation complete!")
print(f"Found {len(results['maxima_pbc'])} unique periodic maxima")
print(f"Found {len(results['unique_regions'])} unique regions")
print(f"Found {len(results['pbc_groups'])} PBC groups")
print(f"Found {len(results['connections'])} region connections")

```

## Extract Results

The periodic segmentation returns a dictionary with all the outputs. Let's extract the key results:


```python
# Extract all results from the dictionary
maxima_pbc = results['maxima_pbc']  # Representative maxima for each PBC group
regions_pbc = results['regions_pbc']  # Segmented regions with PBC applied
regions = results['regions']  # Segmented regions without PBC
connections = results['connections']  # List of region connections
all_maxima = results['all_maxima']  # All periodic maxima images
unique_regions = results['unique_regions']  # Regions influencing the unit cell
maxima_info = results['maxima_info']  # Detailed maxima information
pbc_groups = results['pbc_groups']  # List of PBC groups
close_maxima = results['close_or_inside_maxima']  # Maxima near/in unit cell

# Maxima influencing the unit cell (for visualization and RAG)
maxima = all_maxima[:regions.max()]

print("Results extracted:")
print(f"  - maxima_pbc: {maxima_pbc.shape}")
print(f"  - regions_pbc: {regions_pbc.shape}, max={regions_pbc.max()}")
print(f"  - regions: {regions.shape}, max={regions.max()}")
print(f"  - connections: {len(connections)} pairs")
print(f"  - all_maxima: {all_maxima.shape}")
print(f"  - maxima_info: {maxima_info.shape}")

```


## Create Visualizations

Now let's create visualizations of the results. This figure shows:
- The types of pockets, where all the pockets detected to be of the same type and shown in the same color. 

There is more information available in the figure if you hover over any of the surfaces or the framework atoms.


```python
# For visualization, let's create the figure from scratch
import plotly.graph_objects as go
import plotly.express as px

# Set up colorscale
cscale = px.colors.qualitative.Pastel

# Update structure with covalent bonds
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

# Create traces
cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)

# Plot maxima (those influencing the unit cell)
maxima_trace = mgr.plot_maxima(
    dgrid, maxima, mat_atoms,
    name='Local maxima',
    show_spheres=False,
    marker_colorscale=cscale
)

# Plot surfaces - regions with PBC applied
surface_trace_pbc = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=1.5,
    region_labels=regions_pbc,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name='Regions with PBC',
    stylize=True
)

# Combine into figure
fig = go.Figure(
    data=[
        cell_trace,
        surface_trace_pbc,
        maxima_trace,
        *frame_trace
    ],
    layout=mgr.dark_layout
)

fig.update_layout(scene_aspectmode='data')
fig.show(renderer='png')
# fig.show()  # Uncomment for interactive view
```


    
![png](/Users/arungopalan/repos/mofography/docs/tutorials_markdown/quick_cluster_10_0.png)
    


Any of the plot traces can be recovered using the data variable of the Plotly figure object (`fig.data`). For example, to see the names of the traces:


```python
[d['name'] for d in fig.data]
```




    ['Raw segmentation',
     'PBC groups',
     'Cavity types',
     'Framework',
     'Aromatic rings',
     'Aromatic highlights',
     'Aromatic highlights',
     'Aromatic highlights',
     'Aromatic highlights',
     'Aromatic highlights',
     'Aromatic highlights',
     'Aromatic highlights',
     'Aromatic highlights',
     'Aromatic highlights',
     'Aromatic highlights',
     'Aromatic highlights',
     'Aromatic highlights',
     'cell',
     'local maxima']



To set a different color scheme, say for the `Cavity types` mesh trace to `Pastel`, you can do that also post generation.


```python
import plotly.express as px
# Modify the regions surface trace (index 1 in our figure)
fig.data[1]['colorscale'] = px.colors.qualitative.Pastel
fig.show(renderer='png')
```


    
![png](/Users/arungopalan/repos/mofography/docs/tutorials_markdown/quick_cluster_14_0.png)
    


## Creating the Region Adjacency Graph (RAG)

The periodic segmentation workflow already computed the connections between regions. We can use these to create the pore graph (Region Adjacency Graph):


```python
# Create RAG from the connections we already computed
rag = mgr.rag_from_connections(connections, all_maxima)

print(f"RAG has {rag.number_of_nodes()} nodes and {rag.number_of_edges()} edges")
```


We can add the pore graph to our visualization:


```python
# Add RAG trace to the existing figure
rag_trace = mgr.plot_rag(rag, name='Pore graph')
fig.add_trace(rag_trace)
fig.show(renderer='png')
```



## Where to Next?

Now that you've seen the explicit periodic segmentation workflow, you can:

1. **Explore the Basic Tutorials** to understand each step in detail:
   - `distance_grid.ipynb` - Learn about distance grid computation
   - `region_segmentation.ipynb` - Deep dive into segmentation methods
   - `region_adjacency_graph.ipynb` - Working with pore graphs
   - `connections_and_windows.ipynb` - Detecting windows between pores

2. **Check the Advanced Tutorials** for specific applications:
   - `pockets_detection.ipynb` - Identify different pocket types
   - `pld_and_percolation.ipynb` - Pore limiting diameter and percolation analysis
   - `hts_using_mofography.ipynb` - High-throughput screening

3. **Read the Reference Documentation**:
   - `periodic_segmentation_notes_ai.md` - Complete guide to periodic segmentation
   - `outputs_from_periodic_segmentation.md` - Output reference

**Note:** The old `quick_cluster()` "black box" function has been replaced by this explicit workflow. If you need the legacy approach, see `applying_pbc.ipynb`.


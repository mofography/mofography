mofography 🧪🔬
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/mofography/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/mofography/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/mofography/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/mofography/branch/master)

<!-- Add logo -->

![](./docs/logo.png)

A powerful Python package for the analysis, visualization, featurization, and machine learning on molecular pores in nanoporous materials. 🔍🧠

## Trying out mofography 🚀

*Note: In all the steps below, feel free to replace `conda` with `mamba`, as it is much faster at resolving dependencies and installing packages. To install mamba use `conda install -c conda-forge -n base mamba`*

You can simply try out mofography after installing just the dependencies. In fact, we recommend this approach, since mofography is still in active development.

### 1. Set up the environment 🛠️

Create a new conda environment with Python 3.9:

```bash
conda create -n mgr_env python=3.9
conda activate mgr_env
```

### 2. Install dependencies 📦

Using conda:

```bash
conda install -c conda-forge pymatgen ase dask dask-jobqueue networkx scikit-image tqdm aiohttp icecream sparse scikit-learn rich trimesh pyvista plotly furl requests freud cupy
```

Using pip:

```bash
pip install pymatgen ase dask dask-jobqueue networkx scikit-image tqdm aiohttp icecream sparse scikit-learn rich trimesh pyvista plotly furl requests freud-analysis cupy
```

*Note*: If you don't have a GPU, exclude `cupy` from the list of dependencies.

### 3. Clone the repo 📥

Download the zip file or use git:

```bash
git clone https://github.com/arung-northwestern/mofography.git
```

If you already have a cloned repo, update it:

```bash
cd mofography
git pull
```

### 4. Add mofography to your Python path 🐍

In your scripts and notebooks, add:

```python
import sys
import os
sys.path.append(os.path.abspath('path_to_mofography'))
```

Replace `path_to_mofography` with the relative path to the mofography directory on your local machine.

### 5. Test the package 🧪

Try a full clustering workflow with default parameters on the IRMOF-1 structure:

```python
import mofography as mgr
from ase.io import read
mat_atoms = read('path_to_mofography/examples/TestCIFs/P1-IRMOF-1.cif')
figures = mgr.quick_cluster(mat_atoms)
```

## Tutorials 📚

Find tutorials in the [Docs](https://arung-northwestern.github.io/mofography/docs) section. Here's a recommended order:

### Basic tutorials 🌱

1. [Quick cluster notebook](docs/quick_cluster.ipynb) - Get a feel for the whole thing
2. [Basic structure analysis notebook](docs/basic_structure_io.ipynb)
3. [Basic structure visualization notebook](docs/quick_structure_visualization.ipynb)
4. [Distance grid notebook](docs/distance_grid.ipynb)
5. [Region segmentation notebook](docs/region_segmentation.ipynb)
6. [PBC notebook](docs/applying_pbc.ipynb) - Applying Periodic Boundary Conditions
7. [RAG notebook](docs/region_adjacency_graph.ipynb) - Region Adjacency Graph (Pore graph)
8. [Connections and windows notebook](docs/connections_and_windows.ipynb)
9. [Pocket type detection notebook](docs/pockets_detection.ipynb)

### Advanced tutorials 🚀

These can be completed in any order after finishing the basic tutorials:

- [Advanced structure visualization notebook](docs/framework_visualization_example.ipynb)
- [Framework subset visualization notebook](docs/plot_subset_as_mesh_tutorial.ipynb)
- [PLD and percolation history notebook](docs/pld_and_percolation.ipynb)
- [Compute pore diameters notebook](docs/compute_pore_diameters.ipynb)
- [MOFDB notebook](docs/mofdb_example.ipynb)
- [HTS Dask notebook](docs/hts_using_mofography.ipynb)
- [Functional groups notebook](docs/detect_functional_groups_in_dataset.ipynb)
- [GCN tutorial](docs/graph_convolution_on_rags_ccr.ipynb)

## Copyright ©️

Copyright (c) 2022, Arun Gopalan

### Acknowledgements 🙏

Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.

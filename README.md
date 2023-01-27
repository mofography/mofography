mofography
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/mofography/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/mofography/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/mofography/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/mofography/branch/master)


A powerful python package for the analysis, visualization, featurization and machine learning on molecular pores in nanoporous materials.

### Trying out mofography

*Note: In all the steps below, feel free to replace `conda` with `mamba`, as it is much faster at resolving dependencies and installing packages. To install mamba use `conda install -c conda-forge -n base mamba`*

You can simply try out mofography without actually installing anything. In fact, I would recommend this approach, since mofography is still in active development
- Install the dependencies using conda or pip. You may create a new conda environment with python 3.9 for this, but it is not necessary. 

```bash
conde create -n mgr_env python=3.9
conda activate mgr_env
```

```bash
conda install -c conda-forge pymatgen ase dask dask-jobqueue networkx scikit-image tqdm aiohttp icecream cupy sparse scikit-learn rich numpy scipy copy trimesh pyvista plotly itertools pandas asyncio furl requests
```

- Clone the repo 
  Click on the green `Code` button and download the zip file. Unzip the file to get the repo folder, the path to which we can add and use. Alternatively, you can clone the repo using git. If you are using git, you can clone the repo using the following command

```bash
git clone https://github.com/arung-northwestern/mofography.git .
```

- If you already have a cloned repo, pull the latest changes, since the package is being updated frequently
  
```bash
cd mofography
git pull
```

- Add the path to mofography the your python path in your scripts and notebooks before you import,
  
```python
import sys
import os
sys.path.append(os.path.abspath('path_to_mofography'))
```

of course, replace `path_to_mofography` with the relative path to the mofography directory on your local machine.
`

- Import mofography to see if it works. There are CIF files in the `mofography/examples/TestCIFs` folder, which you can use to test the package. Let's use an full-on clustering workflow with default parameters on the IRMOF-1 structure. The details of the workflow are explained in the [Quick cluster](docs/quick_cluster.ipynb) tutorial.

```python
import mofography as mgr
from ase.io import read
mat_atoms = read('path_to_mofography/examples/TestCIFs/P1-IRMOF-1.cif')
figures = mgr.quick_cluster(mat_atoms)
```

### Where to find some tutorials: 

You can find some tutorials in [Docs](https://arung-northwestern.github.io/mofography/docs)

Currently, there are only a few tutorials, here's a recommended order to go through them:

- Get a feel for the whole thing [Quick cluster notebook](docs/quick_cluster.ipynb) 
- Basic structure analysis [Basic structure analysis notebook](docs/basic_structure_io.ipynb)
- Basic structure visualization [Basic structure visualization notebook](docs/quick_structure_visualization.ipynb)
- Distance grid calculation and visualization [Distance grid notebook](docs/distance_grid.ipynb)
- Region Segmentation [Region segmentation notebook](docs/region_segmentation.ipynb)
- Applying Periodic Boundary Conditions [PBC notebook](docs/applying_pbc.ipynb)
- Region Adjacency Graph (Pore graph) [RAG notebook](docs/region_adjacency_graph.ipynb)
- Detecting the types of pockets [Pocket type detection notebook](docs/pockets_detection.ipynb)
- Finer control over structure visualization [Advanced structure visualization notebook](docs/framework_visualization_example.ipynb)
- Advanced: Detecting channels in IRMOF-1, Complete Workflow [IRMOF channels notebook](docs/irmof_example_23aug2020.ipynb)
  
### Copyright

Copyright (c) 2022, Arun Gopalan


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.

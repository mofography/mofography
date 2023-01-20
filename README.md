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

```bash
clone https://github.com/arung-northwestern/mofography.git .
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

- Import mofography to see if it works. There are CIF files in the `mofography/examples/TestCIFs` folder, which you can use to test the package.

```python
import mofography as mgr


```

### Copyright

Copyright (c) 2022, Arun Gopalan


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.

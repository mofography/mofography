#%%
# make some immports
from pathlib import Path
from ase.io import read
import numpy as np
# set some paths, to the pdb file and the corresponding unit cell. 
path_to_mgr = Path('C:/Users/arung/OneDrive - The University of Manchester/my_cloned_repos/22sep2024/mofography/mofography')
import mofography as mgr
path_to_pdb = path_to_mgr / 'examples/Movie_LTA96_1.1.1_298.150000_0.000000_component_K_0.pdb'
path_to_cif = path_to_mgr / 'examples/TestCIFs/LTA.cif'
# read cif
mat_atoms = read(path_to_cif)
# compute distance grid
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=0.5, precompute= True)
# compute regions and maxima
regions, maxima = mgr.regions_from_dgrid(dgrid, mask_thickness=1.0)
maxima_radii = dgrid[tuple(maxima.T)]
# apply pbc
regions_pbc = mgr.apply_pbc(regions, maxima, maxima_radii, mat_atoms )
print('There are {} unique pbc groups detected (excluding background)', format(len(np.unique(regions_pbc)-1)))
# apply region map to pdb (grou|ing into pbc groups here)
# same will work for individual regions
# same will work for unique pocket types
# just use the corresponding regions array for those.
mgr.apply_region_map_to_raspa_pdb_wrapped(path_to_pdb, regions_pbc, mat_atoms)

# %%

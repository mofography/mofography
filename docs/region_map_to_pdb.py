#%% [markdown] 
# # Read RASPA pdb file and count in each region
# Read tehe PDB file into dictionary of lists using the 'read_raspa_pdb' function. 
# Then convert into a list of ASE atoms objects using the 'get_atoms_list_for_pdb' function. 
# if the symbols in your PDB is not standard atomic symbols, ASE would struggle. to avoid this use the 'symbol_map' argument to point out which custom symbols to map to the standard atomic symbols. If same symbol could mean two things, like 'Hf' for front hydrogen and 'Hf' for hafniu, the one that is defined in the symbol map will get preference. Gives a lot of freedom to use these functions more easily. 
#%%
# make some immports
from pathlib import Path
from ase.io import read
import numpy as np
# set some paths, to the pdb file and the corresponding unit cell. 
# path_to_mgr = Path('C:/Users/arung/OneDrive - The University of Manchester/my_cloned_repos/22sep2024/mofography/mofography')
path_to_pdb = Path('C:/Users/arung/OneDrive - The University of Manchester/my_projects/kaihang/Movie_tobmof-10008_1.1.1_160.000000_500000.000000_allcomponents.pdb')
import mofography as mgr
data = mgr.read_raspa_pdb(path_to_pdb, symbol_map=dict(Hc='H'))
atoms_list = mgr.get_atoms_list_for_pdb(path_to_pdb, symbol_map=dict(Hc='H'))
# path_to_pdb = path_to_mgr / 'examples/Movie_LTA96_1.1.1_298.150000_0.000000_component_K_0.pdb'
# path_to_cif = path_to_mgr / 'examples/TestCIFs/LTA.cif'
#%% [markdown]
# Now that we have the atoms list, and the unit cell is same as the box in the PDB, we can just use the first atoms object as the reference cell for our shifting (which in this case is not needed obviously). Otherwise use a separate unit cell levl atoms object. 
#%%
# read cif
# mat_atoms = read(path_to_cif)
# reference atoms object for shifting to.
mat_atoms = atoms_list[0] # unit cell is same as box in the pdb, so just grab teh first atoms object. 
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
mgr.apply_region_map_to_raspa_pdb_wrapped(path_to_pdb, regions_pbc, mat_atoms, symbol_map=dict(Hc='H'))

#%% [markdown]
# The result shows the following
# - the counts in each configuration under each unique region (pbc_group or cluster id)
# - the mean and std of the counts across configurations for each unique region
# - the total mean
#%% [markdown]

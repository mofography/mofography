# in this notebook let's try to build an interpenetrated MOF on the CuBTC pocket ? 
# or let's try to se if we cna interpenetrate any of the IZA materialsfrom ipdb import set_trace
from pathlib import Path
import sys
import numpy as np
path_to_mgr = Path('/mnt/windows/Users/arung/shared_cloned_repos/repos/mofography_clones/periodic_segmentation/mofography')
sys.path.append(str(path_to_mgr))
import mofography as mgr
from ase.io import read
import plotly.express as px
import plotly.graph_objects as go
from ase.data import vdw_radii

# we can read the pdb in as a dictionary (useful, but not what we shall do today)
# use your pdb file here
path_to_pdb = Path('/home/arung014/Downloads/Movie_LTA_my_1.1.1_303.000000_0.000000_component_Na_1.pdb')
pdb_dict = mgr.read_raspa_pdb(path_to_pdb)
print(pdb_dict.keys())

# let's read in the pdb files as atoms instead 
pdb_atoms= mgr.get_atoms_list_for_pdb(path_to_pdb)

# let's make the trace for the guest molecules
index = 0 # index of the configuration from pdb you want to plot
atomic_number = 11 # not exactly right, should be the cation, but just for demo
spheres = mgr.spheres_for_atoms(pdb_atoms[index], scale=0.5)
spheres_trace = mgr.plot_mesh(spheres, stylize=False, name='Guest molecules')
# OPTIONAL: plot the bonds for the guest if polyatomic
num_atoms_in_guest = 2 # say you have two atoms in guest, 
# we need to make a bonds list (assuming they are in order)
def create_bonds_list(total, num_atoms_in_guest):
    return [list(range(i, i + num_atoms_in_guest)) for i in range(0, total, num_atoms_in_guest)]
bonds_list = create_bonds_list(len(pdb_atoms[index]), num_atoms_in_guest)
cylinders = mgr.cylinders_for_bonds(bonds_list,pdb_atoms[index])
bonds_trace = mgr.plot_mesh(cylinders, stylize=True, name='Guest bonds', showlegend=False)



# let's make the trace for the regions isosurface
path_to_framework = Path('/mnt/windows/Users/arung/shared_cloned_repos/papers/mofography_nature/jupyter_notebooks/cif_files/LTA96.cif')
mat_atoms = read(path_to_framework)
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing = 0.25)
extra_padding = 5
results = mgr.do_periodic_segmentation(dgrid, mat_atoms, spacing=0.25, return_padded_regions=True, mask_thickness=0.25, h=0.2, threshold_abs=1.0,  minimum_window_separation=1.5, wall_thickness=0.08, extra_padding=extra_padding)
regions_pbc = results['regions_pbc']

cscale  = px.colors.qualitative.Pastel # pick a different one if you don't like pastel
surface_value = 1.0 # needs to be greater than mask_thickness
surface_trace_pbc=mgr.plot_surface_dgrid(
    dgrid,
    isovalue=surface_value, 
    region_labels=regions_pbc,
    ase_atoms=mat_atoms,
    cscale=cscale,
    opacity = 0.65, 
    stepsize=1, 
    name='Regions PBC at ' + str(surface_value) + '\u212b',
    stylize=True,
)
cell_trace = mgr.plot_cell(mat_atoms, color='turquoise')

#plt framework atoms
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)
frame_trace [0].name = 'LTA (2x2x2)'
# or a subset of framework to show active sites or soomething   
# frame_subset = mgr.plot_subset_of_atoms(mat_atoms, subset=[1,2,4,5])

# let's make a figure out of these things
fig = go.Figure(data= [*frame_trace, cell_trace, spheres_trace, surface_trace_pbc], layout=mgr.white_layout)
# optional: if you have bonds in guest
# fig = go.Figure(data= [*frame_trace, cell_trace, spheres_trace, bonds_trace, surface_trace_pbc], layout=mgr.white_layout)
# optional: make it orthographic if needed
# fig = mgr.make_fig_orthographic(fig)
fig.show(renderer='firefox')

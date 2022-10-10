#%% [markdown]
# This code takes a cif file, and a pdb file of snapshots from raspa
# and generates cool looking traces, of selecetd adsorbates, near selected nodes, or
# other parts of hte framework.


#%%

import sys, os

#* add mofography to path (put your path here)
path_to_mgr = os.path.abspath('c:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
sys.path.append(path_to_mgr)
import mofography as mgr 
# %%

# Read the framework as a cif file
#Use pdb with index for reading pdb file snapshpts
from ase.io import read 

path_to_cif  = './CAXVUU/CO2/UFF/Movies/System_0/Framework_0_final_4_2_2_P1.cif'

mat_atoms = read(path_to_cif)#
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
mat_graph = mgr.sgraph_from_atoms(mat_updated)
# %% [markdown]
# Read teh adsorbate coordinates from raspa pdb 
# let's convert that snapshot into an ase atoms object, so
# that we can use the built in mofography functions without much hassle
#%%
pdb_file = './CAXVUU/CO2/UFF/Movies/System_0/Movie_CAXVUU_4.2.2_298.000000_10000.000000_component_CO2_0.pdb'
co2_data = mgr.read_raspa_pdb(pdb_file)
index = 8 # select a snapshot
natoms = 3 # number of atoms in the adsorbate
co2_coords = co2_data['coords'][index]
co2_sym = co2_data['symbols'][index]
co2_sym = [str(c) for c in co2_sym]
co2_cell = co2_data['cells'][index]
nmolecules = len(co2_coords)//natoms # numbre of co2 molecules
#conver to ase atoms object
from ase import Atoms 
co2_atoms = Atoms(symbols=co2_sym, positions=co2_coords, cell=co2_cell, pbc=[1,1,1])

# we also need teh co2 bonds for plotting
# in this case this is easy, since we know the bonds are 1-2 and 2-3
co2_bonds = []
for i in range(nmolecules):
    co2_bonds.append([i*natoms, i*natoms+1])
    co2_bonds.append([i*natoms+1, i*natoms+2])


# %% [markdown]
# - Now we have to pickoout the framework parts and adsorbate molecules 
# - the we want to plot. 
# - there is no secret formula for this, but we can do which ever criteria we want.
# one easy manual way is to visualize and then pick out the atoms we want.

#%%
# mesh trace fof the framework
frame_traces= mgr.plot_framework_as_mesh(mat_updated, show_rings= True, stylize =True)#, 
frame_traces[0]['name'] = path_to_cif.split('/')[-1].split('.')[0] #let add a name to the trace

# %%
def plot_snapshot_raspa_pdb(pdb_file, index, natoms, bonds_in_ref=None, trace_name='adsorbates', pbc=[1,1,1] ):
    """
    It reads a RASPA pdb file, and returns a plotly mesh trace
    
    :param pdb_file: the path to the pdb file in raspa output format
    :param index: the index of the snapshot you want to plot
    :param natoms: number of atoms in the molecule
    :param bonds_in_ref: a list of bonds in the reference structure. This is used to plot the bonds in
    the snapshot. If you don't have this, you can set it to None. and bonds will be put in the sequence, 0-1,1-2, ...
    :param trace_name: the name of the trace in the plotly figure, defaults to adsorbates (optional)
    :param pbc: the periodic boundary conditions, this is used to create an ase atoms object from the pdb snapshot.
    :return: A plotly figure object
    """
    import numpy as np 
    from .data import read_raspa_pdb
    config_data = mgr.read_raspa_pdb(pdb_file)
    coords = config_data['coords'][index]
    sym = config_data['symbols'][index]
    sym = [str(c) for c in sym]
    cell = config_data['cells'][index]
    nmolecules = len(coords)//natoms # numbre of co2 molecules
    #conver to ase atoms object
    from ase import Atoms 
    atoms = Atoms(symbols=sym, positions=coords, cell=cell, pbc=[1,1,1])

    # we also need teh co2 bonds for plotting
    # in this case this is easy, since we know the bonds are 1-2 and 2-3
    if bonds_in_ref!= None:
        bonds = []                                          
        for i in range(nmolecules):
            bonds.append([i*natoms, i*natoms+1])
            bonds.append([i*natoms+1, i*natoms+2])
    else:
        #add bonds based on the list of bonds in the reference molecule
        bonds = []
        for i in range(nmolecules):
            for b in bonds_in_ref:
                bonds.append([i*natoms+b[0], ])

    mol_mesh = mesh_for_atoms(ase_atoms =atoms, bonds = bonds)

    mol_mesh_trace = plot_mesh([mol_mesh] , name=trace_name, stylize=True)
    return mol_mesh_trace
                                            #         )
# %%
# co2_mesh = mgr.mesh_for_atoms(ase_atoms = co2_atoms, bonds = co2_bonds)

# co2_mesh_trace = mgr.plot_mesh([co2_mesh] , name='CO<sub>2</sub>', stylize=True)
co2_mesh_trace = plot_snapshot_raspa_pdb(pdb_file, index=8, natoms=3, bonds=None, trace_name='CO<sub>2</sub>')
#%%

#plot the unit cell

cell_trace = mgr.plot_cell(mat_updated, color='purple')
# %%
#add the above traces to plolty and make a figure
import plotly.graph_objects as go
import plotly.io as pio 
pio.renderers.default='browser'
# pio.renderers.default='notebook_connected'

# chose between white dark and epic layouts in mgr 
# don't forget the * in front ofthe framework trace to unpack frameowork and aromatic rings and hightlights.
fig = go.Figure(data=[*frame_traces, co2_mesh_trace]).update(layout=mgr.white_layout).update_layout(legend_groupclick='togglegroup')
fig.show()


# %%
#let add a nice backgroud to the figure for toc graphic
fig.update_layout(images=[mgr.make_bg_image(type='radial', cmap='Spectral_r')])
fig.show()

# %%
fig.write_html('caxvuu_uff_eqeq_co2.html')
# %%

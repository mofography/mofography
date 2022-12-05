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
#%%

# mesh trace fof the framework
frame_traces= mgr.plot_framework_as_mesh(mat_updated, show_rings= True, stylize =True)#, 
# frame_traces[0]['name'] = path_to_cif.split('/')[-1].split('.')[0] #let add a name to the trace
# frame_traces = mgr.update_framework_mesh_with_hoverdata(frame_traces, mat_updated)
# %%
co2_mesh_trace = mgr.plot_snapshot_raspa_pdb(pdb_file, index=8, natoms=3, trace_name='CO<sub>2</sub>')

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

#%%
import networkx as nx 
from ase import Atoms
# nx.cycles.cycle_basis(mat_graph)
# sg = mat_graph.subgraph(nx.cycles.cycle_basis(mat_graph)[402])
subset =[250, 269, 221, 483, 464, 465]
sg = mat_graph.subgraph(subset)
# sg= nx.relabel.convert_node_labels_to_integers(sg, 0)
cell = mat_updated.get_cell()

#make atoms object from the subgraph
atoms = Atoms(symbols=list(nx.get_node_attributes(sg, 'element').values()), positions = list(nx.get_node_attributes(sg, 'position').values()), cell=cell)
atoms_ = mgr.update_with_covalent_bonds(atoms)
atoms_.info['atom ids'] = list(nx.get_node_attributes(sg, 'element').keys())
# bonds = [list(e) for e in sg.edges()]
# atoms.info['covalent_bonds'] = bonds
atoms_trace = mgr.plot_framework_as_mesh(atoms_, show_rings=True, stylize =True)
fig = go.Figure(data=[*atoms_trace, co2_mesh_trace]).update(layout=mgr.white_layout).update_layout(legend_groupclick='togglegroup')
fig.show()

# %%
from ase.data.colors import cpk_colors
print(cpk_colors)
# %%
def plot_subset_of_atoms(ase_atoms, subset, **kwargs):
    import networkx as nx
    mat_graph = mgr.sgraph_from_atoms(ase_atoms)
    sg = mat_graph.subgraph(subset)
    cell = ase_atoms.get_cell()
    #make atoms object from the subgraph
    atoms = Atoms(symbols=list(nx.get_node_attributes(sg, 'element').values()), positions = list(nx.get_node_attributes(sg, 'position').values()), cell=cell)
    atoms_ = mgr.update_with_covalent_bonds(atoms)
    atoms_.info['atom ids'] = list(nx.get_node_attributes(sg, 'element').keys())
    # bonds = [list(e) for e in sg.edges()]
    # atoms.info['covalent_bonds'] = bonds

    # lets put some defualt values in for the kwargs
    stylize = kwargs.get('stylize', True)
    show_rings = kwargs.get('show_rings', True)
    trace_name = kwargs.get('trace_name', 'subset')

    kwargs.pop('stylize', None)
    kwargs.pop('show_rings', None)
    kwargs.pop('trace_name', None)

    atoms_trace = mgr.plot_framework_as_mesh(atoms_, show_rings=show_rings, stylize =stylize, trace_name=trace_name, **kwargs)    

    return atoms_trace

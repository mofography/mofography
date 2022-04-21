#%%
from ast import Str
import sys
import os
import numpy as np 

import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default ='browser'
path_to_mgr = os.path.abspath('c:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
sys.path.append(path_to_mgr)

#%%
from ase.io import read

nu = read('./NU-1000_S1.cif')


from ase.geometry import complete_cell
A = complete_cell(nu.get_cell())
# %%
pos = nu.get_positions()
sc_pos = nu.get_scaled_positions()

pos_2 = np.dot(A.T, sc_pos.T)
# %%
np.max(np.abs(pos-pos_2.T))
# %%
np.array([[2,0,0],[0,2,0],[0,0,2]])

#%%
np.where(np.array(nu.get_chemical_symbols()) == 'C')
# %%
from pymatgen.core import Structure, Molecule
nu_2 = Structure.from_file('./NU-1000_S1.cif')
# nu_3 = Molecule.from_file('./NU-1000_S1.cif')

#%%
nu_2.frac_coords
# %%
from pymatgen.io.ase import AseAtomsAdaptor 

AseAtomsAdaptor.get_structure(nu)
# %%
AseAtomsAdaptor.get_atoms(nu_2)
# %%
import mofography as mgr
go.Figure(data=[*mgr.draw_framework(nu), mgr.plot_cell(nu)], layout=mgr.white_layout)
# %%
from pymatgen.analysis.local_env import CrystalNN, JmolNN
from pymatgen.analysis.graphs import StructureGraph

cnn = CrystalNN(weighted_cn=True, distance_cutoffs=(0.3,0.8))
jnn = JmolNN()

nu_graph = StructureGraph.with_local_env_strategy(nu_2, strategy=cnn)
#%%
nu_graph_2 = StructureGraph.with_local_env_strategy(nu_2, strategy=jnn)

#%%
graph = nu_graph.graph
graph.nodes(data=True)

# %%
graph.edges(data=True) 
# %%
nu_graph.structure.frac_coords
# %%
import mofography as mgr
bonds, bond_orders = mgr.find_covalent_bonds_nl(nu, list(nu_graph.graph.edges()))
go.Figure(data=[*mgr.draw_framework2(nu, bonds= bonds, bond_data={'bond_orders':bond_orders}), mgr.plot_cell(nu)], layout=mgr.dark_layout)
# %%
import networkx as nx 
cycles = nx.cycle_basis(graph, 0)
# %%
mgr.draw_framework2(nu, bonds= bonds, bond_data={'bond_orders':bond_orders})
# %%
def sgraph_from_atoms(ase_atoms, bonds):
    """
    > It takes an ASE atoms object, a list of bonds, and a list of bond orders, and returns a networkx
    graph object. 
    
    The function is pretty simple, but it's a good example of how to use the ASE atoms object to get
    information about the atoms, and how to use the networkx library to create a graph object. 
    
    Let's try it out:
    
    :param ase_atoms: ASE atoms object
    :param bonds: list of tuples of atom indices that are bonded
    :param bond_data: a pandas dataframe with the bond orders of each bond
    :return: A graph object.
    """

    names = ase_atoms.get_chemical_symbols()
    import networkx as nx 
    import pandas as pd 

    G = nx.Graph()
    G.add_nodes_from([tuple((i, {'element':names[i]})) for i in range(ase_atoms.get_global_number_of_atoms())])
    # G.add_nodes_from([(i,j) for i,j in enumerate(mol.as_dict()['sites'])])
    # bond_types = pd.cut(bond_data['bond orders'], bins=[0,1, 2, 3], labels=['single', 'double', 'triple'])
    G.add_edges_from([tuple((bonds[i][0], bonds[i][1])) for i in range(len(bonds))])
    return G

#%%
sgraph = sgraph_from_atoms(nu, bonds=list(nu_graph.graph.edges()))
# %%

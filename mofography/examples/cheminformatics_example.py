#%%
import sys
import os
import numpy as np 

import plotly.graph_objects as go
import plotly.io as pio
from torch import threshold

from mofography.analyze import find_aromatic_cycles_in_sgraph
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


#%%
sgraph = sgraph_from_atoms(nu, bonds=list(nu_graph.graph.edges()))
# %%
import networkx as nx
cycles = nx.cycle_basis(sgraph, 0)
window_cycles = [c for c in cycles if len(c) in [8,12,16]]
# %%
names = np.array(nu_graph.structure.species)
ele_cycles = [names[c] for c in window_cycles]
len_cycles = [len(c) for c in window_cycles]

# %%
sum(np.array(len_cycles)==8)
# %%
grid = mgr.dgrid_from_atoms(nu, spacing =0.25, block_size=20000)
regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=0.5, threshold_abs=1.0, min_distance=4)
connections = mgr.connections_from_regions(region_labels=regions, ase_atoms=nu)
graph = mgr.graph_from_connections(connections)
#%%
cage_dia = grid.compute()[tuple(maxima.T)]
graph = mgr.add_props_to_graph(graph, node_props=[ cage_dia.tolist()], node_prop_names=['node_coords'])#, edge_props=None[[c[2] for c in connections if c[1]]], edge_prop_names=['window_coords'])

# mgr.add_props_to_graph(graph, node_props=[cage_dia], node_prop_names=['node_coords'], edge)

#%%
# node_coords = grid.compute()[maxima]
fig = go.Figure(data= [
                mgr.visualize.plot_surface_dgrid(grid.compute(), isovalue=1, region_labels = regions, ase_atoms = nu, stylize=True, cscale='Viridis'), 
                mgr.visualize.plot_cell(nu),
                mgr.scatter3d(mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=nu), name='local maxima',marker_color = grid.compute()[tuple(maxima.T)], hovertext=list(map(str, np.round(grid.compute()[tuple(maxima.T)],decimals=2))) ),
                *mgr.draw_framework(nu),
                ], layout= mgr.white_layout)# %%

# %%
fig.show()
# %%
import pymatgen as mg
from ase.io import read
# mat_atoms = read('./LTA.cif')
# mat_atoms = read('./NU-1000_S1.cif')
mat_atoms = read('./LTA96.cif')

#%%
import sys
import os
import numpy as np 

import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default ='browser'
path_to_mgr = os.path.abspath('c:/Users/Arun Gopalan/OneDrive - The University of Manchester/my_cloned_repos/mofography')
sys.path.append(path_to_mgr)
#%%
from pymatgen.io.ase import AseAtomsAdaptor 
mat_st = AseAtomsAdaptor.get_structure(mat_atoms)
#%%
# mat_st.get_all_neighbors(7)
#%%
import mofography as mgr
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
mat_graph = mgr.sgraph_from_atoms(mat_updated, use_crystalNN=False)



# * Test the aromatic detection
#%%
spacing=0.25
import mofography as mgr
grid = mgr.dgrid_from_atoms(mat_atoms, spacing =spacing, block_size=20000)
#%%
regions, maxima = mgr.regions_from_dgrid(grid.compute(), mask_thickness=0.5, exclude_border=False, min_distance=int(1/spacing), threshold_abs=2, p_norm=2)
#%%
connections = mgr.connections_from_regions(region_labels=regions, ase_atoms=mat_updated)
pgraph = mgr.graph_from_connections(connections)
#%%
bonds = mat_updated.info['covalent_bonds']
bond_orders = mat_updated.info['covalent_bond_orders']
go.Figure(data= [
                mgr.plot_surface_dgrid(grid.compute(), isovalue=1, region_labels = regions, ase_atoms = mat_updated, stylize=True, cscale=None, random_seed=123), 
                mgr.plot_cell(mat_updated),
                mgr.scatter3d(mgr.find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=mat_updated), name='local maxima',marker_color = grid.compute()[tuple(maxima.T)], hovertext=list(map(str, np.round(grid.compute()[tuple(maxima.T)],decimals=2))) ),
                *mgr.draw_framework2(mat_updated, bonds=bonds, bond_orders=bond_orders, show_aromatics=False, bond_color='grey', bond_opacity=0.5)
                ], layout= mgr.dark_layout)
# %%
hist, bins = np.histogram(grid.compute()[tuple(maxima.T)], bins=100)
binC = (bins[1:] + bins[:-1])/2
go.Figure(data = [go.Scatter(x=binC, y=hist, mode='lines', name='dhist')], layout=mgr.white_layout)
# %%
import networkx as nx 
import pysmiles as ps 
#%%
def get_smiles(cycle, graph):
    return ps.write_smiles(graph.subgraph(cycle))
#%%
sub_graphs = [mat_graph.subgraph(cycle) for cycle in nx.cycle_basis(mat_graph) if len(cycle)==8]

#%%
import itertools 
[(i.nodes(data='element'), j.nodes(data='element')) for (i,j) in list(itertools.combinations(sub_graphs, r=2))]
#%%
from networkx.algorithms import isomorphism as isomo 
em = isomo.numerical_edge_match('bond_order', 1)
nm = isomo.categorical_node_match('element','None')
iso_check = [isomo.is_isomorphic(*g,node_match=nm, edge_match=em) for g in itertools.combinations(sub_graphs, r= 2)]
#%%
[ps.get_smiles(cycle, mat_graph) for cycle in nx.simple_cycles(mat_graph)]
# %%

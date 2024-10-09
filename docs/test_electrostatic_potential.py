# %% [markdown]
# # Testing Poisson Equation Solver for Electrostatic Potential

# %% [markdown]
# ## Importing necessary libraries and modules

# %%
import numpy as np 
from pathlib import Path 
from ase.io import read
import plotly.graph_objects as go
import mofography as mgr
import pacmof as pac

# %% [markdown]
# ## Getting a list of CIFs and selecting a random one

# %% 
path_to_mgr = Path('C:/Users/arung/OneDrive - The University of Manchester/my_cloned_repos/22sep2024/mofography')
list_of_cifs = list((path_to_mgr/ 'examples'/ 'TestCIFs').glob('*.cif'))
path_to_cif = np.random.choice(list_of_cifs)
# path_to_cif = [cif for cif in list_of_cifs if 'NU-1000' in cif.name][0]
print(f"Selected CIF: {path_to_cif.name}")

# %% [markdown]
# ## Reading the CIF file and calculating charges

# %%
path_to_pacmof = Path('C:/Users/arung/OneDrive - The University of Manchester/my_cloned_repos/pacmof_clone/pacmof_oct2024/pacmof')
mat_atoms = read(path_to_cif)
mat_atoms = pac.get_charges_single_serial(path_to_cif, use_default_model=False, path_to_pickle_obj=path_to_pacmof/'pacmof' / 'data' / 'Model_RF_DDEC.pkl')

# %% [markdown]
# ## Solving the Poisson equation using spectral method

# %%
potential_V, grid_shape, cell = mgr.solve_poisson_spectral(mat_atoms, grid_spacing=0.2, epsilon_r=1.0)

# %% [markdown]
# ## Visualizing the electrostatic potential
# let's visualize the electrostatic potential as a series of isosurfaces, 
# carefully chosen to be near influential atoms in the unit cell. 
# %%
mesh_trace=mgr.plot_electrostatic_isosurfaces(potential_V, cell, mat_atoms, step=2)
# volume_trace =mgr. plot_electrostatic_volume(potential_V, cell, mat_atoms)
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
frame_traces = mgr.plot_framework_as_mesh(mat_updated, stylize=True)
cell_trace = mgr.plot_cell(mat_updated)
#%%
# fig = go.Figure(data=[*mesh_trace, *frame_traces, cell_trace], layout=mgr.dark_layout)
fig = go.Figure(data=[*mesh_trace, cell_trace, *frame_traces], layout=mgr.white_layout)
# fig = go.Figure(data=[volume_trace, cell_trace, *frame_traces], layout=mgr.white_layout)
fig.update_layout(scene_aspectmode='data')
fig.show(renderer='browser')
#%% [markdown]

#%%
from networkx import boundary_expansion
import mofography as mgr
import plotly.express as px
from pathlib import Path
from ase.io import read, write
import plotly.graph_objects as go
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
import numpy as np

#%%
path_to_cif = Path('~/repos/mofography/examples/TestCIFs/NU-1000_S1.cif').expanduser()
mat_atoms = read(path_to_cif)
#%%
spacing = 0.5
mask_thickness = 0.5
h = 0.2
extra_padding = 5.0
return_dgrid = False
return_padded_regions = False
threshold_abs = 1.5

dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing).astype(np.float32)
results = mgr.do_periodic_segmentation_rust(dgrid, mat_atoms, spacing=spacing, mask_thickness=mask_thickness, h=h, threshold_abs=threshold_abs, extra_padding=extra_padding, return_dgrid=return_dgrid)
#%%
#
regions = results['regions']
regions_pbc = results['regions_pbc']
maxima_info = results['maxima_info']
maxima_pbc = results['maxima_pbc']
all_maxima = results['all_maxima']
all_maxima_radii = dgrid[tuple(maxima_info[:, 1:4].T)]
maxima_info = results['maxima_info']
unique_regions = results['unique_regions']
close_maxima = results['close_or_inside_maxima']
pbc_groups = results['pbc_groups']
connections = results['connections']
unique_maxima = all_maxima[:regions.max()]
#%%
percolation_history, pld_locations, _ = mgr.compute_percolation_history_winding(results,expand_components=False)
#%%
print(list(results.keys()))
#%%
percolation_history[-20]
connections_trimmed = [c for c in connections if np.all(np.isin(c[0], unique_regions))]
rag = mgr.rag_from_connections(connections_trimmed, unique_maxima)
rag = mgr.add_maxima_to_rag(rag, unique_maxima, all_maxima_radii[:regions.max()], dgrid.shape, mat_atoms)

#%%
cell_trace = mgr.plot_cell(mat_atoms)
rag_trace = mgr.plot_rag(rag)
probe = float(percolation_history[-1][0])  # scalar
probe = 3.25
isovalue = min(probe - 0.25, 1.75)
percolated_cluster_meshes = mgr.plot_percolated_cluster_meshes(
    percolation_history,
    dgrid,
    regions,
    mat_atoms,
    isovalue=isovalue,
    stylize=True,
    probe_radius=probe,
    cscale = px.colors.qualitative.Pastel,
    boundary_line_width=1,
    boundary_line_color='black'
)
fig = go.Figure(data=[cell_trace, rag_trace, *percolated_cluster_meshes], layout=mgr.epic_layout)
# fig = mgr.make_fig_orthographic(fig)
# fig.update_layout(legend_groupclick='togglegroup')
# fig.show()
#%%
fig.data
#%%
fig.show(renderer='browser')
#%%
[p for p in percolation_history if p[0]==3.25]

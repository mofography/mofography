#%%
import mofography as mgr
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from ase.io import read
import plotly.graph_objects as go
import numpy as np
from ase.build import make_supercell
from skimage.measure import regionprops
import plotly.express as px

nx_cells, ny_cells, nz_cells = mgr.get_repetitions_to_meet_cutoff(mat_atoms, 6)
print(nx_cells, ny_cells, nz_cells)
mat_atoms = make_supercell(mat_atoms, [[nx_cells, 0, 0], [0, ny_cells, 0], [0, 0, nz_cells]])

# dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing).astype(np.float32)
dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=spacing)
results = mgr.do_periodic_segmentation_rust(dgrid, mat_atoms, mask_thickness=mask_thickness, h=h,spacing=spacing, return_dgrid =True,  threshold_abs=threshold_abs, min_distance=min_distance)
#%%
# grab separate variables from output
dgrid = results['dgrid']
regions = results['regions']
regions_pbc = results['regions_pbc']
maxima_info = results['maxima_info']
pbc_groups = results['pbc_groups']
connections = results['connections']
all_maxima = results['all_maxima']
all_maxima_radii = dgrid[tuple(maxima_info[:,1:4].T)]
unique_regions = results['unique_regions']
maxima = all_maxima[:regions.max()]
maxima_radii = all_maxima_radii[:regions.max()]
connections_trimmed = [c for c in connections if np.all(np.isin(c[0], unique_regions))]
connections_traces = mgr.plot_connections_as_it_is(connections_trimmed,all_maxima, dgrid.shape, mat_atoms)
rag = mgr.rag_from_connections(connections_trimmed,maxima)

rag = mgr.add_maxima_to_rag(rag,maxima,maxima_radii,dgrid.shape, mat_atoms)
rag = mgr.add_pbc_to_rag(rag, pbc_groups)
props = regionprops(regions)
bins = np.linspace(0,50,30)
rag = mgr.add_psd_to_rag_with_pbc(rag, props, pbc_groups, bins, mat_atoms, dgrid.shape)
rag = mgr.add_clustering_to_rag_geometry_js(rag, props, pbc_groups, mat_atoms, dgrid.shape)
regions_cluster = mgr.clustering_grid_from_regions(regions, rag, cluster_key='cluster_geom_pbc')
#%%
# make figure using some or all of teh dta
isoval = 2.0
cscale = px.colors.qualitative.Pastel
cell_trace = mgr.plot_cell(mat_atoms)
surface_trace_regions = mgr.plot_surface_dgrid(dgrid, isoval, mat_atoms, regions, select='all', cscale = cscale, stylize=True, name = 'Regions')
surface_trace_regions_pbc = mgr.plot_surface_dgrid(dgrid, isoval, mat_atoms,regions_pbc, select='all', cscale = cscale, stylize=True, name='Regions PBC')
surface_trace_clusters= mgr.plot_surface_dgrid(dgrid, isoval, mat_atoms,regions_cluster, select='all', cscale = cscale, stylize=True, name='Regions types')
windows_trace = mgr.plot_windows(connections_trimmed,marker_color='purple')
maxima_trace = mgr.plot_maxima(dgrid, maxima, mat_atoms, )
frame_trace = mgr.plot_framework_as_mesh(mat_atoms, stylize=True)
fig = go.Figure(data=[*frame_trace, cell_trace, surface_trace_regions, surface_trace_regions_pbc, surface_trace_clusters, windows_trace, maxima_trace, *connections_traces], layout=mgr.white_layout)
# fig.write_image('simple_test.svg')

fig.show(renderer='firefox')

#%%

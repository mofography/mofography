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

#%%
path_to_cifs = Path('~/repos/mofography/examples/TestCIFs').expanduser()
material = 'Case2-combined-revised'
# material = 'Case2-combined'
path_to_cif = path_to_cifs / f'{material}.cif'
print(f"{material} selected CIF: {path_to_cif.exists()}")
# %%
spacing = 0.5
mask_thickness = 0.5
h = 0.2
threshold_abs = 1.5
extra_padding = 5

mat_atoms = read(path_to_cif)
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=0.5).astype(np.float32)
results = mgr.do_periodic_segmentation_rust(dgrid, mat_atoms, mask_thickness=mask_thickness, h=h, spacing=spacing, return_dgrid=True, threshold_abs=threshold_abs, min_distance=1.5, return_padded_regions=True, extra_padding=extra_padding, max_smooth_radius=2.0, preserve_threshold=3.0, smooth_threshold=0.0)

# %%
regions = results['regions']
regions_pbc = results['regions_pbc']
maxima_pbc = results['maxima_pbc']
maxima_info = results['maxima_info']
connections = results['connections']
all_maxima = results['all_maxima']
unique_regions = results['unique_regions']
close_or_inside_regions = results['close_or_inside_regions']
connections_trimmed = [c for c in connections if np.any(np.isin(c[0], close_or_inside_regions))]
rag = mgr.rag_from_connections(connections_trimmed, all_maxima)
# %%
# let's do some plotting 
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)
connections_trace = mgr.plot_connections_as_it_is(connections_trimmed, all_maxima, dgrid.shape, mat_atoms, spacing=spacing)
maxima_trace = mgr.plot_maxima(dgrid, maxima_pbc, mat_atoms, show_spheres=False, sphere_opacity=0.5)
cscale = px.colors.qualitative.Pastel
surface_trace = mgr.plot_surface_dgrid(dgrid, isovalue=1.0, region_labels=regions, ase_atoms=mat_atoms, cscale=cscale, stylize=True)
fig = go.Figure(data=[*connections_trace, surface_trace, maxima_trace,  cell_trace], layout=mgr.dark_layout)
fig.show()
# %%
mat_atoms.cell.cellpar()

# %%

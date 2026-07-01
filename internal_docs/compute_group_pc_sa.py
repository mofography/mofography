# compute the surface area for a pbc group
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
from ase.units import _amu

# lload cif and compute mass
path_to_cif=Path('/home/arung014/Downloads/tobmof-6.cif')
mat_atoms=read(str(path_to_cif))
mass_g = mat_atoms.get_masses().sum() * _amu * 1000  # this is in grams

# compute dgrid and segmentaiton
dgrid = mgr.dgrid_from_atoms(mat_atoms, spacing=0.5, block_size=50000).compute()
probe_size = vdw_radii[2]
extra_padding = 10
results= mgr.do_periodic_segmentation(dgrid,mat_atoms= mat_atoms,  wall_thickness=0.08, extra_padding=extra_padding)
regions_pbc = results['regions_pbc']

# compute the isosurface mesh (whole and pieces) and compute areas in m2pg
dist_iso_surf = mgr.get_isosurface_mesh(dgrid, probe_size, mat_atoms)
# dist_iso_submeshes = mgr.get_isosurface_submeshes(dgrid, probe_size, mat_atoms, regions_pbc)
dist_iso_submeshes =mgr.get_isosurface_submeshes(dgrid,probe_size,mat_atoms, regions_pbc, colorscale=None) # no cscale defaults to px.colors.qualitative.Pastel
group_areas = [m.area * 1e-20 / mass_g for m in dist_iso_submeshes]
print('Total number of faces', len(dist_iso_surf.faces))
print('Area sum in m2pg at {} \u212b'.format(probe_size),np.sum(group_areas))
print('Area of whole mesh in m2pg at {} \u212b'.format(probe_size), dist_iso_surf.area * 1e-20 / mass_g)

# refctored function

# let's plot and check
cscale = px.colors.qualitative.Pastel
cell_trace = mgr.plot_cell(mat_atoms)
surface_trace =  mgr.plot_surface_dgrid(
    dgrid,
    isovalue=probe_size,
    region_labels=regions_pbc,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name='Regions at ' + str(probe_size) + '\u212b',
    stylize=True,
)
pbc_group_mesh= mgr.plot_mesh([dist_iso_submeshes[0]], colorscale=cscale, name='pbc_group first', stylize=True)
data = [cell_trace, surface_trace, meshes_traces]
fig = go.Figure(data=data, layout=mgr.dark_layout)
fig = mgr.make_fig_orthographic(fig)
fig.write_html('temp.html')


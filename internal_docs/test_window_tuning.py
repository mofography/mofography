#%%
#+BEGIN_SRC python :session textual_properties_calc :results output
import mofography as mgr
import copy
from mofography.winding_percolation import find_pore_limiting_radius_winding
import numpy as np
from pathlib import Path
import pickle
from dask import delayed
import dask
from dask.distributed import Client, LocalCluster, progress
import logging
import stackprinter
stackprinter.set_excepthook(style='darkbg2')
import pandas as pd
logging.basicConfig(level=logging.INFO)
from tqdm import tqdm
import plotly.graph_objects as go
from scipy import stats
from ase.data import vdw_radii, vdw_alvarez
from ase.io import read
print("imports done")
#+END_SRC
#+BEGIN_SRC python :session textual_properties_calc :results output
# path to the pickled IZA atoms list
path_to_atoms = Path('~/Datasets/iza_atoms.pkl').expanduser()
print("exists:", path_to_atoms.exists())
iza_atoms = pickle.load(open(path_to_atoms, 'rb'))
# iza_atoms is expected to be a list-like where some entries might be 'NaN'
print("loaded iza_atoms, length:", len(iza_atoms))
#+END_SRC
#+BEGIN_SRC python :session textual_properties_calc
spacing = 0.25
radii_source = 'C&P Handbook'   # options used in your script
dgrid_method = 'python AABB'
radii = None
#+END_SRC
#+begin_src python :session textual_properties_calc :results output

# output_iza_results_path = Path('~/repos/mofography_manuscript/resources/mgr_org/mgr_org_scripts/iza_results_0.5spc_0.5mt_0.2h.pkl').expanduser()
output_iza_results_path = Path('~/repos/mofography_manuscript/resources/mgr_org/mgr_org_scripts/iza_results_0.25spc_0.5mt_0.2h_aabb.pkl').expanduser()
# output_iza_results_path = Path('~/repos/mofography_manuscript/resources/mgr_org/mgr_org_scripts/iza_results_0.25spc_0.5mt_0.2h_rust.pkl').expanduser()
iza_results = pickle.load(open(output_iza_results_path, 'rb'))
print('Done')

#+end_src
 #+begin_src python
# FInd the index
name ='SGT-0'
# name ='MFI-0'
mgr_index = [i for i, atoms in enumerate(iza_atoms) if atoms != 'NaN' if atoms.info['name']==name][ 0]
print(mgr_index)
print(iza_atoms[mgr_index].info['name'])
mat_atoms = iza_atoms[mgr_index]
# compute results
spacing = 0.5
h= 0.2
mask_thickness= 0.5
threshold_abs= 1.0
extra_padding= 5
return_dgrid = True
min_distance= 1.5
from ase.build import make_supercell
nx_cells, ny_cells, nz_cells = mgr.get_repetitions_to_meet_cutoff(mat_atoms, 6)
mat_atoms = make_supercell(mat_atoms, [[nx_cells, 0, 0], [0, ny_cells, 0], [0, 0, nz_cells]])
dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=spacing)
# dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=20000).astype(np.float32)
# dgrid = mgr.dgrid_from_atoms(mat_atoms, spacing=spacing, block_size=50000).compute().astype(np.float32)

results = mgr.do_periodic_segmentation_rust(
    dgrid,
    mat_atoms,
    mask_thickness=mask_thickness,
    h=h,
    spacing=spacing,
    return_dgrid=True,
    threshold_abs=threshold_abs,
    min_distance=min_distance,
    return_padded_regions=True,
    extra_padding=extra_padding,
    tune_connections=True,
)

# let's plot using the old saved results
# results = iza_results[mgr_index]

from skimage.measure import regionprops
# grab separate variables from output
dgrid = results["dgrid"]
regions = results["regions"]
regions_pbc = results["regions_pbc"]
maxima_info = results["maxima_info"]
pbc_groups = results["pbc_groups"]
connections = results["connections"]
all_maxima = results["all_maxima"]
all_maxima_radii = dgrid[tuple(maxima_info[:, 1:4].T)]
maxima_pbc = results["maxima_pbc"]
unique_regions = results["unique_regions"]
maxima = all_maxima[: regions.max()]
maxima_radii = all_maxima_radii[: regions.max()]
connections_trimmed = [c for c in connections if np.any(np.isin(c[0], unique_regions))]
# connections_trimmed = [c for c in connections if 5 in c[0]]
maxima_pbc_dict = results['maxima_pbc_dict']
connections_traces = mgr.plot_connections_as_it_is(
    connections_trimmed, all_maxima, dgrid.shape, mat_atoms
)
# regions_padded = results["regions_padded"]
rag = mgr.rag_from_connections(connections_trimmed, maxima)
rag = mgr.add_maxima_to_rag(rag, maxima, maxima_radii, dgrid.shape, mat_atoms)
rag = mgr.add_pbc_to_rag(rag, pbc_groups)
props = regionprops(regions)
bins = np.linspace(0, 50, 30)
rag = mgr.add_psd_to_rag_with_pbc(rag, props, pbc_groups, bins, mat_atoms, dgrid.shape)
rag = mgr.add_clustering_to_rag_geometry_js(
    rag, props, pbc_groups, mat_atoms, dgrid.shape
)
regions_cluster = mgr.clustering_grid_from_regions(
    regions, rag, cluster_key="cluster_geom_pbc"
)
#%%
import plotly.express as px
import plotly.graph_objects as go
cscale = px.colors.qualitative.Pastel
isovals = [0.65,1.5 ]

# mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
# cscale = px.colors.qualitative.Alphabet
isosurface_segmentation = mgr.plot_surface_dgrid(dgrid, stylize=True, isovalue=isovals[0], region_labels=regions, ase_atoms=mat_atoms, name='Regions', opacity=1, cscale=cscale)
isosurface_pbc = mgr.plot_surface_dgrid(dgrid, stylize=True, isovalue=isovals[0], region_labels=regions_pbc, ase_atoms=mat_atoms, name='Regions PBC',  opacity=1, cscale=cscale )
cell_trace = mgr.plot_cell(mat_atoms, color='orange')
maxima_trace = mgr.plot_maxima(dgrid,maxima, mat_atoms, marker_size=10,show_spheres=False, marker_colorscale = cscale, name='Maxima (distinct)')
# frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)
# frame_trace = mgr.draw_framework(mat_updated, stylize=True)

# do the same thing with energy?
# egrid = mgr.egrid_from_atoms2(mat_atoms, spacing=0.5, probe_symbol='He').compute()
# einterpolator = mgr.interpolate_dgrid(egrid) # we can use the same interpolator as that for the distance grid, because at the end of the day, it's just a 3D array and we are interpolating in fractional coordinates

# rag = mgr.add_clustering_to_rag_energy_js(rag, regions, energy_interpolator=einterpolator, js_threshold=0.10, pbc_groups = pbc_groups) # no need for ase_atoms, because we already have the energy interpolator
# regions_clusters = mgr.clustering_grid_from_regions(regions, rag, cluster_key='cluster_energy_pbc')

# isosurface_clusters = mgr.plot_surface_dgrid(dgrid, stylize=True, isovalue=isovals[0], region_labels=regions_clusters, name='Pocket Types by geometry', ase_atoms=mat_atoms, opacity=1, cscale=cscale)
maxima_trace_pbc= mgr.plot_maxima(dgrid,maxima_pbc, mat_atoms,show_spheres=False,  marker_size=7.5, marker_colorscale = cscale,name='Maxima (pbc)')
connections_traces =mgr.plot_connections_as_it_is(connections_trimmed, all_maxima, regions.shape, mat_atoms, cscale=cscale, plot_disconnected=False, plot_periodic=True)
connections_traces [-2].marker.size = 1.5
windows_trace = mgr.plot_windows(connections_trimmed, name='Windows', marker_color='yellow', marker_size=10, marker_opacity=0.8)

fig = go.Figure([cell_trace, *connections_traces, isosurface_pbc, isosurface_segmentation, windows_trace, maxima_trace, maxima_trace_pbc], layout=mgr.white_layout)
fig.update_layout(legend=dict(groupclick='toggleitem'), title=dict(text=iza_atoms[mgr_index].info['name']))
fig.show()
 #+end_src
#%%
mgr.find_pore_limiting_radius2(dgrid), mgr.find_pore_limiting_radius_winding(results)

# %%
import mofography as mgr
from pathlib import Path
from ase.io import read
import plotly.graph_objects as go
import numpy as np
import time
# path_to_cif = Path('~/repos/mofography/examples/TestCIFs/NU-1000_S1.cif').expanduser()
# path_to_cif = Path('~/repos/mofography/examples/TestCIFs/NU-1000_S1.cif').expanduser()
path_to_cif = Path('~/repos/mofography_manuscript/resources/mgr_org/mgr_org_scripts/pim-1r9_5.cif').expanduser()
spacing = 0.25
mat_atoms = read(path_to_cif)
start_time = time.time()
dgrid_cpu = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=20000)
cpu_time = time.time() - start_time
start_time = time.time()
dgrid_rust_aabb = mgr.dgrid_from_atoms_aabb_rust(mat_atoms, spacing=spacing)
aabb_time = time.time() - start_time
start_time = time.time()
dgrid_rust = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=spacing)
rust_time = time.time() - start_time
print(f"CPU time: {cpu_time}, AABB time: {aabb_time}, Rust time: {rust_time}")
# Plot absolute deviation differences between the three distance fields as bar plots in stacked subplots (3 rows, 1 column)
from plotly.subplots import make_subplots

# Compute flattened absolute deviations
dev_aabb = np.abs(dgrid_cpu - dgrid_rust_aabb).flatten()
dev_rust = np.abs(dgrid_cpu - dgrid_rust).flatten()
dev_aabb_rust = np.abs(dgrid_rust_aabb - dgrid_rust).flatten()

fig = make_subplots(
    rows=3, cols=1,
    shared_xaxes=False,
    subplot_titles=[
        "Abs. deviation: CPU vs AABB",
        "Abs. deviation: CPU vs Rust",
        "Abs. deviation: AABB vs Rust"
    ]
)

# Compute histogram data for manual bar plotting
def absdev_histogram(absdev, nbins=50):
    hist, bin_edges = np.histogram(absdev, bins=nbins)
    # bar centers
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return bin_centers, hist

nbins = 50
bin_centers1, hist1 = absdev_histogram(dev_aabb, nbins)
bin_centers2, hist2 = absdev_histogram(dev_rust, nbins)
bin_centers3, hist3 = absdev_histogram(dev_aabb_rust, nbins)

fig.add_trace(
    go.Bar(x=bin_centers1, y=hist1, name="CPU vs AABB"),
    row=1, col=1
)
fig.add_trace(
    go.Bar(x=bin_centers2, y=hist2, name="CPU vs Rust"),
    row=2, col=1
)
fig.add_trace(
    go.Bar(x=bin_centers3, y=hist3, name="AABB vs Rust"),
    row=3, col=1
)

fig.update_layout(
    height=900,
    showlegend=False,
    xaxis_title="Absolute Deviation [Å]",
    yaxis_title="Count",
    xaxis2_title="Absolute Deviation [Å]",
    yaxis2_title="Count",
    xaxis3_title="Absolute Deviation [Å]",
    yaxis3_title="Count",
    title="Comparison of Absolute Deviation Differences (Stacked Bar Plots)"
)

fig.show()
# %%
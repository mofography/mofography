"""
Diagnostic script to confirm whether watershed produces background (label 0) voxels
in the center of regions - i.e. void voxels that should have been labeled but weren't.

Uses the exact params from watershed_bucket.md and Case4-combined-revised.cif.
"""
import os
import numpy as np
import pandas as pd
from ase.io import read
from pathlib import Path
import plotly.graph_objects as go

import mofography as mgr
from mofography.helper import find_coord_from_indices
from mofography.visualize import (
    plot_surface_dgrid,
    plot_cell,
    scatter3d,
    plot_maxima,
    plot_connections_as_it_is,
    dark_layout,
)

# Match watershed_bucket.md params
os.environ.setdefault("RAYON_NUM_THREADS", "8")

spacing = 0.1
# mask_thickness: 0 = original bug report; 0.5 = test if bug persists with thicker mask
mask_thickness = 0
h = 0.3
threshold_abs = 0.605
preserve_threshold = 103.0
smooth_threshold = 106.0
max_smooth_radius = 2.0
extra_padding = 2
return_dgrid = False
min_distance = 1.0
min_region_size = 0.0001
minimum_window_separation = 1.5

# Load CIF from watershed_bucket_bug
# script_dir = Path(__file__).resolve().parent
script_dir = Path(__file__).resolve().parent
cif_path = script_dir / "Case4-combined-revised.cif"
mat_atoms = read(str(cif_path))

radius = pd.Series({"C": 172})
radii = (radius[mat_atoms.get_chemical_symbols()].values / 100.0)
dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=spacing, radii=radii)

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
    preserve_threshold=preserve_threshold,
    smooth_threshold=smooth_threshold,
    max_smooth_radius=max_smooth_radius,
    extra_padding=extra_padding,
    compute_connections=True,
    min_region_size=min_region_size,
    minimum_window_separation=minimum_window_separation,
)

regions = results["regions"]
dgrid = results["dgrid"]
void_mask = dgrid > mask_thickness

# Background = label 0. Bug: any void voxel (dgrid > mask_thickness) with label 0
# shouldn't exist — watershed should assign a region to all masked (void) voxels.
n_void = np.sum(void_mask)
n_background_in_void = np.sum((regions == 0) & void_mask)
n_background_total = np.sum(regions == 0)

print("=" * 60)
print("Watershed background voxel diagnostic")
print("=" * 60)
print(f"mask_thickness: {mask_thickness}")
print(f"Total void voxels (dgrid > {mask_thickness}): {n_void}")
print(f"Total background (label 0) voxels: {n_background_total}")
print(f"Background voxels IN void (BUG if > 0): {n_background_in_void}")
print(f"Fraction of void with label 0: {n_background_in_void/n_void*100:.4f}%")
print()

# Of those, "interior" = have at least one labeled neighbor (clearly inside a pore)
def has_labeled_neighbor(regions, i, j, k):
    """Check if voxel (i,j,k) has at least one 26-neighbor with non-zero label."""
    nx, ny, nz = regions.shape
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                if di == 0 and dj == 0 and dk == 0:
                    continue
                ni, nj, nk = i + di, j + dj, k + dk
                if 0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz:
                    if regions[ni, nj, nk] != 0:
                        return True
    return False


background_in_void = (regions == 0) & void_mask
interior_background = np.zeros_like(regions, dtype=bool)

indices = np.argwhere(background_in_void)
interior_indices = []
for idx in indices:
    i, j, k = idx[0], idx[1], idx[2]
    if has_labeled_neighbor(regions, i, j, k):
        interior_background[i, j, k] = True
        interior_indices.append(idx)

interior_indices = np.array(interior_indices) if interior_indices else np.empty((0, 3))
n_interior_background = len(interior_indices)

print("Background voxels in void that have at least one labeled neighbor")
print("(i.e. 'interior' to a region - the bug):")
print(f"  Count: {n_interior_background}")
print()

# Spatial distribution: PBC pairwise distance and clustering to confirm two-end clusters
if n_interior_background > 0:
    import ase.geometry
    cell = ase.geometry.complete_cell(mat_atoms.get_cell()).T  # columns = cell vectors
    shape = np.array(regions.shape, dtype=float)
    frac = interior_indices.astype(float) / shape

    def pbc_dist_cart(f1, f2):
        d_frac = (f1 - f2) - np.round(f1 - f2)
        cart = d_frac @ cell
        return np.linalg.norm(cart)

    # Max pairwise PBC distance
    n_int = len(interior_indices)
    max_pbc_dist = 0.0
    for i in range(n_int):
        for j in range(i + 1, n_int):
            d = pbc_dist_cart(frac[i], frac[j])
            if d > max_pbc_dist:
                max_pbc_dist = d
    cell_diag = np.linalg.norm(np.sum(cell, axis=0))
    print("Spatial distribution (PBC):")
    print(f"  Max pairwise PBC distance: {max_pbc_dist:.3f} Å")
    print(f"  Cell diagonal (sum of vectors): {cell_diag:.3f} Å")

    # Agglomerative clustering with PBC (2 clusters), single linkage
    from scipy.cluster.hierarchy import linkage, fcluster
    n_pairs = n_int * (n_int - 1) // 2
    condensed = np.zeros(n_pairs)
    idx = 0
    for i in range(n_int):
        for j in range(i + 1, n_int):
            condensed[idx] = pbc_dist_cart(frac[i], frac[j])
            idx += 1
    Z = linkage(condensed, method="single")
    labels_2 = fcluster(Z, t=2, criterion="maxclust")
    cluster_sizes = [np.sum(labels_2 == k) for k in np.unique(labels_2)]
    mean_frac_y = [np.mean(frac[labels_2 == k], axis=0)[1] for k in np.unique(labels_2)]
    print(f"  Cluster sizes (n_clusters=2): {cluster_sizes}")
    print(f"  Mean fractional y per cluster: {[f'{y:.3f}' for y in mean_frac_y]}")
    print()

if n_interior_background > 0:
    print("*** BUG CONFIRMED: Background voxels exist in the center of regions ***")
    # Get Cartesian coords for plotting (subsample for manageable HTML)
    max_plot_points = 120
    interior_coords = find_coord_from_indices(interior_indices, regions.shape, mat_atoms)
    if len(interior_coords) > max_plot_points:
        step = max(1, len(interior_coords) // max_plot_points)
        interior_coords = interior_coords[::step][:max_plot_points]
        print(f"  (plotting {len(interior_coords)} of {n_interior_background} interior points)")

    # Build Plotly figure
    import plotly.graph_objects as go

    traces = []

    # Isosurface at mask_thickness + 0.25, colored by regions
    isovalue = mask_thickness + 0.25
    surf_trace = plot_surface_dgrid(
        dgrid,
        isovalue,
        mat_atoms,
        region_labels=regions,
        select="all",
        opacity=0.7,
    )
    traces.append(surf_trace)

    # Cell
    cell_trace = plot_cell(mat_atoms, color="gray")
    traces.append(cell_trace)

    # Maxima and pore graph (unique regions only)
    unique_regions = results["unique_regions"]
    all_maxima = results["all_maxima"]
    unique_maxima = all_maxima[: regions.max()]  # reindexed priority
    connections = results["connections"]
    conns_unique = [c for c in connections if np.all(np.isin(c[0], unique_regions))]

    maxima_traces = plot_maxima(dgrid, unique_maxima, mat_atoms, name="Maxima")
    if isinstance(maxima_traces, (list, tuple)):
        traces.extend(maxima_traces)
    else:
        traces.append(maxima_traces)

    conn_traces = plot_connections_as_it_is(
        conns_unique,
        unique_maxima,
        regions.shape,
        mat_atoms,
        plot_periodic=True,
    )
    if isinstance(conn_traces, (list, tuple)):
        traces.extend(conn_traces)
    else:
        traces.append(conn_traces)

    # Interior background voxels (the bug) - visible on dark background with yellow outline
    interior_trace = scatter3d(
        interior_coords,
        name="Background voxels interior to regions (BUG)",
        marker_color="white",
        marker_size=8,
        marker_line_color="yellow",
        marker_line_width=2,
    )
    traces.append(interior_trace)

    # Boundary background: subsample if many
    interior_set = {tuple(idx) for idx in interior_indices}
    all_bg_indices = np.argwhere(background_in_void)
    boundary_indices = np.array(
        [idx for idx in all_bg_indices if tuple(idx) not in interior_set]
    )
    if len(boundary_indices) > 0:
        boundary_coords = find_coord_from_indices(
            boundary_indices, regions.shape, mat_atoms
        )
        if len(boundary_coords) > max_plot_points:
            step = max(1, len(boundary_coords) // max_plot_points)
            boundary_coords = boundary_coords[::step][:max_plot_points]
        boundary_trace = scatter3d(
            boundary_coords,
            name="Background at boundary (expected)",
            marker_color="rgba(200,200,200,0.8)",
            marker_size=5,
            marker_line_color="yellow",
            marker_line_width=1,
        )
        traces.append(boundary_trace)

    fig = go.Figure(data=traces)
    fig.update_layout(
        **dark_layout,
    )
    fig.update_layout(
        title="Watershed bug: background (label 0) voxels in center of regions",
        scene=dict(
            aspectmode="data",
            xaxis_title="X (Å)",
            yaxis_title="Y (Å)",
            zaxis_title="Z (Å)",
        ),
        width=900,
        height=700,
    )

    out_path = script_dir / "background_voxels_visualization.html"
    fig.write_html(str(out_path))
    print(f"Plot saved to: {out_path}")
    print("Open in browser to inspect. White markers (yellow outline) = interior background voxels (bug).")
    # Spatial distribution: plotting uses the same grid→Cartesian transform as the rest of the figure,
    # so positions are correct. If you see points only at the channel ends, it can be view/occlusion;
    # the fractional coords below show the full extent (they often span the long axis).
    interior_frac = interior_indices.astype(float) / np.array(regions.shape)
    print("Interior bug voxels: fractional coords min [i,j,k] / max [i,j,k]:", interior_frac.min(axis=0), "/", interior_frac.max(axis=0))
else:
    print("No interior background voxels found - bug not reproduced with these params.")

print()
print("For comparison with regions_from_dgrid (old method):")
try:
    regions_old, maxima_old = mgr.regions_from_dgrid(dgrid, mask_thickness, h)
    print(f"  regions (watershed): np.sum(regions==0) = {np.sum(regions == 0)}")
    print(f"  regions_old:        np.sum(regions_old==0) = {np.sum(regions_old == 0)}")
except Exception as e:
    print(f"  (skipped: {e})")

# Same run with mask_thickness=0.5 to see if bug persists
print()
print("=" * 60)
print("Re-run with mask_thickness=0.5")
print("=" * 60)
results_05 = mgr.do_periodic_segmentation_rust(
    dgrid,
    mat_atoms,
    mask_thickness=0.5,
    h=h,
    spacing=spacing,
    return_dgrid=False,
    threshold_abs=threshold_abs,
    min_distance=min_distance,
    return_padded_regions=False,
    preserve_threshold=preserve_threshold,
    smooth_threshold=smooth_threshold,
    max_smooth_radius=max_smooth_radius,
    extra_padding=extra_padding,
    compute_connections=True,
    min_region_size=min_region_size,
    minimum_window_separation=minimum_window_separation,
)
regions_05 = results_05["regions"]
void_mask_05 = dgrid > 0.5
n_void_05 = np.sum(void_mask_05)
n_bg_in_void_05 = np.sum((regions_05 == 0) & void_mask_05)
print(f"Total void voxels (dgrid > 0.5): {n_void_05}")
print(f"Background voxels IN void (BUG if > 0): {n_bg_in_void_05}")
print(f"Fraction of void with label 0: {n_bg_in_void_05/n_void_05*100:.4f}%")
if n_bg_in_void_05 > 0:
    print("*** Same bug with mask_thickness=0.5 ***")
else:
    print("No background-in-void with mask_thickness=0.5.")

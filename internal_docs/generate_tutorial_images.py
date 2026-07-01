#!/usr/bin/env python3
"""
Generate images for the Periodic Segmentation Quick Start tutorial.

This script creates PNG visualizations of IRMOF-1 using the periodic segmentation workflow.
Images are saved with white backgrounds for documentation use.
"""

import os
import sys
import numpy as np
from pathlib import Path

# Setup paths
script_path = Path(__file__).resolve()
if script_path.parent.name == 'docs':
    repo_root = script_path.parent.parent
else:
    repo_root = script_path.parent

sys.path.insert(0, str(repo_root))

import mofography as mgr
from ase.io import read
import plotly.graph_objects as go
import plotly.express as px

# =============================================================================
# Load Structure
# =============================================================================

print("Loading IRMOF-1 structure...")
path_to_cif = repo_root / 'examples' / 'TestCIFs' / 'P1-IRMOF-1.cif'
if not path_to_cif.exists():
    print(f"ERROR: CIF file not found at {path_to_cif}")
    sys.exit(1)

mat_atoms = read(str(path_to_cif))
print(f"Loaded: {mat_atoms.get_chemical_formula()}")

# =============================================================================
# Run Periodic Segmentation
# =============================================================================

print("\nRunning periodic segmentation...")

spacing = 0.5
mask_thickness = 1.0
h = 0.2
threshold_abs = 1.5
minimum_window_separation = 1.5
wall_thickness = 0.08
extra_padding = 10

# Compute distance grid
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=50000)
dgrid = mgr.zoom_grid(dgrid, zoom_factor=1).compute()
dgrid = mgr.periodic_gaussian_filter(dgrid, sigma=spacing*2)
dgrid = dgrid.astype(np.float32)

# Run periodic segmentation
results = mgr.do_periodic_segmentation_rust(
    dgrid,
    mat_atoms,
    spacing=spacing,
    mask_thickness=mask_thickness,
    h=h,
    threshold_abs=threshold_abs,
    minimum_window_separation=minimum_window_separation,
    wall_thickness=wall_thickness,
    extra_padding=extra_padding,
    return_padded_regions=True,
    compute_connections=True
)

# Extract results
maxima_pbc = results['maxima_pbc']
regions_pbc = results['regions_pbc']
regions = results['regions']
connections = results['connections']
all_maxima = results['all_maxima']
unique_regions = results['unique_regions']
pbc_groups = results['pbc_groups']

maxima = all_maxima[:regions.max()]

print(f"Found {len(maxima_pbc)} unique pores")

# =============================================================================
# Create output directory
# =============================================================================

output_dir = Path(__file__).parent / 'tutorial_images'
output_dir.mkdir(exist_ok=True)
print(f"\nSaving images to: {output_dir}")

# =============================================================================
# Generate Image 1: Split vs Merged Regions
# =============================================================================

print("\nGenerating image: split_vs_merged_regions.png")

mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
cscale = px.colors.qualitative.Pastel

cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)

# Plot split regions
surface_trace_split = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=1.5,
    region_labels=regions,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name='Split regions',
    stylize=True,
    opacity=0.8
)

fig1 = go.Figure(
    data=[cell_trace, surface_trace_split, *frame_trace],
    layout=mgr.white_layout
)

fig1.update_layout(
    scene_aspectmode='data',
    title=dict(
        text="BEFORE PBC: Split Pores at Boundaries",
        x=0.5,
        xanchor='center',
        font=dict(size=14, color='black')
    ),
    scene=dict(
        xaxis=dict(title='X (Å)'),
        yaxis=dict(title='Y (Å)'),
        zaxis=dict(title='Z (Å)')
    )
)

fig1.write_image(output_dir / 'split_regions.png', width=800, height=600, scale=2)
print("  Saved: split_regions.png")

# =============================================================================
# Generate Image 2: Merged Regions
# =============================================================================

print("Generating image: merged_regions.png")

surface_trace_merged = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=1.5,
    region_labels=regions_pbc,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name='Merged regions',
    stylize=True,
    opacity=0.8
)

fig2 = go.Figure(
    data=[cell_trace, surface_trace_merged, *frame_trace],
    layout=mgr.white_layout
)

fig2.update_layout(
    scene_aspectmode='data',
    title=dict(
        text="AFTER PBC: Correctly Merged Pores",
        x=0.5,
        xanchor='center',
        font=dict(size=14, color='black')
    ),
    scene=dict(
        xaxis=dict(title='X (Å)'),
        yaxis=dict(title='Y (Å)'),
        zaxis=dict(title='Z (Å)')
    )
)

fig2.write_image(output_dir / 'merged_regions.png', width=800, height=600, scale=2)
print("  Saved: merged_regions.png")

# =============================================================================
# Generate Image 3: With Connections
# =============================================================================

print("Generating image: with_connections.png")

connections_trimmed = np.array([
    c for c in connections
    if np.any(np.isin(c[0], unique_regions))
])

maxima_trace = mgr.plot_maxima(
    dgrid, maxima, mat_atoms,
    name='Pore centers',
    show_spheres=False,
    marker_colorscale=cscale
)

connections_traces = mgr.plot_connections_as_it_is(
    connections_trimmed,
    all_maxima,
    dgrid.shape,
    mat_atoms,
    cscale=cscale,
    plot_disconnected=False,
    show_spheres=False
)

fig3 = go.Figure(
    data=[
        cell_trace,
        surface_trace_merged,
        maxima_trace,
        *connections_traces,
        *frame_trace
    ],
    layout=mgr.white_layout
)

fig3.update_layout(
    scene_aspectmode='data',
    title=dict(
        text="IRMOF-1: Pore Network with Connections",
        x=0.5,
        xanchor='center',
        font=dict(size=14, color='black')
    ),
    scene=dict(
        xaxis=dict(title='X (Å)'),
        yaxis=dict(title='Y (Å)'),
        zaxis=dict(title='Z (Å)')
    )
)

fig3.write_image(output_dir / 'with_connections.png', width=800, height=600, scale=2)
print("  Saved: with_connections.png")


print("\n✅ All images generated successfully!")
print(f"Images saved to: {output_dir}")

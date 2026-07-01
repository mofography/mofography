#!/usr/bin/env python3
"""
Generate images for the Distance Grid tutorial.

This script creates PNG visualizations of IRMOF-1 using the distance grid workflow.
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
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
print(f"Loaded: {mat_atoms.get_chemical_formula()}")

# =============================================================================
# Compute Distance Grid using Rust Backend
# =============================================================================

print("\nComputing distance grid (Rust backend)...")

spacing = 0.5

import mgr_rust
print(f"Using {mgr_rust.rayon_num_threads()} threads")

dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=spacing)

print(f"Grid shape: {dgrid.shape}")
print(f"Grid range: {dgrid.min():.2f} to {dgrid.max():.2f} Angstroms")

# =============================================================================
# Create output directory
# =============================================================================

output_dir = Path(__file__).parent / 'tutorial_images'
output_dir.mkdir(exist_ok=True)
print(f"\nSaving images to: {output_dir}")

# =============================================================================
# Generate Image 1: Distance Grid Isosurface (2.5 Angstroms)
# =============================================================================

print("\nGenerating image: distance_grid_isosurface.png")

cscale = px.colors.qualitative.Pastel

# Create traces
cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)

# Plot isosurface at 2.5 Angstroms
isosurface_trace = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=2.5,
    ase_atoms=mat_atoms,
    color='turquoise',
    opacity=0.6,
    stylize=True,
    name='2.5 Å isosurface'
)

fig1 = go.Figure(
    data=[cell_trace, isosurface_trace, *frame_trace],
    layout=mgr.white_layout
)

fig1.update_layout(
    scene_aspectmode='data',
    title=dict(
        text="Distance Grid Isosurface (2.5 Å)",
        x=0.5,
        xanchor='center',
        font=dict(size=16, color='black')
    ),
    scene=dict(
        xaxis=dict(title='X (Å)'),
        yaxis=dict(title='Y (Å)'),
        zaxis=dict(title='Z (Å)')
    )
)

fig1.write_image(output_dir / 'distance_grid_isosurface.png', width=1200, height=900, scale=2)
print("  Saved: distance_grid_isosurface.png")

# =============================================================================
# Generate Image 2: Multiple Isosurfaces (1.5, 2.5, 3.5 Angstroms)
# =============================================================================

print("\nGenerating image: distance_grid_multiple_isosurfaces.png")

# Plot multiple isosurfaces - use discrete colors
iso_traces = []
iso_colors = ['purple', 'turquoise', 'orange']
for iso_value, color in zip([1.5, 2.5, 3.5], iso_colors):
    trace = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=iso_value,
        ase_atoms=mat_atoms,
        color=color,
        opacity=0.4,
        name=f'{iso_value} Å'
    )
    iso_traces.append(trace)

fig2 = go.Figure(
    data=[cell_trace, *iso_traces, *frame_trace],
    layout=mgr.white_layout
)

fig2.update_layout(
    scene_aspectmode='data',
    title=dict(
        text="Multiple Distance Isosurfaces (1.5, 2.5, 3.5 Å)",
        x=0.5,
        xanchor='center',
        font=dict(size=16, color='black')
    ),
    scene=dict(
        xaxis=dict(title='X (Å)'),
        yaxis=dict(title='Y (Å)'),
        zaxis=dict(title='Z (Å)')
    )
)

fig2.write_image(output_dir / 'distance_grid_multiple_isosurfaces.png', width=1200, height=900, scale=2)
print("  Saved: distance_grid_multiple_isosurfaces.png")

# =============================================================================
# Generate Image 3: Dark Layout Version (for reference)
# =============================================================================

print("\nGenerating image: distance_grid_isosurface_dark.png")

fig3 = go.Figure(
    data=[cell_trace, isosurface_trace, *frame_trace],
    layout=mgr.dark_layout
)

fig3.update_layout(
    scene_aspectmode='data',
    title=dict(
        text="Distance Grid Isosurface (2.5 Å) - Dark Theme",
        x=0.5,
        xanchor='center',
        font=dict(size=16, color='white')
    ),
    scene=dict(
        xaxis=dict(title='X (Å)'),
        yaxis=dict(title='Y (Å)'),
        zaxis=dict(title='Z (Å)')
    )
)

fig3.write_image(output_dir / 'distance_grid_isosurface_dark.png', width=1200, height=900, scale=2)
print("  Saved: distance_grid_isosurface_dark.png")


print("\n✅ All images generated successfully!")
print(f"Images saved to: {output_dir}")

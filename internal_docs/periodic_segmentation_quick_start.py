#!/usr/bin/env python3
"""
Periodic Segmentation Quick Start - Standalone Test Script

This script demonstrates MOFography's periodic segmentation workflow
for analyzing pore networks in crystalline nanoporous materials.

Run this script to test the workflow before using it in your own projects.
"""

import os
import sys
import numpy as np
from pathlib import Path

# =============================================================================
# STEP 1: Setup and Imports
# =============================================================================

print("=" * 60)
print("MOFography - Periodic Segmentation Quick Start")
print("=" * 60)

# Get the script's directory and find the repository root
script_path = Path(__file__).resolve()
script_dir = script_path.parent

# If script is in docs folder, parent is repo root
if script_dir.name == 'docs':
    repo_root = script_dir.parent
else:
    repo_root = script_dir

# Add mofography package to path
mofography_path = repo_root
if str(mofography_path) not in sys.path:
    sys.path.insert(0, str(mofography_path))

import mofography as mgr
from ase.io import read

# Set up path to test structure (try multiple locations)
cif_name = 'P1-IRMOF-1.cif'
possible_cif_paths = [
    repo_root / 'examples' / 'TestCIFs' / cif_name,
    repo_root / 'mofography' / 'examples' / 'TestCIFs' / cif_name,
    repo_root / 'docs' / 'TestCIFs' / cif_name,
]

path_to_cif = None
for cif_path in possible_cif_paths:
    expanded_path = cif_path.expanduser()
    if expanded_path.exists():
        path_to_cif = expanded_path
        break

# Allow user to override with environment variable
if 'MOFOGRAPHY_TEST_CIF' in os.environ:
    path_to_cif = Path(os.environ['MOFOGRAPHY_TEST_CIF']).expanduser()

print(f"\nTest structure: {path_to_cif}")
if path_to_cif is None or not path_to_cif.exists():
    print(f"ERROR: CIF file not found")
    print(f"Tried paths:")
    for p in possible_cif_paths:
        print(f"  - {p}")
    print(f"\nYou can:")
    print(f"  1. Run from the repository root")
    print(f"  2. Set MOFOGRAPHY_TEST_CIF environment variable")
    sys.exit(1)

path_to_cif = str(path_to_cif)  # Convert to string for ASE

# =============================================================================
# STEP 2: Load Structure
# =============================================================================

print("\n" + "-" * 60)
print("STEP 2: Loading Structure")
print("-" * 60)

mat_atoms = read(path_to_cif)

print(f"Loaded structure: {mat_atoms.get_chemical_formula()}")
print(f"  Number of atoms: {len(mat_atoms)}")
print(f"  Unit cell: {mat_atoms.cell.cellpar()[:3]} Angstroms")

# =============================================================================
# STEP 3: Set Parameters
# =============================================================================

print("\n" + "-" * 60)
print("STEP 3: Setting Parameters")
print("-" * 60)

# Grid parameters
spacing = 0.5  # Grid spacing in Angstroms

# Segmentation parameters
mask_thickness = 1.0  # Void space threshold (Angstroms)
h = 0.2  # h-maxima persistence threshold
threshold_abs = 1.5  # Absolute threshold for maxima detection

# Periodic segmentation parameters
minimum_window_separation = 1.5  # Minimum distance between windows
wall_thickness = 0.08  # Boundary detection thickness (fraction of cell)
extra_padding = 10  # Extra padding in Angstroms

print(f"  spacing = {spacing} Å")
print(f"  mask_thickness = {mask_thickness} Å")
print(f"  h = {h}")
print(f"  threshold_abs = {threshold_abs} Å")
print(f"  wall_thickness = {wall_thickness}")
print(f"  extra_padding = {extra_padding} Å")

# =============================================================================
# STEP 4: Compute Distance Grid
# =============================================================================

print("\n" + "-" * 60)
print("STEP 4: Computing Distance Grid")
print("-" * 60)

# Compute distance grid using CPU backend
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=50000)

# Apply zoom and smoothing
dgrid = mgr.zoom_grid(dgrid, zoom_factor=1).compute()
dgrid = mgr.periodic_gaussian_filter(dgrid, sigma=spacing*2)

# CRITICAL: Convert to float32 for Rust backend
original_dtype = dgrid.dtype
dgrid = dgrid.astype(np.float32)

print(f"  Original dtype: {original_dtype}")
print(f"  Converted dtype: {dgrid.dtype}")
print(f"  Grid shape: {dgrid.shape}")
print(f"  Grid range: {dgrid.min():.2f} to {dgrid.max():.2f} Angstroms")
print(f"  Memory: {dgrid.nbytes / 1024**2:.1f} MB")

# =============================================================================
# STEP 5: Run Periodic Segmentation
# =============================================================================

print("\n" + "-" * 60)
print("STEP 5: Running Periodic Segmentation")
print("-" * 60)
print("  (This may take a minute...)")

import time
start_time = time.time()

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

elapsed_time = time.time() - start_time

print(f"  Completed in {elapsed_time:.2f} seconds")
print(f"\n  Results:")
print(f"    maxima_pbc:       {len(results['maxima_pbc'])} unique periodic maxima")
print(f"    unique_regions:   {len(results['unique_regions'])} unique regions")
print(f"    pbc_groups:       {len(results['pbc_groups'])} PBC groups")
print(f"    connections:      {len(results['connections'])} region connections")

# =============================================================================
# STEP 6: Extract Results
# =============================================================================

print("\n" + "-" * 60)
print("STEP 6: Extracting Results")
print("-" * 60)

maxima_pbc = results['maxima_pbc']
regions_pbc = results['regions_pbc']
regions = results['regions']
connections = results['connections']
all_maxima = results['all_maxima']
unique_regions = results['unique_regions']
maxima_info = results['maxima_info']
pbc_groups = results['pbc_groups']
close_maxima = results['close_or_inside_maxima']

# Maxima influencing the unit cell
maxima = all_maxima[:regions.max()]

print(f"  maxima_pbc:        {maxima_pbc.shape}")
print(f"  regions_pbc:       {regions_pbc.shape}, max={regions_pbc.max()}")
print(f"  regions:           {regions.shape}, max={regions.max()}")
print(f"  all_maxima:        {all_maxima.shape}")
print(f"  maxima (in unit cell): {maxima.shape}")

# =============================================================================
# STEP 7: Trim Connections (Recommended!)
# =============================================================================

print("\n" + "-" * 60)
print("STEP 7: Trimming Connections")
print("-" * 60)
print("  Note: Trimming to connections that influence the unit cell")
print("  This avoids fictitious pocket types from padded regions")

# Option 1: Use unique_regions (all regions that have voxels in unit cell)
# Option 2: Use close_or_inside_regions (more restrictive, regions in/near unit cell)
# Here we use unique_regions for broader coverage
connections_trimmed = np.array([
    c for c in connections
    if np.any(np.isin(c[0], unique_regions))
])

print(f"  Original connections: {len(connections)}")
print(f"  Trimmed connections: {len(connections_trimmed)}")

# =============================================================================
# STEP 8: Create Pore Graph (RAG) from Trimmed Connections
# =============================================================================

print("\n" + "-" * 60)
print("STEP 8: Creating Pore Graph (RAG)")
print("-" * 60)

# Create RAG from trimmed connections
# Use all_maxima since connections reference all periodic images
rag = mgr.rag_from_connections(connections_trimmed, all_maxima)

print(f"  Nodes (pores): {rag.number_of_nodes()}")
print(f"  Edges (windows): {rag.number_of_edges()}")

# =============================================================================
# STEP 9: Create Comprehensive Visualization
# =============================================================================

print("\n" + "-" * 60)
print("STEP 9: Creating Comprehensive Visualization")
print("-" * 60)

import plotly.graph_objects as go
import plotly.express as px

# Set up colorscale
cscale = px.colors.qualitative.Pastel

# Update structure with covalent bonds
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

# Create base traces
cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)

# Plot maxima (pore centers)
maxima_trace = mgr.plot_maxima(
    dgrid, maxima, mat_atoms,
    name='Pore centers',
    show_spheres=False,
    marker_colorscale=cscale
)

# Plot regions with PBC merged
surface_trace = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=1.5,
    region_labels=regions_pbc,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name='Pores (PBC merged)',
    stylize=True,
    opacity=0.8
)

# Combine into figure
fig = go.Figure(
    data=[cell_trace, surface_trace, maxima_trace, *frame_trace],
    layout=mgr.dark_layout
)

# Add trimmed connections
connections_traces = mgr.plot_connections_as_it_is(
    connections_trimmed,
    all_maxima,
    dgrid.shape,
    mat_atoms,
    cscale=cscale,
    plot_disconnected=False,
    show_spheres=False
)

for trace in connections_traces:
    fig.add_trace(trace)

fig.update_layout(
    scene_aspectmode='data',
    title=dict(
        text="IRMOF-1: Pore Network with Connections",
        x=0.5,
        xanchor='center',
        font=dict(size=16, color='white')
    ),
    scene=dict(
        xaxis=dict(title='X (Å)'),
        yaxis=dict(title='Y (Å)'),
        zaxis=dict(title='Z (Å)')
    ),
    legend=dict(
        orientation="v",
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01
    )
)

print(f"  Created figure with {len(fig.data)} traces:")
for i, trace in enumerate(fig.data):
    print(f"    {i+1}. {trace.name}")

# =============================================================================
# STEP 10: Save and Display Results
# =============================================================================

print("\n" + "-" * 60)
print("STEP 10: Saving Results")
print("-" * 60)

# Save as HTML
output_html = 'periodic_segmentation_quick_start_output.html'
fig.write_html(output_html)
print(f"  Saved interactive HTML to: {output_html}")

# Also save a static PNG version
output_png = 'periodic_segmentation_quick_start_output.png'
try:
    fig.write_image(output_png, width=1200, height=900, scale=1)
    print(f"  Saved static PNG to: {output_png}")
except Exception as e:
    print(f"  Could not save PNG (kaleido package required): {e}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"✓ Successfully analyzed {mat_atoms.get_chemical_formula()}")
print(f"✓ Found {len(maxima_pbc)} unique pores")
print(f"✓ Found {len(connections_trimmed)} pore connections (trimmed from {len(connections)})")
print(f"✓ Computation time: {elapsed_time:.2f} seconds")
print(f"✓ Output saved to: {output_html}")
print("\nVisualization includes:")
print(f"  • Pores (PBC merged) - correctly merged across boundaries")
print(f"  • Pore centers - local maxima in distance field")
print(f"  • Connections - windows between pores")
print("\nOpen the HTML file in your browser to explore the 3D visualization!")
print("\nNote: Pocket type detection is covered in the pockets_detection.ipynb tutorial")
print("=" * 60)

# Optional: Display in browser automatically
import webbrowser
import subprocess

try:
    # Try to open the HTML file
    if sys.platform == 'darwin':  # macOS
        subprocess.run(['open', output_html])
    elif sys.platform == 'win32':  # Windows
        os.startfile(output_html)
    else:  # Linux
        subprocess.run(['xdg-open', output_html])
except Exception as e:
    print(f"\nNote: Could not auto-open browser. Please open {output_html} manually.")

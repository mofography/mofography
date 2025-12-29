#!/usr/bin/env python3
"""
Distance Grid Quick Start - Standalone Test Script

This script demonstrates MOFography's distance grid computation methods.
Tests all three methods (CPU/Dask, GPU, Rust) and compares performance.

Run this script to test the workflow before using it in your own projects.
"""

import os
import sys
import numpy as np
import time
from pathlib import Path

# =============================================================================
# STEP 1: Setup and Imports
# =============================================================================

print("=" * 70)
print("MOFography - Distance Grid Quick Start")
print("=" * 70)

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

print("\n" + "-" * 70)
print("STEP 2: Loading Structure")
print("-" * 70)

mat_atoms = read(path_to_cif)
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

print(f"Loaded structure: {mat_atoms.get_chemical_formula()}")
print(f"  Number of atoms: {len(mat_atoms)}")
print(f"  Unit cell: {mat_atoms.cell.cellpar()[:3]} Angstroms")

# =============================================================================
# STEP 3: Set Parameters
# =============================================================================

print("\n" + "-" * 70)
print("STEP 3: Setting Parameters")
print("-" * 70)

spacing = 0.5  # Grid spacing in Angstroms

print(f"  spacing = {spacing} Å")

# =============================================================================
# STEP 4: Method 1 - CPU/Dask
# =============================================================================

print("\n" + "-" * 70)
print("STEP 4: Computing Distance Grid (CPU/Dask)")
print("-" * 70)

start_time = time.time()

dgrid_cpu = mgr.dgrid_from_atoms_cpu(
    mat_atoms,
    spacing=spacing,
    block_size=50000
)

# Convert from Dask array to NumPy
dgrid_cpu = dgrid_cpu.compute()

cpu_time = time.time() - start_time

print(f"  Completed in {cpu_time:.2f} seconds")
print(f"  Grid shape: {dgrid_cpu.shape}")
print(f"  Grid range: {dgrid_cpu.min():.2f} to {dgrid_cpu.max():.2f} Angstroms")
print(f"  Memory: {dgrid_cpu.nbytes / 1024**2:.1f} MB")
print(f"  Dtype: {dgrid_cpu.dtype}")

# =============================================================================
# STEP 5: Method 2 - GPU (if available)
# =============================================================================

print("\n" + "-" * 70)
print("STEP 5: Computing Distance Grid (GPU)")
print("-" * 70)

# Check if GPU is available
try:
    import cupy as cp
    has_gpu = True
    print("  GPU available!")
except ImportError:
    has_gpu = False
    print("  GPU not available (install cupy: pip install cupy-cuda1xx)")

dgrid_gpu = None
gpu_time = None

if has_gpu:
    start_time = time.time()

    try:
        dgrid_gpu = mgr.dgrid_from_atoms_gpu(
            mat_atoms,
            spacing=spacing,
            block_size=10000
        )

        print(f"  Device: {dgrid_gpu.device}")

        # Convert to NumPy for comparison
        dgrid_gpu_np = cp.asnumpy(dgrid_gpu)
        gpu_time = time.time() - start_time

        print(f"  Completed in {gpu_time:.2f} seconds")
        print(f"  Grid shape: {dgrid_gpu.shape}")
        print(f"  Grid range: {dgrid_gpu_np.min():.2f} to {dgrid_gpu_np.max():.2f} Angstroms")

    except Exception as e:
        print(f"  ERROR: {e}")
        print("  GPU computation failed, continuing...")
        has_gpu = False
        dgrid_gpu = None

# =============================================================================
# STEP 6: Method 3 - Rust (RECOMMENDED)
# =============================================================================

print("\n" + "-" * 70)
print("STEP 6: Computing Distance Grid (Rust - Recommended)")
print("-" * 70)

# Show thread count
import mgr_rust
print(f"  Using {mgr_rust.rayon_num_threads()} threads for parallel processing")

start_time = time.time()

dgrid_rust = mgr.dgrid_from_atoms_rust(
    mat_atoms,
    spacing=spacing
)

rust_time = time.time() - start_time

print(f"  Completed in {rust_time:.2f} seconds")
print(f"  Grid shape: {dgrid_rust.shape}")
print(f"  Grid range: {dgrid_rust.min():.2f} to {dgrid_rust.max():.2f} Angstroms")
print(f"  Memory: {dgrid_rust.nbytes / 1024**2:.1f} MB")
print(f"  Dtype: {dgrid_rust.dtype}")

# =============================================================================
# STEP 7: Performance Comparison
# =============================================================================

print("\n" + "-" * 70)
print("STEP 7: Performance Comparison")
print("-" * 70)

print(f"\n  Method          | Time (s) | Memory (MB)")
print(f"  " + "-" * 45)
print(f"  CPU/Dask        | {cpu_time:8.2f} | {dgrid_cpu.nbytes / 1024**2:10.1f}")
if has_gpu and gpu_time is not None:
    print(f"  GPU             | {gpu_time:8.2f} | (GPU VRAM)")
print(f"  Rust (Recommended)| {rust_time:8.2f} | {dgrid_rust.nbytes / 1024**2:10.1f}")

# =============================================================================
# STEP 8: Verify Results Match
# =============================================================================

print("\n" + "-" * 70)
print("STEP 8: Verifying Results")
print("-" * 70)

# Compare CPU and Rust
diff_cpu_rust = np.abs(dgrid_cpu - dgrid_rust)
print(f"\n  CPU vs Rust:")
print(f"    Max difference: {diff_cpu_rust.max():.6f} Å")
print(f"    Mean difference: {diff_cpu_rust.mean():.6f} Å")
print(f"    Results match: {np.allclose(dgrid_cpu, dgrid_rust, rtol=1e-5)}")

if has_gpu and dgrid_gpu is not None:
    dgrid_gpu_np = cp.asnumpy(dgrid_gpu)
    diff_cpu_gpu = np.abs(dgrid_cpu - dgrid_gpu_np)
    print(f"\n  CPU vs GPU:")
    print(f"    Max difference: {diff_cpu_gpu.max():.6f} Å")
    print(f"    Mean difference: {diff_cpu_gpu.mean():.6f} Å")
    print(f"    Results match: {np.allclose(dgrid_cpu, dgrid_gpu_np, rtol=1e-4)}")

# =============================================================================
# STEP 9: Create Visualization
# =============================================================================

print("\n" + "-" * 70)
print("STEP 9: Creating Visualization")
print("-" * 70)

import plotly.graph_objects as go
import plotly.express as px

# Use Rust result (fastest)
dgrid = dgrid_rust

# Set up colorscale
cscale = px.colors.qualitative.Pastel

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

# Create figure
cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)

fig = go.Figure(
    data=[cell_trace, isosurface_trace, *frame_trace],
    layout=mgr.dark_layout
)

fig.update_layout(
    scene_aspectmode='data',
    title=dict(
        text=f"Distance Grid Isosurface ({mat_atoms.get_chemical_formula()})",
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

print("\n" + "-" * 70)
print("STEP 10: Saving Results")
print("-" * 70)

# Save as HTML
output_html = 'distance_grid_quick_start_output.html'
fig.write_html(output_html)
print(f"  Saved interactive HTML to: {output_html}")

# Also save a static PNG version
output_png = 'distance_grid_quick_start_output.png'
try:
    fig.write_image(output_png, width=1200, height=900, scale=1)
    print(f"  Saved static PNG to: {output_png}")
except Exception as e:
    print(f"  Could not save PNG (kaleido package required): {e}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"✓ Successfully analyzed {mat_atoms.get_chemical_formula()}")
print(f"✓ Distance grid computed at {spacing} Å spacing")
print(f"✓ Grid size: {dgrid.shape}")
print(f"✓ Performance:")
print(f"    - CPU/Dask: {cpu_time:.2f} s")
if has_gpu and gpu_time is not None:
    print(f"    - GPU: {gpu_time:.2f} s")
print(f"    - Rust: {rust_time:.2f} s (recommended)")
print(f"\nRecommendation:")
print(f"  → Use dgrid_from_atoms_rust() for most cases")
print(f"  → Use dgrid_from_atoms_cpu() for very large systems")
if has_gpu:
    print(f"  → Use dgrid_from_atoms_gpu() if you have a good GPU")
print(f"\nOutput saved to: {output_html}")
print("\nOpen the HTML file in your browser to explore the 3D visualization!")
print("=" * 70)

# Optional: Display in browser automatically
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

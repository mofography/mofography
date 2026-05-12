#!/usr/bin/env python3
"""
Region Segmentation Tutorial - Standalone Script

This script demonstrates the region segmentation workflow in MOFography.
Run this from the docs directory:
    python region_segmentation_quick_start.py
"""

import os
import sys
import numpy as np
from pathlib import Path

# Setup paths
script_path = Path(__file__).resolve()
if script_path.parent.name == "docs":
    repo_root = script_path.parent.parent
else:
    repo_root = script_path.parent

sys.path.insert(0, str(repo_root))

import mofography as mgr
from ase.io import read
import plotly.graph_objects as go
import plotly.express as px

# =============================================================================
# Configuration
# =============================================================================

cif_name = "P1-IRMOF-1.cif"
possible_cif_paths = [
    repo_root / "examples" / "TestCIFs" / cif_name,
    repo_root / "mofography" / "examples" / "TestCIFs" / cif_name,
]

path_to_cif = None
for cif_path in possible_cif_paths:
    expanded_path = cif_path.expanduser()
    if expanded_path.exists():
        path_to_cif = expanded_path
        break

if path_to_cif is None:
    print(f"Error: Could not find {cif_name}")
    print(f"Searched in:")
    for p in possible_cif_paths:
        print(f"  {p}")
    sys.exit(1)

# Parameters
spacing = 0.5
mask_thickness = 1.0
h = 0.2
threshold_abs = 1.5
minimum_window_separation = 1.5
wall_thickness = 0.08
extra_padding = 10
isosurface_value = 1.5

# =============================================================================
# Step 1: Load Structure
# =============================================================================

print("\n" + "=" * 60)
print("STEP 1: Loading Structure")
print("=" * 60)

mat_atoms = read(str(path_to_cif))
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

print(f"Loaded: {mat_atoms.get_chemical_formula()}")
print(f"  Atoms: {len(mat_atoms)}")
print(f"  Cell: {mat_atoms.cell.cellpar()[:3]} Å")

# =============================================================================
# Step 2: Compute Distance Grid
# =============================================================================

print("\n" + "=" * 60)
print("STEP 2: Computing Distance Grid")
print("=" * 60)

print(f"Computing distance grid (spacing={spacing} Å)...")
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=50000)
dgrid = mgr.zoom_grid(dgrid, zoom_factor=1).compute()
dgrid = mgr.periodic_gaussian_filter(dgrid, sigma=spacing * 2)

# CRITICAL: Convert to float32 for Rust backend
original_dtype = dgrid.dtype
dgrid = dgrid.astype(np.float32)

print(f"  Grid shape: {dgrid.shape}")
print(f"  Original dtype: {original_dtype}")
print(f"  Converted dtype: {dgrid.dtype}")
print(f"  Distance range: {dgrid.min():.2f} to {dgrid.max():.2f} Å")

# =============================================================================
# Step 3: Run Periodic Segmentation
# =============================================================================

print("\n" + "=" * 60)
print("STEP 3: Running Periodic Segmentation")
print("=" * 60)

print("Running do_periodic_segmentation_rust()...")
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
    compute_connections=True,
)

print("\n✅ Segmentation complete!")
print(f"  maxima_pbc (unique pores): {len(results['maxima_pbc'])}")
print(f"  unique_regions: {len(results['unique_regions'])}")
print(f"  pbc_groups: {len(results['pbc_groups'])}")
print(f"  connections: {len(results['connections'])}")

# =============================================================================
# Step 4: Extract Results
# =============================================================================

print("\n" + "=" * 60)
print("STEP 4: Extracting Results")
print("=" * 60)

# Extract key results
maxima_pbc = results["maxima_pbc"]
regions_pbc = results["regions_pbc"]
regions = results["regions"]
all_maxima = results["all_maxima"]
unique_regions = results["unique_regions"]
close_maxima = results["close_or_inside_maxima"]

# Get maxima influencing the unit cell (clever indexing trick!)
maxima = all_maxima[: regions.max()]

print(f"\nRegion Arrays:")
print(f"  regions_pbc (merged): shape={regions_pbc.shape}, max={regions_pbc.max()}")
print(f"  regions (split): shape={regions.shape}, max={regions.max()}")
print(f"\nNumber of unique pores: {regions_pbc.max()}")

print(f"\nMaxima Arrays:")
print(f"  maxima_pbc: {len(maxima_pbc)} (one per unique pore)")
print(f"  all_maxima: {len(all_maxima)} (including periodic images)")
print(f"  maxima (in unit cell): {len(maxima)}")

# =============================================================================
# Step 5: Trim Connections (Optional)
# =============================================================================

print("\n" + "=" * 60)
print("STEP 5: Trimming Connections")
print("=" * 60)

connections = results["connections"]
connections_trimmed = np.array(
    [c for c in connections if np.any(np.isin(c[0], unique_regions))]
)

print(f"Original connections: {len(connections)}")
print(f"Trimmed connections: {len(connections_trimmed)}")

# =============================================================================
# Step 6: Create RAG (Optional)
# =============================================================================

print("\n" + "=" * 60)
print("STEP 6: Creating Pore Graph (RAG)")
print("=" * 60)

rag = mgr.rag_from_connections(connections_trimmed, all_maxima)
print(f"RAG nodes (pores): {rag.number_of_nodes()}")
print(f"RAG edges (connections): {rag.number_of_edges()}")

# =============================================================================
# Step 7: Visualization
# =============================================================================

print("\n" + "=" * 60)
print("STEP 7: Creating Visualization")
print("=" * 60)

print("Generating visualization...")

# Set up colorscale
cscale = px.colors.qualitative.Pastel

# Create traces
cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)

# Plot merged regions (PBC applied)
surface_trace = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=isosurface_value,
    region_labels=regions_pbc,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name="Pores (PBC merged)",
    stylize=True,
    opacity=0.8,
)

# Plot pore centers
maxima_trace = mgr.plot_maxima(
    dgrid,
    maxima,
    mat_atoms,
    name="Pore centers",
    show_spheres=False,
    marker_colorscale=cscale,
)

# Plot connections
connections_traces = mgr.plot_connections_as_it_is(
    connections_trimmed,
    all_maxima,
    dgrid.shape,
    mat_atoms,
    cscale=cscale,
    plot_disconnected=False,
    show_spheres=False,
)

# Combine into figure
fig = go.Figure(
    data=[cell_trace, surface_trace, maxima_trace, *connections_traces, *frame_trace],
    layout=mgr.dark_layout,
)

fig.update_layout(
    scene_aspectmode="data",
    title=dict(
        text=f"IRMOF-1: Region Segmentation ({regions_pbc.max()} unique pores)",
        x=0.5,
        xanchor="center",
        font=dict(size=16, color="white"),
    ),
    scene=dict(
        xaxis=dict(title="X (Å)"), yaxis=dict(title="Y (Å)"), zaxis=dict(title="Z (Å)")
    ),
)

# Save outputs
output_file = "region_segmentation_results.html"
fig.write_html(output_file)
print(f"\n✅ Visualization saved to: {output_file}")

# Also save as PNG for documentation
try:
    png_file = "region_segmentation_results.png"
    fig.write_image(png_file, width=1200, height=900, scale=2)
    print(f"✅ PNG saved to: {png_file}")
except Exception as e:
    print(f"⚠️  Could not save PNG (kaleido may not be installed): {e}")

# =============================================================================
# Summary
# =============================================================================

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)

print(f"""
Results:
  - Structure: {mat_atoms.get_chemical_formula()}
  - Grid: {dgrid.shape} voxels at {spacing} Å spacing
  - Unique pores found: {regions_pbc.max()}
  - PBC groups: {len(results["pbc_groups"])}
  - Connections (trimmed): {len(connections_trimmed)}
  - RAG nodes: {rag.number_of_nodes()}
  - RAG edges: {rag.number_of_edges()}

Output files:
  - {output_file}
""")

print("\n🎉 Region segmentation complete!")
print("\nNext steps:")
print("  1. Open the HTML file in a browser to explore the visualization")
print("  2. Try adjusting parameters (h, mask_thickness, etc.)")
print("  3. Check out the region_segmentation.md tutorial for more details")

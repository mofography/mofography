#!/usr/bin/env python3
"""
Generate images for the Region Segmentation tutorial.

This script creates PNG visualizations using the region segmentation workflow.
Images are saved with white backgrounds for documentation use.
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
from plotly.subplots import make_subplots

# =============================================================================
# Load Structure
# =============================================================================

print("Loading structure...")
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
    sys.exit(1)

mat_atoms = read(str(path_to_cif))
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

print(f"Loaded: {mat_atoms.get_chemical_formula()}")

# =============================================================================
# Parameters
# =============================================================================

print("\nSetting parameters...")
spacing = 0.5
mask_thickness = 1.0
h = 0.2
threshold_abs = 1.5
minimum_window_separation = 1.5
wall_thickness = 0.08
extra_padding = 10

# =============================================================================
# Compute Distance Grid
# =============================================================================

print("\nComputing distance grid...")
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=50000)
dgrid = mgr.zoom_grid(dgrid, zoom_factor=1).compute()
dgrid = mgr.periodic_gaussian_filter(dgrid, sigma=spacing * 2)
dgrid = dgrid.astype(np.float32)

print(f"  Grid shape: {dgrid.shape}")

# =============================================================================
# Run Periodic Segmentation
# =============================================================================

print("\nRunning periodic segmentation...")
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

# Extract results
maxima_pbc = results["maxima_pbc"]
regions_pbc = results["regions_pbc"]
regions = results["regions"]
all_maxima = results["all_maxima"]
connections = results["connections"]
unique_regions = results["unique_regions"]
close_maxima = results["close_or_inside_maxima"]
pbc_groups = results["pbc_groups"]

maxima = all_maxima[: regions.max()]

print(f"  Found {len(maxima_pbc)} unique periodic maxima")
print(f"  Found {len(unique_regions)} unique regions")

# =============================================================================
# Create output directory
# =============================================================================

output_dir = script_path.parent / "tutorial_images"
output_dir.mkdir(exist_ok=True)
print(f"\nSaving images to: {output_dir}")

# =============================================================================
# Generate Image 1: Regions with PBC (merged)
# =============================================================================

print("\nGenerating image: regions_pbc_merged.png")

cscale = px.colors.qualitative.Pastel
surface_value = 1.5

surface_trace_pbc = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=surface_value,
    region_labels=regions_pbc,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name="Regions (PBC merged)",
    stylize=True,
    opacity=0.8,
)

maxima_trace_pbc = mgr.plot_maxima(
    dgrid,
    maxima_pbc,
    mat_atoms,
    name="Pore centers (PBC)",
    show_spheres=False,
    marker_colorscale=cscale,
)

cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)

fig1 = go.Figure(
    data=[cell_trace, surface_trace_pbc, maxima_trace_pbc, *frame_trace],
    layout=mgr.white_layout,
)

fig1.update_layout(
    scene_aspectmode="data",
    title=dict(
        text="IRMOF-1: Pore Regions with Periodic Boundaries Merged",
        x=0.5,
        xanchor="center",
        font=dict(size=14, color="black"),
    ),
    scene=dict(
        xaxis=dict(title="X (Å)"), yaxis=dict(title="Y (Å)"), zaxis=dict(title="Z (Å)")
    ),
)

fig1.write_image(output_dir / "regions_pbc_merged.png", width=800, height=600, scale=2)
print("  Saved: regions_pbc_merged.png")

# =============================================================================
# Generate Image 2: Regions without PBC (split at boundaries)
# =============================================================================

print("\nGenerating image: regions_split_no_pbc.png")

surface_trace = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=surface_value,
    region_labels=regions,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name="Regions (no PBC)",
    stylize=True,
    opacity=0.8,
)

maxima_trace = mgr.plot_maxima(
    dgrid,
    maxima,
    mat_atoms,
    name="Pore centers (split)",
    show_spheres=False,
    marker_colorscale=cscale,
)

fig2 = go.Figure(
    data=[cell_trace, surface_trace, maxima_trace, *frame_trace],
    layout=mgr.white_layout,
)

fig2.update_layout(
    scene_aspectmode="data",
    title=dict(
        text="IRMOF-1: Pore Regions Split at Unit Cell Boundaries",
        x=0.5,
        xanchor="center",
        font=dict(size=14, color="black"),
    ),
    scene=dict(
        xaxis=dict(title="X (Å)"), yaxis=dict(title="Y (Å)"), zaxis=dict(title="Z (Å)")
    ),
)

fig2.write_image(
    output_dir / "regions_split_no_pbc.png", width=800, height=600, scale=2
)
print("  Saved: regions_split_no_pbc.png")

# =============================================================================
# Generate Image 3: All maxima types (close, unique, pbc)
# =============================================================================

print("\nGenerating image: maxima_types_comparison.png")

# Combine surface trace with all maxima types
close_maxima_trace = mgr.plot_maxima(
    dgrid,
    close_maxima,
    mat_atoms,
    name="Close/inside maxima",
    show_spheres=False,
    marker_colorscale=cscale,
)

unique_maxima_trace = mgr.plot_maxima(
    dgrid,
    all_maxima[: len(unique_regions)],
    mat_atoms,
    name="Unique region maxima",
    show_spheres=False,
    marker_colorscale=cscale,
)

pbc_maxima_trace = mgr.plot_maxima(
    dgrid,
    maxima_pbc,
    mat_atoms,
    name="PBC representatives",
    show_spheres=False,
    marker_colorscale=cscale,
)

fig3 = go.Figure(
    data=[
        cell_trace,
        surface_trace_pbc,
        close_maxima_trace,
        unique_maxima_trace,
        pbc_maxima_trace,
        *frame_trace,
    ],
    layout=mgr.white_layout,
)

fig3.update_layout(
    scene_aspectmode="data",
    title=dict(
        text="IRMOF-1: Types of Maxima Detected",
        x=0.5,
        xanchor="center",
        font=dict(size=14, color="black"),
    ),
    scene=dict(
        xaxis=dict(title="X (Å)"), yaxis=dict(title="Y (Å)"), zaxis=dict(title="Z (Å)")
    ),
)

fig3.write_image(
    output_dir / "maxima_types_comparison.png", width=800, height=600, scale=2
)
print("  Saved: maxima_types_comparison.png")

# =============================================================================
# Generate Image 4: Cross-section showing segmentation
# =============================================================================

print("\nGenerating image: segmentation_cross_section.png")

# Take a middle slice to show segmentation
z_slice = dgrid.shape[2] // 2

import plotly.express as px
import plotly.graph_objects as go

fig4 = go.Figure()

# Add cross-section image
fig4.add_trace(
    go.Heatmap(
        z=dgrid[:, :, z_slice].T,
        colorscale="Viridis",
        showscale=True,
        colorbar=dict(title="Distance (Å)"),
    )
)

# Add maxima on cross-section
maxima_z = maxima[:, 2]
maxima_on_slice = maxima_z > z_slice * spacing
maxima_on_slice = maxima_on_slice & (maxima_z < (z_slice + 1) * spacing)

if np.any(maxima_on_slice):
    fig4.add_trace(
        go.Scatter(
            x=maxima[maxima_on_slice, 0],
            y=maxima[maxima_on_slice, 1],
            mode="markers",
            marker=dict(size=15, color="red", line=dict(width=2, color="white")),
            name="Maxima on slice",
        )
    )

fig4.update_layout(
    title=dict(
        text=f"Distance Grid Cross-Section at z={z_slice * spacing:.1f} Å",
        x=0.5,
        xanchor="center",
        font=dict(size=14, color="black"),
    ),
    xaxis=dict(title="X (grid indices)", scaleanchor="y"),
    yaxis=dict(title="Y (grid indices)"),
    plot_bgcolor="white",
)

fig4.write_image(
    output_dir / "segmentation_cross_section.png", width=800, height=800, scale=2
)
print("  Saved: segmentation_cross_section.png")

# =============================================================================
# Generate Image 5: Parameter effect - different h values
# =============================================================================

print("\nGenerating image: parameter_effect_h_values.png")

h_values = [0.1, 0.2, 0.5]
fig5 = make_subplots(
    rows=1,
    cols=3,
    subplot_titles=[f"h={h} Å" for h in h_values],
    specs=[[{"type": "surface"}, {"type": "surface"}, {"type": "surface"}]],
)

for i, h_val in enumerate(h_values, 1):
    # Re-run segmentation with different h value
    results_h = mgr.do_periodic_segmentation_rust(
        dgrid,
        mat_atoms,
        spacing=spacing,
        mask_thickness=mask_thickness,
        h=h_val,
        threshold_abs=threshold_abs,
        minimum_window_separation=minimum_window_separation,
        wall_thickness=wall_thickness,
        extra_padding=extra_padding,
        return_padded_regions=True,
        compute_connections=False,
    )

    regions_h = results_h["regions_pbc"]

    # Get simple surface visualization
    surface_trace = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=surface_value,
        region_labels=regions_h,
        ase_atoms=mat_atoms,
        cscale=cscale,
        opacity=0.8,
    )

    # Add to subplot (using only one surface trace)
    fig5.add_trace(surface_trace, row=1, col=i)

fig5.update_layout(
    title=dict(
        text="Effect of h Parameter on Maxima Detection",
        x=0.5,
        xanchor="center",
        font=dict(size=14, color="black"),
    ),
    height=400,
    showlegend=False,
)

for i in range(1, 4):
    fig5.update_scenes(
        xaxis=dict(title="X (Å)"),
        yaxis=dict(title="Y (Å)"),
        zaxis=dict(title="Z (Å)"),
        row=1,
        col=i,
    )

fig5.write_image(
    output_dir / "parameter_effect_h_values.png", width=1200, height=400, scale=2
)
print("  Saved: parameter_effect_h_values.png")

print("\n✅ All images generated successfully!")
print(f"Images saved to: {output_dir}")

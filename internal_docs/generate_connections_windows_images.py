#!/usr/bin/env python3
"""
Generate 3D visualization images for the Connections and Windows tutorial.

This script creates PNG visualizations from real IRMOF-1 data:
1. Main pore network with all windows
2. Window size histogram
"""

import os
import sys
import numpy as np
from pathlib import Path

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

print("=" * 60)
print("CONNECTIONS AND WINDOWS - IMAGE GENERATION")
print("=" * 60)

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

print(f"\nLoading structure: {cif_name}")
mat_atoms = read(str(path_to_cif))
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
print(f"  Loaded: {mat_atoms.get_chemical_formula()}")

output_dir = script_path.parent / "tutorial_images"
output_dir.mkdir(exist_ok=True)
print(f"\nOutput directory: {output_dir}")

spacing = 0.5
print(f"\nComputing distance grid (spacing={spacing})...")
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=50000)
dgrid = mgr.zoom_grid(dgrid, zoom_factor=1).compute()
dgrid = mgr.periodic_gaussian_filter(dgrid, sigma=spacing * 2)
dgrid = dgrid.astype(np.float32)
print(f"  Grid shape: {dgrid.shape}")

params = {
    "mask_thickness": 1.0,
    "h": 0.2,
    "threshold_abs": 1.5,
    "minimum_window_separation": 1.5,
    "wall_thickness": 0.08,
    "extra_padding": 10,
}

print(f"\nRunning periodic segmentation...")
results = mgr.do_periodic_segmentation_rust(
    dgrid, mat_atoms, spacing=spacing, **params, compute_connections=True
)

maxima_pbc = results["maxima_pbc"]
regions_pbc = results["regions_pbc"]
connections = results["connections"]
all_maxima = results["all_maxima"]
unique_regions = results["unique_regions"]

connections_trimmed = np.array(
    [c for c in connections if np.any(np.isin(c[0], unique_regions))]
)

print(f"  Unique pores: {len(maxima_pbc)}")
print(f"  Total windows: {len(connections_trimmed)}")

cscale = px.colors.qualitative.Pastel

print("\n" + "-" * 60)
print("Generating Image 1: windows_3d.png (Full pore network)")
print("-" * 60)

cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)

surface_trace = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=1.5,
    region_labels=regions_pbc,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name="Pores",
    stylize=True,
    opacity=0.6,
)

connections_traces = mgr.plot_connections_as_it_is(
    connections_trimmed,
    all_maxima,
    dgrid.shape,
    mat_atoms,
    cscale=cscale,
    plot_disconnected=False,
    show_spheres=True,
)

maxima_trace = mgr.plot_maxima(
    dgrid,
    maxima_pbc,
    mat_atoms,
    name="Pore centers",
    show_spheres=False,
    marker_colorscale=cscale,
)

fig1 = go.Figure(
    data=[
        cell_trace,
        surface_trace,
        maxima_trace,
        *connections_traces,
        *frame_trace,
    ],
    layout=mgr.white_layout,
)

fig1.update_layout(
    scene_aspectmode="data",
    title=dict(
        text=f"IRMOF-1: Pore Network with Windows ({len(maxima_pbc)} pores, {len(connections_trimmed)} windows)",
        x=0.5,
        xanchor="center",
        font=dict(size=14, color="black"),
    ),
    scene=dict(
        xaxis=dict(title="X (A)"),
        yaxis=dict(title="Y (A)"),
        zaxis=dict(title="Z (A)"),
    ),
)

fig1.write_image(output_dir / "windows_3d.png", width=900, height=700, scale=2)
print(f"  Saved: windows_3d.png")

print("\n" + "-" * 60)
print("Generating Image 2: windows_histogram.png (Window size distribution)")
print("-" * 60)

all_radii = np.array([c[4][0] for c in connections_trimmed])
all_diameters = 2 * all_radii

fig2 = px.histogram(
    x=all_diameters,
    nbins=15,
    labels={"x": "Window Diameter (A)"},
    title="IRMOF-1: Distribution of Window Sizes",
    color_discrete_sequence=["steelblue"],
)

fig2.add_vline(x=all_diameters.min(), line_dash="dash", line_color="red", line_width=2)
fig2.add_annotation(
    x=all_diameters.min(),
    y=max(np.histogram(all_diameters, bins=15)[0]) * 0.8,
    text=f"PLD = {all_diameters.min():.2f} A",
    showarrow=False,
    font=dict(color="red", size=12),
)

fig2.update_layout(
    showlegend=False,
    height=350,
    width=700,
    title=dict(x=0.5, font=dict(size=16)),
    xaxis=dict(title_font=dict(size=14), tickfont=dict(size=12), range=[4, 10]),
    yaxis=dict(title="Count", title_font=dict(size=14), tickfont=dict(size=12)),
    plot_bgcolor="white",
)

fig2.update_xaxes(showline=True, linecolor="black", mirror=True)
fig2.update_yaxes(showline=True, linecolor="black", mirror=True)

fig2.write_image(output_dir / "windows_histogram.png", width=700, height=350, scale=2)
print(f"  Saved: windows_histogram.png")

print(f"\nWindow diameter statistics:")
print(f"  PLD (min): {all_diameters.min():.3f} A")
print(f"  Max: {all_diameters.max():.3f} A")
print(f"  Mean: {all_diameters.mean():.3f} A")
print(f"  Median: {np.median(all_diameters):.3f} A")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Images generated in: {output_dir}")
print(
    f"  1. windows_3d.png      - Full pore network ({len(connections_trimmed)} windows)"
)
print(f"  2. windows_histogram.png - Window size distribution")
print(f"\nStatistics:")
print(f"  Structure: {mat_atoms.get_chemical_formula()}")
print(f"  Unique pores: {len(maxima_pbc)}")
print(f"  PLD: {all_diameters.min():.3f} A")
print("=" * 60)
print("\nAll images generated successfully!")

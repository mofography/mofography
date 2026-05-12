#!/usr/bin/env python3
"""
Generate images for the Pore Graph tutorial.

This script creates PNG visualizations including:
1. Conceptual 2D diagram of pore space vs pore graph
2. Actual 3D pore network visualization
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
connections = results["connections"]
all_maxima = results["all_maxima"]
unique_regions = results["unique_regions"]

# Maxima for unit cell
maxima = all_maxima[: regions.max()]

print(f"  Found {len(maxima_pbc)} unique periodic maxima")
print(f"  Found {len(connections)} total connections")

# =============================================================================
# Create output directory
# =============================================================================

output_dir = script_path.parent / "tutorial_images"
output_dir.mkdir(exist_ok=True)
print(f"\nSaving images to: {output_dir}")

# =============================================================================
# Generate Image 1: Conceptual 2D diagram
# =============================================================================

print("\nGenerating image: pore_graph_concept.png")

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Circle, Rectangle

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Pore space (physical representation)
ax1 = axes[0]
ax1.set_xlim(-0.5, 3.5)
ax1.set_ylim(-0.5, 3.5)
ax1.set_aspect("equal")
ax1.set_title(
    "Physical Pore Space\n(Segmented Regions)", fontsize=12, fontweight="bold"
)

# Draw unit cell (thick black lines)
rect = patches.Rectangle(
    (0, 0), 1, 1, linewidth=3, edgecolor="black", facecolor="lightblue", alpha=0.3
)
ax1.add_patch(rect)
ax1.text(0.5, 1.1, "Unit Cell", ha="center", fontsize=10)

# Draw periodic images (dotted lines)
for i in [-1, 0, 1, 2]:
    for j in [-1, 0, 1, 2]:
        if i == 0 and j == 0:
            continue
        rect_dotted = patches.Rectangle(
            (i, j),
            1,
            1,
            linewidth=1,
            edgecolor="gray",
            facecolor="none",
            linestyle="--",
            alpha=0.5,
        )
        ax1.add_patch(rect_dotted)

# Draw pore compartments (simplified - 4 pores in a 2x2 arrangement)
pore_positions = [(0.25, 0.25), (0.75, 0.25), (0.25, 0.75), (0.75, 0.75)]
pore_labels = ["1 (0,0,0)", "2 (0,0,0)", "3 (0,0,0)", "4 (0,0,0)"]
colors = ["red", "blue", "green", "orange"]

for i, (pos, label, color) in enumerate(zip(pore_positions, pore_labels, colors)):
    circle = patches.Circle(pos, 0.2, color=color, alpha=0.6, ec="black", linewidth=2)
    ax1.add_patch(circle)
    ax1.text(
        pos[0],
        pos[1],
        str(i + 1),
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        color="white",
    )
    ax1.text(pos[0], pos[1] - 0.35, label, ha="center", va="top", fontsize=8)

# Add periodic images
for i, (pos, color) in enumerate(zip(pore_positions, colors)):
    # Right periodic image
    ax1.add_patch(
        patches.Circle(
            (pos[0] + 1, pos[1]),
            0.2,
            color=color,
            alpha=0.3,
            ec="black",
            linestyle="--",
            linewidth=1,
            fill=False,
        )
    )
    ax1.text(
        pos[0] + 1,
        pos[1],
        str(i + 1),
        ha="center",
        va="center",
        fontsize=8,
        color=color,
    )

    # Top periodic image
    ax1.add_patch(
        patches.Circle(
            (pos[0], pos[1] + 1),
            0.2,
            color=color,
            alpha=0.3,
            ec="black",
            linestyle="--",
            linewidth=1,
            fill=False,
        )
    )
    ax1.text(
        pos[0],
        pos[1] + 1,
        str(i + 1),
        ha="center",
        va="center",
        fontsize=8,
        color=color,
    )

    # Diagonal periodic image
    ax1.add_patch(
        patches.Circle(
            (pos[0] + 1, pos[1] + 1),
            0.2,
            color=color,
            alpha=0.3,
            ec="black",
            linestyle="--",
            linewidth=1,
            fill=False,
        )
    )
    ax1.text(
        pos[0] + 1,
        pos[1] + 1,
        str(i + 1),
        ha="center",
        va="center",
        fontsize=8,
        color=color,
    )

ax1.set_xlabel("X (unit cell coordinates)")
ax1.set_ylabel("Y (unit cell coordinates)")

# Right: Pore graph
ax2 = axes[1]
ax2.set_xlim(-0.5, 3.5)
ax2.set_ylim(-0.5, 3.5)
ax2.set_aspect("equal")
ax2.set_title("Pore Graph\n(Network Representation)", fontsize=12, fontweight="bold")

# Draw unit cell outline (dashed)
rect2 = patches.Rectangle(
    (0, 0), 1, 1, linewidth=2, edgecolor="black", facecolor="none", linestyle="--"
)
ax2.add_patch(rect2)

# Draw periodic cell outlines
for i in [-1, 0, 1, 2]:
    for j in [-1, 0, 1, 2]:
        if i == 0 and j == 0:
            continue
        rect_dotted = patches.Rectangle(
            (i, j),
            1,
            1,
            linewidth=1,
            edgecolor="gray",
            facecolor="none",
            linestyle=":",
            alpha=0.3,
        )
        ax2.add_patch(rect_dotted)

# Draw nodes (pore centers)
for i, (pos, label, color) in enumerate(zip(pore_positions, pore_labels, colors)):
    # Main node
    ax2.scatter(
        pos[0], pos[1], s=200, c=color, edgecolors="black", linewidth=2, zorder=5
    )
    ax2.text(
        pos[0],
        pos[1],
        str(i + 1),
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )

    # Periodic nodes (empty circles)
    for dx, dy in [(1, 0), (0, 1), (1, 1)]:
        ax2.scatter(
            pos[0] + dx,
            pos[1] + dy,
            s=100,
            c="none",
            edgecolors=color,
            linewidth=1.5,
            linestyle="--",
            zorder=4,
        )
        ax2.text(
            pos[0] + dx,
            pos[1] + dy,
            str(i + 1),
            ha="center",
            va="center",
            fontsize=7,
            color=color,
        )

# Draw edges (connections between pores in the unit cell)
edge_pairs = [(0, 1), (0, 2), (1, 3), (2, 3)]  # Adjacent pores
for i1, i2 in edge_pairs:
    p1 = pore_positions[i1]
    p2 = pore_positions[i2]
    ax2.plot([p1[0], p2[0]], [p1[1], p2[1]], "k-", linewidth=2, zorder=3)

# Draw periodic edges (dotted)
for i1, i2 in edge_pairs:
    p1 = pore_positions[i1]
    p2 = pore_positions[i2]
    # Horizontal periodic
    ax2.plot(
        [p1[0] + 1, p2[0] + 1], [p1[1], p2[1]], "k--", linewidth=1, alpha=0.5, zorder=2
    )
    # Vertical periodic
    ax2.plot(
        [p1[0], p2[0]], [p1[1] + 1, p2[1] + 1], "k--", linewidth=1, alpha=0.5, zorder=2
    )

# Add legend
ax2.scatter([], [], s=200, c="gray", edgecolors="black", label="Node (physical)")
ax2.scatter(
    [],
    [],
    s=100,
    c="none",
    edgecolors="gray",
    linestyle="--",
    label="Node (periodic image)",
)
ax2.plot([], [], "k-", linewidth=2, label="Edge (internal)")
ax2.plot([], [], "k--", linewidth=1, alpha=0.5, label="Edge (periodic)")
ax2.legend(loc="upper right", fontsize=8)

ax2.set_xlabel("X (unit cell coordinates)")
ax2.set_ylabel("Y (unit cell coordinates)")

plt.tight_layout()
fig.savefig(
    output_dir / "pore_graph_concept.png",
    dpi=150,
    bbox_inches="tight",
    facecolor="white",
)
plt.close()

print("  Saved: pore_graph_concept.png")

# =============================================================================
# Generate Image 2: Actual 3D pore network
# =============================================================================

print("\nGenerating image: pore_network_3d.png")

cscale = px.colors.qualitative.Pastel
surface_value = 1.5

# Create base traces
cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)

surface_trace = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=surface_value,
    region_labels=regions_pbc,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name="Pores",
    stylize=True,
    opacity=0.6,
)

# Trim connections
connections_trimmed = np.array(
    [c for c in connections if np.any(np.isin(c[0], unique_regions))]
)

# Use plot_connections_as_it_is for visualization
connections_traces = mgr.plot_connections_as_it_is(
    connections_trimmed,
    all_maxima,
    dgrid.shape,
    mat_atoms,
    cscale=cscale,
    plot_disconnected=False,
    show_spheres=False,
)

# Plot pore centers
maxima_trace = mgr.plot_maxima(
    dgrid,
    maxima_pbc,
    mat_atoms,
    name="Pore centers",
    show_spheres=False,
    marker_colorscale=cscale,
)

# Create 3D figure
fig2 = go.Figure(
    data=[cell_trace, surface_trace, maxima_trace, *connections_traces, *frame_trace],
    layout=mgr.white_layout,
)

fig2.update_layout(
    scene_aspectmode="data",
    title=dict(
        text=f"IRMOF-1: Pore Network ({len(maxima_pbc)} unique pores, {len(connections_trimmed)} connections)",
        x=0.5,
        xanchor="center",
        font=dict(size=14, color="black"),
    ),
    scene=dict(
        xaxis=dict(title="X (Å)"), yaxis=dict(title="Y (Å)"), zaxis=dict(title="Z (Å)")
    ),
)

fig2.write_image(output_dir / "pore_network_3d.png", width=800, height=600, scale=2)
print("  Saved: pore_network_3d.png")

# =============================================================================
# Generate Image 3: Connection types
# =============================================================================

print("\nGenerating image: connection_types.png")

internal_connections = [c for c in connections_trimmed if c[2] == "internal"]
periodic_connections = [c for c in connections_trimmed if c[2] == "periodic"]

print(f"  Internal: {len(internal_connections)}")
print(f"  Periodic: {len(periodic_connections)}")

# Create separate traces
internal_traces = mgr.plot_connections_as_it_is(
    np.array(internal_connections),
    all_maxima,
    dgrid.shape,
    mat_atoms,
    cscale=cscale,
    plot_disconnected=False,
    show_spheres=False,
    name="Internal",
)

# Use a proper colorscale for periodic
red_cscale = ["rgb(255, 200, 200)", "rgb(255, 50, 50)"]
periodic_traces = mgr.plot_connections_as_it_is(
    np.array(periodic_connections),
    all_maxima,
    dgrid.shape,
    mat_atoms,
    cscale=red_cscale,
    plot_disconnected=False,
    show_spheres=False,
    name="Periodic",
)

fig3 = go.Figure(
    data=[cell_trace, maxima_trace, *internal_traces, *periodic_traces, *frame_trace],
    layout=mgr.white_layout,
)

fig3.update_layout(
    scene_aspectmode="data",
    title=dict(
        text=f"IRMOF-1: Internal ({len(internal_connections)}) vs Periodic ({len(periodic_connections)}) Connections",
        x=0.5,
        xanchor="center",
        font=dict(size=14, color="black"),
    ),
    scene=dict(
        xaxis=dict(title="X (Å)"), yaxis=dict(title="Y (Å)"), zaxis=dict(title="Z (Å)")
    ),
)

fig3.write_image(output_dir / "connection_types.png", width=800, height=600, scale=2)
print("  Saved: connection_types.png")

# =============================================================================
# Print summary
# =============================================================================

print("\n" + "=" * 50)
print("PORE GRAPH TUTORIAL IMAGES")
print("=" * 50)
print(f"1. pore_graph_concept.png - Conceptual 2D diagram")
print(f"2. pore_network_3d.png - Actual 3D pore network")
print(f"3. connection_types.png - Internal vs periodic connections")
print(f"\nStatistics:")
print(f"  Unique pores: {len(maxima_pbc)}")
print(f"  Total connections: {len(connections_trimmed)}")
print(f"  Internal: {len(internal_connections)}")
print(f"  Periodic: {len(periodic_connections)}")
print("=" * 50)

print("\n✅ All images generated successfully!")
print(f"Images saved to: {output_dir}")

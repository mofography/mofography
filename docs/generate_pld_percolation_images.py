#!/usr/bin/env python3
"""
Generate images for Pore Limiting Diameter (PLD) and Percolation tutorial.
Creates visualizations of PLD concepts, LCD/PLD/DFSD comparison, percolating clusters,
and percolation history.
"""

import os
import sys
import numpy as np
from pathlib import Path

# Add repo root to path
script_path = Path(__file__).resolve()
docs_dir = script_path.parent
repo_root = docs_dir.parent

sys.path.insert(0, str(repo_root))

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch, FancyArrowPatch
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import matplotlib.patches as mpatches

print("=" * 70)
print("PLD AND PERCOLATION TUTORIAL - IMAGE GENERATION")
print("=" * 70)

# =============================================================================
# Setup paths
# =============================================================================

output_dir = docs_dir / "tutorial_images"
output_dir.mkdir(exist_ok=True)

# =============================================================================
# Helper functions for drawing
# =============================================================================

def draw_atom(ax, x, y, r, color='lightblue', edge_color='black', alpha=0.7):
    """Draw an atom (circle)"""
    circle = Circle((x, y), r, facecolor=color, edgecolor=edge_color,
                    linewidth=2, alpha=alpha)
    ax.add_patch(circle)
    return circle

def draw_probe(ax, x, y, r, color='red', alpha=0.4, label=None):
    """Draw a probe molecule"""
    circle = Circle((x, y), r, facecolor=color, edgecolor='darkred',
                    linewidth=2, alpha=alpha, linestyle='--')
    ax.add_patch(circle)
    return circle

# =============================================================================
# Image 1: PLD Concept
# =============================================================================

print("\n[1/7] Generating pld_concept.png...")

fig, ax = plt.subplots(1, 1, figsize=(12, 8))

# Draw unit cell boundary
cell = Rectangle((0, 0), 10, 10, linewidth=3, edgecolor='black',
                 facecolor='#f8f8fa', zorder=0)
ax.add_patch(cell)

# Draw atoms to create a pore network
atom_positions = [
    (1, 1, 0.4), (1, 9, 0.4), (9, 1, 0.4), (9, 9, 0.4),  # Corners
    (3, 3, 0.5), (3, 7, 0.5), (7, 3, 0.5), (7, 7, 0.5),  # Inner ring
    (5, 5, 0.4),  # Center
]

for x, y, r in atom_positions:
    draw_atom(ax, x, y, r, color='lightblue', edge_color='navy', alpha=0.6)

# Draw the pore limiting diameter sphere
pld_radius = 1.0
draw_probe(ax, 5, 5, pld_radius, color='red', alpha=0.3)

# Draw path arrows showing the probe can pass through
for y_pos in [1, 3, 5, 7, 9]:
    if y_pos not in [1, 3, 5, 7, 9]:
        continue
    # X-direction arrow
    ax.annotate('', xy=(9.5, y_pos), xytext=(0.5, y_pos),
                arrowprops=dict(arrowstyle='->', color='green', lw=2))

# Add labels
ax.text(5, -0.5, 'Unit Cell', ha='center', fontsize=14, fontweight='bold')
ax.text(5, 0.8, 'PLD = 2.0 Å', ha='center', fontsize=12,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
ax.text(5, 10.5, 'Probe can traverse from one face to opposite',
        ha='center', fontsize=12)

# Legend
legend_elements = [
    mpatches.Patch(color='lightblue', label='Atoms'),
    mpatches.Patch(color='red', alpha=0.5, label='Probe (r = 1.0 Å)'),
    plt.Line2D([0], [0], color='green', lw=2, label='Path')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

ax.set_xlim(-1, 11)
ax.set_ylim(-2, 12)
ax.set_aspect('equal')
ax.set_title('Pore Limiting Diameter (PLD)\nLargest spherical probe passing from one face to opposite',
             fontsize=14, fontweight='bold')
ax.axis('off')

plt.tight_layout()
fig.savefig(output_dir / "pld_concept.png", dpi=150, bbox_inches='tight',
            facecolor='white')
plt.close()
print(f"  ✅ Saved: {output_dir / 'pld_concept.png'}")

# =============================================================================
# Image 2: LCD vs PLD vs DFSD comparison
# =============================================================================

print("\n[2/7] Generating pld_lcd_dfsd_comparison.png...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# LCD - Largest cavity anywhere
ax = axes[0]
cell = Rectangle((0, 0), 10, 10, linewidth=3, edgecolor='black',
                 facecolor='#f8f8fa', zorder=0)
ax.add_patch(cell)

# Atoms with one large cavity
atom_positions = [
    (2, 2, 0.4), (2, 8, 0.4), (8, 2, 0.4), (8, 8, 0.4),
    (5, 5, 0.3),  # Small obstruction
]
for x, y, r in atom_positions:
    draw_atom(ax, x, y, r, color='lightblue', edge_color='navy', alpha=0.6)

# LCD sphere - largest cavity
draw_probe(ax, 2, 2, 2.5, color='blue', alpha=0.3)
draw_probe(ax, 8, 8, 2.5, color='blue', alpha=0.3)

ax.text(5, -0.5, 'LCD', ha='center', fontsize=14, fontweight='bold', color='blue')
ax.text(5, 10.5, 'Largest Cavity Diameter', ha='center', fontsize=11)
ax.text(5, 0.5, 'max(distance_grid) × 2', ha='center', fontsize=10, style='italic')
ax.set_xlim(-1, 11)
ax.set_ylim(-2, 12)
ax.set_aspect('equal')
ax.axis('off')

# PLD - largest sphere that can traverse
ax = axes[1]
cell = Rectangle((0, 0), 10, 10, linewidth=3, edgecolor='black',
                 facecolor='#f8f8fa', zorder=0)
ax.add_patch(cell)

# Atoms with narrow constriction
atom_positions = [
    (2, 2, 0.4), (2, 8, 0.4), (8, 2, 0.4), (8, 8, 0.4),
    (5, 4.5, 0.6), (5, 5.5, 0.6),  # Narrow window
]
for x, y, r in atom_positions:
    draw_atom(ax, x, y, r, color='lightblue', edge_color='navy', alpha=0.6)

# PLD sphere - limited by narrow window
draw_probe(ax, 2, 2, 1.8, color='red', alpha=0.3)
draw_probe(ax, 5, 5, 1.0, color='red', alpha=0.4)  # Narrowest point

ax.text(5, -0.5, 'PLD', ha='center', fontsize=14, fontweight='bold', color='red')
ax.text(5, 10.5, 'Pore Limiting Diameter', ha='center', fontsize=11)
ax.text(5, 0.5, 'max percolating radius × 2', ha='center', fontsize=10, style='italic')
ax.set_xlim(-1, 11)
ax.set_ylim(-2, 12)
ax.set_aspect('equal')
ax.axis('off')

# DFSD - largest sphere in percolating regions
ax = axes[2]
cell = Rectangle((0, 0), 10, 10, linewidth=3, edgecolor='black',
                 facecolor='#f8f8fa', zorder=0)
ax.add_patch(cell)

# Atoms with large non-percolating cavity
atom_positions = [
    (2, 2, 0.4), (2, 8, 0.4), (8, 2, 0.4), (8, 8, 0.4),
    (2, 5, 0.6), (3, 5, 0.6),  # Narrow connection blocking access to one cavity
]
for x, y, r in atom_positions:
    draw_atom(ax, x, y, r, color='lightblue', edge_color='navy', alpha=0.6)

# DFSD sphere - only in percolating regions
draw_probe(ax, 2, 2, 1.8, color='green', alpha=0.3)
draw_probe(ax, 2, 8, 1.8, color='green', alpha=0.3)

ax.text(5, -0.5, 'DFSD', ha='center', fontsize=14, fontweight='bold', color='green')
ax.text(5, 10.5, 'Directional Free Sphere Diameter', ha='center', fontsize=11)
ax.text(5, 0.5, 'max in percolating regions × 2', ha='center', fontsize=10, style='italic')
ax.set_xlim(-1, 11)
ax.set_ylim(-2, 12)
ax.set_aspect('equal')
ax.axis('off')

plt.suptitle('Three Pore Diameter Metrics', fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout()
fig.savefig(output_dir / "pld_lcd_dfsd_comparison.png", dpi=150, bbox_inches='tight',
            facecolor='white')
plt.close()
print(f"  ✅ Saved: {output_dir / 'pld_lcd_dfsd_comparison.png'}")

# =============================================================================
# Image 3: Percolation History visualization
# =============================================================================

print("\n[3/7] Generating percolation_history.png...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

probe_sizes = [1.0, 3.0, 5.0]  # Small, medium, large
titles = ['Small Probe\n(1.0 Å)', 'Medium Probe\n(3.0 Å)', 'Large Probe\n(5.0 Å)\n(PLD)']

for idx, (ax, probe_r, title) in enumerate(zip(axes, probe_sizes, titles)):
    # Unit cell
    cell = Rectangle((0, 0), 10, 10, linewidth=3, edgecolor='black',
                     facecolor='#f8f8fa', zorder=0)
    ax.add_patch(cell)

    # Atoms creating a network with varying window sizes
    atom_positions = [
        (1, 1, 0.4), (1, 9, 0.4), (9, 1, 0.4), (9, 9, 0.4),
        (3, 3, 0.4), (3, 7, 0.4), (7, 3, 0.4), (7, 7, 0.4),
        (5, 5, 0.4),
        # Narrow windows
        (2, 5, 0.7), (8, 5, 0.7), (5, 2, 0.7), (5, 8, 0.7),
    ]

    for x, y, r in atom_positions:
        draw_atom(ax, x, y, r, color='lightblue', edge_color='navy', alpha=0.6)

    # Draw accessible regions based on probe size
    if probe_r <= 1.5:
        # Small probe: all regions accessible
        accessible_regions = [(2, 2, 1.5), (8, 8, 1.5), (2, 8, 1.5), (8, 2, 1.5),
                             (5, 5, 1.5), (3, 3, 1.2), (7, 7, 1.2)]
        color = 'lightgreen'
        n_clusters = 1
    elif probe_r <= 3.5:
        # Medium probe: some regions inaccessible
        accessible_regions = [(2, 2, 1.5), (8, 8, 1.5), (3, 3, 1.2), (7, 7, 1.2)]
        color = 'yellow'
        n_clusters = 2
    else:
        # Large probe: only narrow path
        accessible_regions = [(2, 2, 1.5), (8, 8, 1.5)]
        color = 'orange'
        n_clusters = 2

    for x, y, r in accessible_regions:
        draw_probe(ax, x, y, r, color=color, alpha=0.4)

    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.text(5, -0.8, f'{n_clusters} cluster(s)', ha='center', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_xlim(-1, 11)
    ax.set_ylim(-2, 12)
    ax.set_aspect('equal')
    ax.axis('off')

plt.suptitle('Percolation History: How Accessibility Changes with Probe Size',
             fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
fig.savefig(output_dir / "percolation_history.png", dpi=150, bbox_inches='tight',
            facecolor='white')
plt.close()
print(f"  ✅ Saved: {output_dir / 'percolation_history.png'}")

# =============================================================================
# Image 4: Winding number concept
# =============================================================================

print("\n[4/7] Generating winding_number_concept.png...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

winding_examples = [
    ((0, 0, 0), 'Contractible\nLoop\nNo Percolation', 'gray'),
    ((1, 0, 0), 'Winding (1,0,0)\nX Percolation', 'green'),
    ((0, 1, 1), 'Winding (0,1,1)\nY & Z Percolation', 'blue')
]

for ax, (winding, title, color) in zip(axes, winding_examples):
    # Unit cell with grid
    for i in range(3):
        for j in range(3):
            cell = Rectangle((i*3, j*3), 3, 3, linewidth=1, edgecolor='lightgray',
                             facecolor='white', zorder=0)
            ax.add_patch(cell)

    # Draw nodes
    nodes = [(1.5, 1.5), (4.5, 1.5), (7.5, 1.5),
             (1.5, 4.5), (4.5, 4.5), (7.5, 4.5),
             (1.5, 7.5), (4.5, 7.5), (7.5, 7.5)]

    for x, y in nodes:
        draw_atom(ax, x, y, 0.3, color='lightblue', edge_color='navy', alpha=0.8)

    # Draw edges based on winding
    if winding == (0, 0, 0):
        # Contractible loop
        path_nodes = [(1.5, 1.5), (4.5, 1.5), (4.5, 4.5), (1.5, 4.5), (1.5, 1.5)]
    elif winding == (1, 0, 0):
        # Wraps in X
        path_nodes = [(1.5, 1.5), (4.5, 1.5), (7.5, 1.5), (7.5, 4.5),
                     (4.5, 4.5), (1.5, 4.5), (1.5, 1.5)]
    else:  # (0, 1, 1)
        # Wraps in Y and Z
        path_nodes = [(1.5, 1.5), (4.5, 1.5), (7.5, 1.5), (7.5, 4.5),
                     (7.5, 7.5), (4.5, 7.5), (1.5, 7.5), (1.5, 4.5), (1.5, 1.5)]

    # Draw path
    for i in range(len(path_nodes) - 1):
        x1, y1 = path_nodes[i]
        x2, y2 = path_nodes[i + 1]
        ax.plot([x1, x2], [y1, y2], color=color, linewidth=3, alpha=0.7)

    # Mark start/end
    ax.plot(path_nodes[0][0], path_nodes[0][1], 'ro', markersize=10, zorder=10)

    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.text(4.5, -0.8, f'Winding: {winding}', ha='center', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 9)
    ax.set_aspect('equal')
    ax.axis('off')

plt.suptitle('Winding Numbers: Detecting Non-Contractible Loops',
             fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
fig.savefig(output_dir / "winding_number_concept.png", dpi=150, bbox_inches='tight',
            facecolor='white')
plt.close()
print(f"  ✅ Saved: {output_dir / 'winding_number_concept.png'}")

# =============================================================================
# Image 5: Percolating clusters visualization (schematic)
# =============================================================================

print("\n[5/7] Generating percolating_clusters_example.png...")

fig, ax = plt.subplots(1, 1, figsize=(10, 10))

# Unit cell
cell = Rectangle((0, 0), 10, 10, linewidth=3, edgecolor='black',
                 facecolor='#f8f8fa', zorder=0)
ax.add_patch(cell)

# Create different regions (clusters)
# Cluster 1 (blue): percolating in X
cluster1_regions = [(2, 5), (3, 5), (4, 5)]
for x, y in cluster1_regions:
    draw_probe(ax, x, y, 1.0, color='blue', alpha=0.5)
    draw_atom(ax, x, y, 0.2, color='blue', edge_color='darkblue', alpha=0.3)

# Cluster 2 (green): percolating in Y
cluster2_regions = [(5, 2), (5, 3), (5, 4)]
for x, y in cluster2_regions:
    draw_probe(ax, x, y, 1.0, color='green', alpha=0.5)
    draw_atom(ax, x, y, 0.2, color='green', edge_color='darkgreen', alpha=0.3)

# Cluster 3 (red): percolating in Z (out of page)
cluster3_regions = [(7, 7), (8, 7), (7, 8)]
for x, y in cluster3_regions:
    draw_probe(ax, x, y, 1.0, color='red', alpha=0.5)
    draw_atom(ax, x, y, 0.2, color='red', edge_color='darkred', alpha=0.3)

# Inaccessible regions (gray)
inaccessible = [(2, 2), (8, 8)]
for x, y in inaccessible:
    draw_atom(ax, x, y, 0.5, color='gray', edge_color='black', alpha=0.4)

# Labels
ax.text(3, 6, 'Cluster 1\n(X-percolating)', ha='center', fontsize=10,
        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
ax.text(6, 3, 'Cluster 2\n(Y-percolating)', ha='center', fontsize=10,
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
ax.text(7.5, 8.5, 'Cluster 3\n(Z-percolating)', ha='center', fontsize=10,
        bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))

ax.set_xlim(-1, 11)
ax.set_ylim(-1, 11)
ax.set_aspect('equal')
ax.set_title('Percolating Clusters at PLD Threshold\nEach color = distinct percolating cluster',
             fontsize=14, fontweight='bold')
ax.axis('off')

plt.tight_layout()
fig.savefig(output_dir / "percolating_clusters_example.png", dpi=150, bbox_inches='tight',
            facecolor='white')
plt.close()
print(f"  ✅ Saved: {output_dir / 'percolating_clusters_example.png'}")

# =============================================================================
# Image 6: Percolation types (1D, 2D, 3D)
# =============================================================================

print("\n[6/7] Generating percolation_types.png...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

configs = [
    ('1D Channels', [('Y', True), ('X', False), ('Z', False)], 'green'),
    ('2D Layers', [('X', True), ('Y', True), ('Z', False)], 'blue'),
    ('3D Network', [('X', True), ('Y', True), ('Z', True)], 'purple')
]

for idx, (ax, (title, axes_status, color)) in enumerate(zip(axes, configs)):
    # Unit cell
    cell = Rectangle((0, 0), 10, 10, linewidth=3, edgecolor='black',
                     facecolor='#f8f8fa', zorder=0)
    ax.add_patch(cell)

    # Draw atoms based on percolation type
    if title == '1D Channels':
        # Y-direction channels
        for y in [2, 5, 8]:
            draw_probe(ax, 5, y, 1.2, color=color, alpha=0.4)
        atom_positions = [(3, 3, 0.6), (7, 3, 0.6), (3, 7, 0.6), (7, 7, 0.6)]
    elif title == '2D Layers':
        # XY plane
        draw_probe(ax, 5, 5, 2.0, color=color, alpha=0.4)
        atom_positions = [(5, 2, 0.6), (5, 8, 0.6)]
    else:  # 3D
        # All directions
        for x, y in [(2, 2), (2, 8), (8, 2), (8, 8), (5, 5)]:
            draw_probe(ax, x, y, 1.5, color=color, alpha=0.4)
        atom_positions = []

    for x, y, r in atom_positions:
        draw_atom(ax, x, y, r, color='lightblue', edge_color='navy', alpha=0.6)

    # Percolation indicators
    y_offset = -1.5
    for axis, status in axes_status:
        symbol = '✅' if status else '❌'
        ax.text(2 + idx*0.1, y_offset, f'{symbol} {axis}',
                fontsize=12, transform=ax.transData)
        y_offset -= 0.6

    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlim(-1, 11)
    ax.set_ylim(-4, 11)
    ax.set_aspect('equal')
    ax.axis('off')

plt.suptitle('Types of Percolation in Porous Materials',
             fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
fig.savefig(output_dir / "percolation_types.png", dpi=150, bbox_inches='tight',
            facecolor='white')
plt.close()
print(f"  ✅ Saved: {output_dir / 'percolation_types.png'}")

# =============================================================================
# Image 7: Summary diagram showing the workflow
# =============================================================================

print("\n[7/7] generating pld_workflow.png...")

fig, ax = plt.subplots(1, 1, figsize=(14, 8))

# Define boxes for workflow
boxes = [
    (0.5, 0.85, 'Load Structure\n(CIF file)', 'lightblue'),
    (0.5, 0.70, 'Compute Distance Grid', 'lightgreen'),
    (0.5, 0.55, 'Periodic Segmentation\n→ Pore Graph', 'lightyellow'),
    (0.5, 0.40, 'Percolation Analysis\n(Winding Numbers)', 'lightcoral'),
    (0.5, 0.25, 'Compute LCD, PLD, DFSD', 'lightpink'),
    (0.5, 0.10, 'Visualize & Interpret', 'lavender'),
]

# Draw boxes
for x, y, text, color in boxes:
    box = FancyBboxPatch((x - 0.15, y - 0.06), 0.3, 0.12,
                         boxstyle="round,pad=0.02",
                         facecolor=color, edgecolor='black', linewidth=2)
    ax.add_patch(box)
    ax.text(x, y, text, ha='center', va='center', fontsize=11, fontweight='bold')

# Draw arrows
for i in range(len(boxes) - 1):
    y_start = boxes[i][1] - 0.06
    y_end = boxes[i + 1][1] + 0.06
    arrow = FancyArrowPatch((0.5, y_start), (0.5, y_end),
                           arrowstyle='->', mutation_scale=20, linewidth=2,
                           color='gray')
    ax.add_patch(arrow)

# Add side notes
ax.text(0.75, 0.70, 'Atoms object\nUnit cell', fontsize=9, style='italic',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
ax.text(0.75, 0.55, 'Regions, Maxima\nConnections', fontsize=9, style='italic',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
ax.text(0.75, 0.40, 'Percolation history\nWinding loops', fontsize=9, style='italic',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
ax.text(0.75, 0.25, 'LCD = max(dgrid) × 2\nPLD = max percolating × 2\nDFSD = max(percolating) × 2',
        fontsize=9, style='italic',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
ax.set_title('PLD Computation Workflow', fontsize=16, fontweight='bold', y=0.98)

plt.tight_layout()
fig.savefig(output_dir / "pld_workflow.png", dpi=150, bbox_inches='tight',
            facecolor='white')
plt.close()
print(f"  ✅ Saved: {output_dir / 'pld_workflow.png'}")

# =============================================================================
# Summary
# =============================================================================

print("\n" + "=" * 70)
print("IMAGE GENERATION COMPLETE")
print("=" * 70)
print(f"\nAll images saved to: {output_dir}")
print("\nGenerated images:")
print("  1. pld_concept.png - What is PLD?")
print("  2. pld_lcd_dfsd_comparison.png - Three pore diameter metrics")
print("  3. percolation_history.png - How accessibility changes with probe size")
print("  4. winding_number_concept.png - How winding numbers detect percolation")
print("  5. percolating_clusters_example.png - Visualizing percolating clusters")
print("  6. percolation_types.png - 1D, 2D, and 3D percolation")
print("  7. pld_workflow.png - Complete workflow diagram")
print("\n" + "=" * 70)
#%%
import mofography as mgr
from ase.io import read
import plotly.graph_objects as go
from pathlib import Path
import numpy as np
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
material = 'NU-1000_S1.cif'
path_to_cif = Path('~/repos/mofography/examples/TestCIFs/{}'.format(material)).expanduser()
print(path_to_cif.exists())
# %%
spacing = 0.5
h=0.2
mask_thickness = 0.5
threshold_abs = 1.5
extra_padding = 5.0
return_dgrid = False

mat_atoms = read(path_to_cif)
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing).astype(np.float32)

results = mgr.do_periodic_segmentation_rust(dgrid,  mat_atoms, spacing = spacing, mask_thickness=mask_thickness, h=h, threshold_abs=threshold_abs, extra_padding=extra_padding, return_dgrid=return_dgrid)

# %%
lcd , pld, 
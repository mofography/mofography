#!/usr/bin/env python3
"""
Generate images for Pore Size Distributions (Gelb & Gubbins) tutorial.
Creates conceptual diagrams and output visualizations for both global and pocket-wise PSD.
"""

import os
import sys
import numpy as np
from pathlib import Path

# Add repo root to path
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
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch

print("=" * 60)
print("PORE SIZE DISTRIBUTIONS - IMAGE GENERATION")
print("=" * 60)

# =============================================================================
# Setup paths and parameters
# =============================================================================

output_dir = script_path.parent / "tutorial_images"
output_dir.mkdir(exist_ok=True)

path_to_cif = repo_root / "examples" / "TestCIFs" / "NU-1000_S1.cif"

# Check if CIF exists
if not path_to_cif.exists():
    print(f"ERROR: CIF not found: {path_to_cif}")
    sys.exit(1)

# =============================================================================
# Load structure and run workflow
# =============================================================================

print("\nLoading structure...")
mat_atoms = read(str(path_to_cif))

spacing = 0.25  # Fine spacing for accurate PSD
dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=spacing)

print(f"Distance grid shape: {dgrid.shape}")
print(f"  Spacing: {spacing} Å")

# Run periodic segmentation
print("\nRunning periodic segmentation...")
results = mgr.do_periodic_segmentation_rust(
    dgrid=dgrid,
    mat_atoms=mat_atoms,
    spacing=spacing,
    mask_thickness=1.0,
    h=0.2,
    threshold_abs=1.5,
    minimum_window_separation=1.5,
    wall_thickness=0.08,
    extra_padding=5.0,
    compute_connections=True,
    return_padded_regions=False,
)

regions_pbc = results["regions_pbc"]
print(f"Found {len(np.unique(regions_pbc))} unique regions")

# =============================================================================
# Compute PSDs
# =============================================================================

print("\nComputing global PSD (uniform sampling)...")
num_samples_total = 10000
sample_points_uniform, pore_sizes_uniform = mgr.compute_pore_sizes_uniform_rust(
    dgrid, mat_atoms, num_samples_total=num_samples_total
)
print(f"  Computed {len(pore_sizes_uniform)} pore sizes")

print("\nComputing pocket-wise PSD...")
num_samples_per_region = 1000
sample_points_pocketwise, pore_sizes_pocketwise, region_labels = (
    mgr.compute_pore_sizes_pocketwise_rust(
        dgrid, regions_pbc, mat_atoms, samples_per_region=num_samples_per_region
    )
)
print(f"  Computed pore sizes for {len(np.unique(region_labels))} pockets")

# =============================================================================
# Generate Image 1: Conceptual diagram of Gelb & Gubbins method
# =============================================================================

print("\nGenerating Image 1: psd_concept.png (Conceptual diagram)")

fig1, ax1 = plt.subplots(1, 1, figsize=(12, 8))

# Draw unit cell
unit_cell = FancyBboxPatch(
    (2, 2),
    8,
    8,
    boxstyle="round,pad=0.1",
    linewidth=3,
    edgecolor="black",
    facecolor="#f8f8fa",
)
ax1.add_patch(unit_cell)

# Draw framework atoms (spheres)
np.random.seed(42)
n_framework_atoms = 8
atom_positions = []
for i in range(n_framework_atoms):
    x = 3 + np.random.rand() * 6
    y = 3 + np.random.rand() * 6
    atom_positions.append((x, y))
    atom_circle = patches.Circle(
        (x, y), 0.3, facecolor="#e74c3c", edgecolor="black", linewidth=2, alpha=0.8
    )
    ax1.add_patch(atom_circle)

# Draw sample points and their largest inscribed spheres
np.random.seed(123)
n_samples = 5
sample_positions = []
for i in range(n_samples):
    x = 3 + np.random.rand() * 6
    y = 3 + np.random.rand() * 6

    # Check if point is in void space (away from atoms)
    min_dist = min(
        [np.sqrt((x - ax[0]) ** 2 + (y - ax[1]) ** 2) for ax in atom_positions]
    )

    if min_dist > 1.5:  # In void space
        sample_positions.append((x, y, min_dist))

        # Draw largest inscribed sphere
        sphere = patches.Circle(
            (x, y),
            min_dist,
            facecolor="#66c2a5",
            edgecolor="black",
            linewidth=1.5,
            alpha=0.4,
        )
        ax1.add_patch(sphere)
        ax1.text(
            x,
            y,
            f"r={min_dist:.1f}",
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
        )
        ax1.plot(x, y, "rx", color="#ff6b6b", markersize=8, markeredgewidth=2)

# Add annotations
ax1.text(
    6,
    -0.5,
    "Framework atoms (excluded regions)",
    fontsize=11,
    ha="center",
    fontweight="bold",
)
ax1.text(
    6, 11, "Sample points in void space", fontsize=11, ha="center", fontweight="bold"
)
ax1.text(2.5, 5.5, "Largest inscribed\nsphere", fontsize=9, ha="center", rotation=45)

ax1.set_xlim(0, 12)
ax1.set_ylim(0, 12)
ax1.set_aspect("equal")
ax1.set_title(
    "Gelb & Gubbins Method: Largest Inscribed Spheres", fontsize=14, fontweight="bold"
)
ax1.axis("off")

plt.tight_layout()
fig1.savefig(
    output_dir / "psd_concept.png", dpi=150, bbox_inches="tight", facecolor="white"
)
plt.close()
print(f"  Saved: {output_dir / 'psd_concept.png'}")

# =============================================================================
# Generate Image 2: Global PSD output
# =============================================================================

print("\nGenerating Image 2: psd_global_output.png (Global PSD histogram)")

bins = np.linspace(0, 50, 200)
hist, bin_edges = np.histogram(pore_sizes_uniform * 2, bins=bins, density=False)
bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2

fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6))

ax2.bar(
    bin_centers,
    hist,
    width=bins[1] - bins[0],
    edgecolor="black",
    alpha=0.8,
    color="#66c2a5",
)
ax2.set_xlabel("Pore Diameter (Å)", fontsize=13)
ax2.set_ylabel("Counts", fontsize=13)
ax2.set_title(
    f"Global PSD: Uniform Sampling Across Unit Cell\n(NU-1000, {num_samples_total} samples)",
    fontsize=14,
    fontweight="bold",
)
ax2.grid(True, alpha=0.3, linestyle="--")
ax2.tick_params(labelsize=11)

# Add statistics
mean_pore = np.mean(pore_sizes_uniform * 2)
std_pore = np.std(pore_sizes_uniform * 2)
max_pore = np.max(pore_sizes_uniform * 2)
ax2.axvline(
    mean_pore,
    color="red",
    linestyle="--",
    linewidth=2,
    alpha=0.7,
    label=f"Mean: {mean_pore:.1f} Å",
)
ax2.axvline(mean_pore + std_pore, color="blue", linestyle=":", linewidth=1.5, alpha=0.6)
ax2.axvline(mean_pore - std_pore, color="blue", linestyle=":", linewidth=1.5, alpha=0.6)
ax2.legend(fontsize=10, loc="upper right")

plt.tight_layout()
fig2.savefig(
    output_dir / "psd_global_output.png",
    dpi=150,
    bbox_inches="tight",
    facecolor="white",
)
plt.close()
print(f"  Saved: {output_dir / 'psd_global_output.png'}")

# =============================================================================
# Generate Image 3: Pocket-wise concept
# =============================================================================

print("\nGenerating Image 3: pocketwise_concept.png (Pocket-wise sampling)")

fig3, ax3 = plt.subplots(1, 1, figsize=(12, 8))

# Draw unit cell
unit_cell3 = FancyBboxPatch(
    (2, 2),
    8,
    8,
    boxstyle="round,pad=0.1",
    linewidth=3,
    edgecolor="black",
    facecolor="#f8f8fa",
)
ax3.add_patch(unit_cell3)

# Draw colored regions (conceptual)
colors_pockets = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"]
np.random.seed(456)
n_pockets = 4
pocket_centers = []
for i in range(n_pockets):
    x = 3 + np.random.rand() * 6
    y = 3 + np.random.rand() * 6
    pocket_centers.append((x, y))
    color = colors_pockets[i % len(colors_pockets)]

    # Draw pocket region (rough shape)
    pocket = patches.Ellipse(
        (x, y),
        1.8,
        1.4,
        angle=np.random.rand() * 180,
        facecolor=color,
        edgecolor="black",
        alpha=0.3,
    )
    ax3.add_patch(pocket)

    # Draw sample points in this pocket
    for _ in range(3):
        sx = x + (np.random.rand() - 0.5) * 2
        sy = y + (np.random.rand() - 0.5) * 2
        ax3.plot(
            sx,
            sy,
            "o",
            color=color,
            markersize=6,
            markeredgecolor="black",
            markeredgewidth=1.5,
        )

# Add labels
ax3.text(
    6,
    -0.5,
    "Different pocket types (colored regions)",
    fontsize=11,
    ha="center",
    fontweight="bold",
)
ax3.text(6, 11, "Sample points per pocket", fontsize=11, ha="center", fontweight="bold")

ax3.set_xlim(0, 12)
ax3.set_ylim(0, 12)
ax3.set_aspect("equal")
ax3.set_title(
    "Pocket-wise PSD: Sampling Within Each Pore Region", fontsize=14, fontweight="bold"
)
ax3.axis("off")

plt.tight_layout()
fig3.savefig(
    output_dir / "pocketwise_concept.png",
    dpi=150,
    bbox_inches="tight",
    facecolor="white",
)
plt.close()
print(f"  Saved: {output_dir / 'pocketwise_concept.png'}")

# =============================================================================
# Generate Image 4: Pocket-wise PSD output
# =============================================================================

print("\nGenerating Image 4: psd_pocketwise_output.png (Pocket-wise PSD curves)")

# Compute PSDs for each pocket
from mofography.helper import compute_psds_from_points_pocketwise

pocket_psds = compute_psds_from_points_pocketwise(
    region_labels, pore_sizes_pocketwise * 2, bins=bins
)

fig4, ax4 = plt.subplots(1, 1, figsize=(12, 7))

# Plot overall histogram (light gray)
ax4.bar(
    bin_centers,
    hist,
    width=bins[1] - bins[0],
    edgecolor="black",
    alpha=0.3,
    color="gray",
    label="Overall (uniform sampling)",
)

# Plot individual pocket PSDs
colors_curve = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"]
for i, (label, psd_data) in enumerate(pocket_psds.items()):
    if i >= len(colors_curve):
        break
    ax4.plot(
        psd_data["bin_centers"],
        psd_data["hist"],
        linewidth=2.5,
        color=colors_curve[i],
        label=f"Pocket {label}",
        linestyle="-",
        alpha=0.8,
    )

ax4.set_xlabel("Pore Diameter (Å)", fontsize=13)
ax4.set_ylabel("Counts", fontsize=13)
ax4.set_title(
    f"Pocket-wise PSD: Individual Pocket Distributions\n(NU-1000, {len(pocket_psds)} pockets, {num_samples_per_region} samples each)",
    fontsize=14,
    fontweight="bold",
)
ax4.set_xlim(0, 50)
ax4.legend(fontsize=10, loc="upper right")
ax4.grid(True, alpha=0.3, linestyle="--")
ax4.tick_params(labelsize=11)

plt.tight_layout()
fig4.savefig(
    output_dir / "psd_pocketwise_output.png",
    dpi=150,
    bbox_inches="tight",
    facecolor="white",
)
plt.close()
print(f"  Saved: {output_dir / 'psd_pocketwise_output.png'}")

# =============================================================================
# Generate Image 5: Global vs Pocket-wise comparison
# =============================================================================

print("\nGenerating Image 5: psd_comparison.png (Global vs Pocket-wise)")

fig5, (ax5a, ax5b) = plt.subplots(1, 2, figsize=(14, 6))

# Left: Global PSD
ax5a.bar(
    bin_centers,
    hist,
    width=bins[1] - bins[0],
    edgecolor="black",
    alpha=0.7,
    color="#66c2a5",
)
ax5a.set_xlabel("Pore Diameter (Å)", fontsize=12)
ax5a.set_ylabel("Counts", fontsize=12)
ax5a.set_title(
    "Global PSD\n(Uniform sampling across unit cell)", fontsize=13, fontweight="bold"
)
ax5a.grid(True, alpha=0.3)
ax5a.tick_params(labelsize=10)

# Right: Pocket-wise PSD
ax5b.bar(
    bin_centers,
    hist,
    width=bins[1] - bins[0],
    edgecolor="black",
    alpha=0.3,
    color="gray",
    label="Overall",
)
for i, (label, psd_data) in enumerate(pocket_psds.items()):
    if i >= len(colors_curve):
        break
    ax5b.plot(
        psd_data["bin_centers"],
        psd_data["hist"],
        linewidth=2,
        color=colors_curve[i],
        label=f"Pocket {label}",
        linestyle="-",
        alpha=0.8,
    )

ax5b.set_xlabel("Pore Diameter (Å)", fontsize=12)
ax5b.set_title(
    "Pocket-wise PSD\n(Per-pocket distributions)", fontsize=13, fontweight="bold"
)
ax5b.legend(fontsize=9, loc="upper right")
ax5b.grid(True, alpha=0.3)
ax5b.tick_params(labelsize=10)

plt.suptitle(
    "Global vs Pocket-wise PSD Comparison", fontsize=15, fontweight="bold", y=0.98
)
plt.tight_layout()
fig5.savefig(
    output_dir / "psd_comparison.png", dpi=150, bbox_inches="tight", facecolor="white"
)
plt.close()
print(f"  Saved: {output_dir / 'psd_comparison.png'}")

# =============================================================================
# Summary
# =============================================================================

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"\nImages generated in: {output_dir}")
print(f"  1. psd_concept.png - Gelb & Gubbins method conceptual diagram")
print(f"  2. psd_global_output.png - Global PSD histogram (NU-1000)")
print(f"  3. pocketwise_concept.png - Pocket-wise sampling conceptual diagram")
print(f"  4. psd_pocketwise_output.png - Individual pocket PSD curves")
print(f"  5. psd_comparison.png - Global vs Pocket-wise comparison")
print("\nAll images generated successfully!")
print("=" * 60)

#!/usr/bin/env python3
"""Deep-dive: compare watershed label arrays between backends at multiple mask thresholds.
Reports pixel-level differences and their spatial locations."""

import sys, time, numpy as np

sys.path.insert(0, "/Users/arungopalan/repos/mofography")
import mofography as mgr
import mgr_rust
from ase.io import read
from ase.data import covalent_radii
from ase.geometry import cellpar_to_cell
from skimage.segmentation import watershed as sk_watershed

# ── Constants ───────────────────────────────────────────────────
CIF_PATH = "/Users/arungopalan/repos/mofography/mofography/examples/TestCIFs/P1-Cu-BTC.cif"
SPACING = 0.5
H = 0.2
THRESHOLD_ABS = 1.5
MIN_DISTANCE = 1.5

# ── Load & compute ──────────────────────────────────────────────
print("Loading & computing dgrid ...")
mat_atoms = read(CIF_PATH)
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
radii = covalent_radii[mat_updated.get_atomic_numbers()]
grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing=SPACING, block_size=20000)
dgrid = grid.compute().astype(np.float32)
cell = cellpar_to_cell(mat_atoms.cell.cellpar()).T

# Detect maxima once
merged, _, _ = mgr_rust.detect_maxima_selective_smoothing_parallel(
    dgrid, preserve_threshold=3.0, smooth_threshold=6.0, max_smooth_radius=2.0,
    h=H, threshold_abs=THRESHOLD_ABS, min_distance=MIN_DISTANCE, cell=cell,
)
markers = np.zeros_like(dgrid, dtype=np.int32)
markers[tuple(merged.T)] = np.arange(1, len(merged) + 1)
print(f"Maxima: {len(merged)}, dgrid range: [{dgrid.min():.3f}, {dgrid.max():.3f}]")

padding = 5
padded_grid = np.pad(dgrid, padding, mode="wrap")
padded_markers = np.pad(markers, padding, mode="wrap")

# ── Test at multiple mask thresholds ────────────────────────────
mask_values = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
print(f"\n{'='*80}")
print(f"{'mask':>6s} | {'diffs rust_bucket vs rust_det':>30s} | {'diffs rust_det vs skimage_c2':>30s} | {'diffs rust_det vs skimage_c3':>30s}")
print(f"{'':>6s} | {'npix %}    nlabels':>30s} | {'npix %}    nlabels':>30s} | {'npix %}    nlabels':>30s}")
print(f"{'-'*80}")

for mask_thick in mask_values:
    padded_mask_u8 = np.pad((dgrid > mask_thick).astype(np.uint8), padding, mode="wrap")
    padded_mask_bool = np.pad(dgrid > mask_thick, padding, mode="wrap")

    # rust_det
    reg_det = mgr_rust.seeded_watershed_det(
        padded_grid.astype(np.float32), padded_markers.astype(np.uint32), padded_mask_u8,
    )
    # rust_bucket
    reg_bkt = mgr_rust.seeded_watershed_bucket(
        padded_grid.astype(np.float32), padded_markers.astype(np.uint32), padded_mask_u8,
    )
    # skimage conn=2
    reg_sk2 = sk_watershed(-padded_grid.astype(float), padded_markers.astype(np.int32),
                           mask=padded_mask_bool, compactness=0, connectivity=2).astype(np.uint32)
    # skimage conn=3
    reg_sk3 = sk_watershed(-padded_grid.astype(float), padded_markers.astype(np.int32),
                           mask=padded_mask_bool, compactness=0, connectivity=3).astype(np.uint32)

    # Compare
    det_vs_bkt = np.sum(reg_det != reg_bkt)
    det_vs_sk2 = np.sum(reg_det != reg_sk2)
    det_vs_sk3 = np.sum(reg_det != reg_sk3)
    total = reg_det.size

    line = f"{mask_thick:6.2f} | "
    for diff_count, label_cmp in [(det_vs_bkt, "bkt_vs_det"), (det_vs_sk2, "sk2_vs_det"), (det_vs_sk3, "sk3_vs_det")]:
        pct = 100 * diff_count / total
        n_l = len(np.unique(eval(f"reg_{label_cmp.split('_vs_')[0]}"))) - 1
        line += f"{diff_count:6d} {pct:5.2f}% {n_l:4d}lbl    "
    print(line)

# ── Detailed comparison at mask where differences appear ────────
print(f"\n{'='*80}")
print("DETAILED: spatial map of differences at mask=0.15")
padded_mask_u8 = np.pad((dgrid > 0.15).astype(np.uint8), padding, mode="wrap")
padded_mask_bool = np.pad(dgrid > 0.15, padding, mode="wrap")

reg_det = mgr_rust.seeded_watershed_det(
    padded_grid.astype(np.float32), padded_markers.astype(np.uint32), padded_mask_u8,
)
reg_sk2 = sk_watershed(-padded_grid.astype(float), padded_markers.astype(np.int32),
                       mask=padded_mask_bool, compactness=0, connectivity=2).astype(np.uint32)
reg_sk3 = sk_watershed(-padded_grid.astype(float), padded_markers.astype(np.int32),
                       mask=padded_mask_bool, compactness=0, connectivity=3).astype(np.uint32)

for label, reg_other in [("skimage_conn2", reg_sk2), ("skimage_conn3", reg_sk3)]:
    diff_mask = reg_det != reg_other
    n_diff = np.sum(diff_mask)
    if n_diff == 0:
        print(f"\n  rust_det vs {label}: IDENTICAL (0 diffs)")
        continue
    print(f"\n  rust_det vs {label}: {n_diff} diffs")
    diff_coords = np.vstack(np.where(diff_mask)).T
    if len(diff_coords) > 20:
        print(f"  First 20 diff locations (z,y,x):")
        for c in diff_coords[:20]:
            print(f"    {tuple(c)}: rust_det={reg_det[tuple(c)]}, {label}={reg_other[tuple(c)]}, dgrid={padded_grid[tuple(c)]:.3f}")

# ── Also try with higher grid resolution ────────────────────────
print(f"\n{'='*80}")
print("Testing with finer grid (spacing=0.25) ...")
grid2 = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing=0.25, block_size=20000)
dgrid2 = grid2.compute().astype(np.float32)
print(f"  dgrid shape: {dgrid2.shape}")

merged2, _, _ = mgr_rust.detect_maxima_selective_smoothing_parallel(
    dgrid2, preserve_threshold=3.0, smooth_threshold=6.0, max_smooth_radius=2.0,
    h=H, threshold_abs=THRESHOLD_ABS, min_distance=MIN_DISTANCE, cell=cell,
)
markers2 = np.zeros_like(dgrid2, dtype=np.int32)
markers2[tuple(merged2.T)] = np.arange(1, len(merged2) + 1)
print(f"  Maxima: {len(merged2)}, dgrid range: [{dgrid2.min():.3f}, {dgrid2.max():.3f}]")

padded_grid2 = np.pad(dgrid2, padding, mode="wrap")
padded_markers2 = np.pad(markers2, padding, mode="wrap")

for mask_thick in [0.05, 0.15]:
    padded_mask_u8 = np.pad((dgrid2 > mask_thick).astype(np.uint8), padding, mode="wrap")
    padded_mask_bool = np.pad(dgrid2 > mask_thick, padding, mode="wrap")

    reg_det = mgr_rust.seeded_watershed_det(
        padded_grid2.astype(np.float32), padded_markers2.astype(np.uint32), padded_mask_u8,
    )
    reg_sk2 = sk_watershed(-padded_grid2.astype(float), padded_markers2.astype(np.int32),
                           mask=padded_mask_bool, compactness=0, connectivity=2).astype(np.uint32)
    reg_sk3 = sk_watershed(-padded_grid2.astype(float), padded_markers2.astype(np.int32),
                           mask=padded_mask_bool, compactness=0, connectivity=3).astype(np.uint32)

    diff_sk2 = np.sum(reg_det != reg_sk2)
    diff_sk3 = np.sum(reg_det != reg_sk3)
    print(f"  mask={mask_thick:.2f}: rust_det vs skimage_c2 = {diff_sk2} diffs, vs skimage_c3 = {diff_sk3} diffs")

print("\nDone.")

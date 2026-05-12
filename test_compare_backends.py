#!/usr/bin/env python3
"""Compare watershed backends (rust_bucket, rust_det, skimage) on CuBTC.

Pars: spacing=0.5, h=0.2, mask_thickness=0.05, threshold_abs=1.5
Tests skimage with connectivity=2 and connectivity=3.
Reports number of connections found and identifies missing ones.
"""

import sys
import time
import numpy as np

sys.path.insert(0, "/Users/arungopalan/repos/mofography")

import mofography as mgr
import mgr_rust
from ase.io import read
from ase.data import covalent_radii

from skimage.segmentation import watershed as sk_watershed
from mofography.helper import find_touching_region_pairs

# ── Parameters ──────────────────────────────────────────────────
CIF_PATH = "/Users/arungopalan/repos/mofography/mofography/examples/TestCIFs/P1-Cu-BTC.cif"
SPACING = 0.5
H = 0.2
MASK_THICKNESS = 0.05
THRESHOLD_ABS = 1.5
MIN_DISTANCE = 1.5  # default
COMPACTNESS = 0

# ── Load & compute dgrid ────────────────────────────────────────
print("=" * 70)
print("Loading CuBTC ...")
mat_atoms = read(CIF_PATH)
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
radii = covalent_radii[mat_updated.get_atomic_numbers()]

t0 = time.time()
grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing=SPACING, block_size=20000)
dgrid = grid.compute().astype(np.float32)
print(f"  dgrid shape: {dgrid.shape}, time: {time.time() - t0:.1f}s")

# ── Detect maxima (same for all runs) ───────────────────────────
from ase.geometry import cellpar_to_cell

cell = cellpar_to_cell(mat_atoms.cell.cellpar()).T

merged, unmerged, raw = mgr_rust.detect_maxima_selective_smoothing_parallel(
    dgrid,
    preserve_threshold=3.0,
    smooth_threshold=6.0,
    max_smooth_radius=2.0,
    h=H,
    threshold_abs=THRESHOLD_ABS,
    min_distance=MIN_DISTANCE,
    cell=cell,
)
print(f"  Maxima detected: {len(merged)} merged, {len(unmerged)} unmerged")

if len(merged) == 0:
    print("ERROR: No maxima detected — aborting.")
    sys.exit(1)

# ── Build markers and mask ──────────────────────────────────────
markers = np.zeros_like(dgrid, dtype=np.int32)
coords = tuple(merged.T)
markers[coords] = np.arange(1, len(merged) + 1)

void_mask = (dgrid > MASK_THICKNESS).astype(np.uint8)
void_mask_bool = dgrid > MASK_THICKNESS

# Use a moderate padding for periodic segmentation
padding = 5
padded_grid = np.pad(dgrid, padding, mode="wrap")
padded_markers = np.pad(markers, padding, mode="wrap")
padded_mask_u8 = np.pad(void_mask, padding, mode="wrap")
padded_mask_bool = np.pad(void_mask_bool, padding, mode="wrap")

print(f"\n  padded_grid shape: {padded_grid.shape}")
print(f"  mask voxels (in padded): {padded_mask_bool.sum()} / {padded_mask_bool.size}")

# ── Helper: run a backend and return padded regions ─────────────
def run_backend(name, **kw):
    t0 = time.time()
    if name == "rust_bucket":
        reg = mgr_rust.seeded_watershed_bucket(
            padded_grid.astype(np.float32),
            padded_markers.astype(np.uint32),
            mask=padded_mask_u8.astype(np.uint8),
        )
    elif name == "rust_det":
        reg = mgr_rust.seeded_watershed_det(
            padded_grid.astype(np.float32),
            padded_markers.astype(np.uint32),
            padded_mask_u8.astype(np.uint8),
        )
    elif name.startswith("skimage_conn"):
        connectivity = kw.get("connectivity", 2)
        reg = sk_watershed(
            -padded_grid.astype(float),
            padded_markers.astype(np.int32),
            mask=padded_mask_bool,
            compactness=COMPACTNESS,
            connectivity=connectivity,
        ).astype(np.uint32)
    else:
        raise ValueError(f"Unknown backend: {name}")
    elapsed = time.time() - t0
    n_labels = len(np.unique(reg)) - 1
    return reg, n_labels, elapsed


# ── Run all backends ────────────────────────────────────────────
backends = [
    ("rust_det", "rust_det", {}),
    ("rust_bucket", "rust_bucket", {}),
    ("skimage_conn2", "skimage_conn2", {"connectivity": 2}),
    ("skimage_conn3", "skimage_conn3", {"connectivity": 3}),
]

results = {}
for label, name, kw in backends:
    reg, n_labels, elapsed = run_backend(name, **kw)
    results[label] = {"regions": reg, "n_labels": n_labels, "time": elapsed}
    print(f"  {label:20s} → {n_labels:4d} labels, {elapsed:.1f}s")

# ── Compute touching pairs for each ─────────────────────────────
print("\n" + "=" * 70)
print("Computing touching region pairs (6-NN boundary detection) ...")

all_pairs = {}
for label in results:
    reg_padded = results[label]["regions"]
    t0 = time.time()
    pairs = find_touching_region_pairs(reg_padded, connectivity=26)
    all_pairs[label] = pairs
    n_conn = len(pairs)
    results[label]["n_connections"] = n_conn
    print(f"  {label:20s} → {n_conn:4d} connections (boundary: 6-NN)")

# ── Compare: find connections in rust_det missing from skimage ──
print("\n" + "=" * 70)
print("COMPARISON: rus_det vs skimage_conn2\n")

# Convert to sets of frozensets for comparison
def pairs_to_set(pairs):
    return {frozenset(row) for row in pairs}

ref = pairs_to_set(all_pairs["rust_det"])

for target in ["skimage_conn2", "skimage_conn3"]:
    comp = pairs_to_set(all_pairs[target])
    missing_from_target = ref - comp
    extra_in_target = comp - ref

    print(f"  {'─' * 50}")
    print(f"  rust_det  connections: {len(ref)}")
    print(f"  {target:15s} connections: {len(comp)}")
    print(f"  Missing from {target}: {len(missing_from_target)}")
    print(f"  Extra in {target}:    {len(extra_in_target)}")

    if missing_from_target:
        print(f"\n  Missing connections (in rust_det not in {target}):")
        for pair in sorted(missing_from_target):
            a, b = tuple(sorted(pair))
            print(f"    regions {a} ↔ {b}")

    if extra_in_target:
        print(f"\n  Extra connections (in {target} not in rust_det):")
        for pair in sorted(extra_in_target):
            a, b = tuple(sorted(pair))
            print(f"    regions {a} ↔ {b}")

# ── Compare rust_det vs rust_bucket ─────────────────────────────
print("\n" + "=" * 70)
print("COMPARISON: rust_det vs rust_bucket\n")

bucket_set = pairs_to_set(all_pairs["rust_bucket"])
missing_in_bucket = ref - bucket_set
extra_in_bucket = bucket_set - ref

print(f"  rust_det    connections: {len(ref)}")
print(f"  rust_bucket connections: {len(bucket_set)}")
print(f"  Missing from rust_bucket: {len(missing_in_bucket)}")
print(f"  Extra in rust_bucket:    {len(extra_in_bucket)}")

if missing_in_bucket:
    print(f"\n  Missing connections (in rust_det not in bucket):")
    for pair in sorted(missing_in_bucket):
        a, b = tuple(sorted(pair))
        print(f"    regions {a} ↔ {b}")

# ── Summary ─────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY")
print(f"  Maxima: {len(merged)}")
for label in results:
    r = results[label]
    print(f"  {label:20s} → {r['n_labels']:4d} labels, {r['n_connections']:4d} connections ({r['time']:.1f}s)")

print("\nDone.")

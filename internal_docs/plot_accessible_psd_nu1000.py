"""
Plot occupiable and accessible PSDs for NU-1000 and open in browser.

Run from repo root: python docs/plot_accessible_psd_nu1000.py

Segmentation results are cached to NU-1000_S1_segmentation.npz in this directory.
If the cache exists, segmentation is skipped and only PSD plotting runs.
Delete the npz file to force recomputation.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np
from ase.io import read

# Use mofography from the repo
import mofography as mgr

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
CIF_PATH = REPO_ROOT / "examples" / "TestCIFs" / "NU-1000_S1.cif"
CACHE_PATH = SCRIPT_DIR / "NU-1000_S1_segmentation.npz"

# Segmentation params (from docs/test_percolated_cluster_coloring.py)
SPACING = 0.25
MASK_THICKNESS = 0.5
H = 0.2
EXTRA_PADDING = 5.0
THRESHOLD_ABS = 1.5

# PSD defaults (match interactive app)
NUM_SAMPLES_TOTAL = 10000
NUM_SAMPLES_PER_REGION = 1000
NUM_BINS = 200
MAX_DIAMETER = 50.0
BINS = np.linspace(0, MAX_DIAMETER, NUM_BINS + 1)

MATERIAL_NAME = "NU-1000"
PROBE_RADII = [3.25, 10.0]


def _serializable_value(v):
    """Return a value suitable for np.savez (array or object array for pickle)."""
    if isinstance(v, np.ndarray):
        return v
    if np.isscalar(v):
        return np.array(v)
    return np.array([pickle.dumps(v)], dtype=object)


def save_segmentation_results(results: dict, path: Path) -> None:
    """Save results dict to npz, excluding ase_atoms; pickle non-array values."""
    to_save = {}
    for k, v in results.items():
        if k == "ase_atoms":
            continue
        to_save[k] = _serializable_value(v)
    np.savez_compressed(path, **to_save)


def load_segmentation_results(path: Path) -> dict:
    """Load results dict from npz; unpickle object arrays."""
    data = np.load(path, allow_pickle=True)
    results = {}
    for k in data.files:
        v = data[k]
        if v.dtype == object and v.size == 1:
            try:
                results[k] = pickle.loads(v.item())
            except Exception:
                results[k] = v.item()
        else:
            results[k] = v
    return results


def get_or_compute_segmentation(cif_path: Path, cache_path: Path):
    """Load cached segmentation or compute and cache. Returns (results, mat_atoms)."""
    mat_atoms = read(str(cif_path))

    if cache_path.exists():
        print(f"Loading cached segmentation from {cache_path}")
        results = load_segmentation_results(cache_path)
        results["ase_atoms"] = mat_atoms
        if "dgrid" not in results:
            raise ValueError("Cached npz missing 'dgrid'; delete cache and recompute with return_dgrid=True.")
        return results, mat_atoms

    print("Computing dgrid (Rust)...")
    dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=SPACING)
    print("Computing segmentation...")
    results = mgr.do_periodic_segmentation_rust(
        dgrid,
        mat_atoms,
        spacing=SPACING,
        mask_thickness=MASK_THICKNESS,
        h=H,
        threshold_abs=THRESHOLD_ABS,
        extra_padding=EXTRA_PADDING,
        return_dgrid=True,
        compute_connections=True,
    )
    results["ase_atoms"] = mat_atoms
    if "dgrid" not in results:
        results["dgrid"] = dgrid

    print(f"Saving segmentation cache to {cache_path}")
    save_segmentation_results(results, cache_path)
    return results, mat_atoms


def main():
    print("CIF:", CIF_PATH)
    if not CIF_PATH.exists():
        raise FileNotFoundError(f"CIF not found: {CIF_PATH}")

    results, mat_atoms = get_or_compute_segmentation(CIF_PATH, CACHE_PATH)
    dgrid = results["dgrid"]
    regions_pbc = results["regions_pbc"]

    # --- Occupiable global PSD ---
    print("Opening occupiable global PSD in browser...")
    _, uniform_fig = mgr.compute_psd_uniform_rust(
        dgrid,
        mat_atoms,
        num_samples_total=NUM_SAMPLES_TOTAL,
        bins=BINS,
        density=False,
        return_fig=True,
        material=MATERIAL_NAME,
    )
    uniform_fig.show(renderer="browser")

    # --- Occupiable pocket-wise PSD ---
    print("Opening occupiable pocket-wise PSD in browser...")
    _, pocket_fig = mgr.compute_psds_pocketwise_rust(
        dgrid,
        regions_pbc,
        mat_atoms,
        num_samples_per_region=NUM_SAMPLES_PER_REGION,
        bins=BINS,
        return_fig=True,
        material=MATERIAL_NAME,
    )
    pocket_fig.show(renderer="browser")

    # --- Accessible global and pocket-wise at 3.25 and 10.0 ---
    for probe_radius in PROBE_RADII:
        print(f"Opening accessible global PSD (probe {probe_radius} Å) in browser...")
        _, acc_global_fig = mgr.compute_psd_global_accessible(
            results,
            mat_atoms,
            probe_radius,
            bins=BINS,
            samples_per_region=NUM_SAMPLES_PER_REGION,
            return_fig=True,
            material=MATERIAL_NAME,
        )
        acc_global_fig.show(renderer="browser")

        print(f"Opening accessible pocket-wise PSD (probe {probe_radius} Å) in browser...")
        _, acc_pocket_fig = mgr.compute_psds_pocketwise_accessible(
            results,
            mat_atoms,
            probe_radius,
            num_samples_per_region=NUM_SAMPLES_PER_REGION,
            bins=BINS,
            return_fig=True,
            material=MATERIAL_NAME,
        )
        acc_pocket_fig.show(renderer="browser")

    print("Done. Six figures should have opened in the browser.")


if __name__ == "__main__":
    main()

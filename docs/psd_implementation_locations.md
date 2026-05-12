# PSD implementation locations

Short reference for finding pore size distribution (PSD) code: occupiable vs accessible, Rust vs Python, and how to test.

## Rust backend

**File:** `mgr_rust/src/lib.rs`

| Function | Approx. line | Role |
|----------|--------------|------|
| `compute_pore_sizes_uniform` | ~5311 | Uniform sampling over all pore voxels (dgrid > 0); used for **global occupiable** PSD. |
| `compute_pore_sizes_for_psd_pocketwise` | ~5384 | Sample within each region; accepts `regions_pbc` (or masked). Used for **pocket-wise occupiable** and for **accessible** (Python passes masked regions). |

Pore-space filter in uniform: `if dgrid_val > 0.00` (no region/percolation filter in Rust).

## Python layer

**File:** `mofography/rust_fun.py`

| Function | Approx. line | Role |
|----------|--------------|------|
| `compute_pore_sizes_uniform_rust` | ~1358 | Wraps Rust uniform sampler; returns `(sample_points, pore_sizes)`. |
| `compute_psd_uniform_rust` | ~1384 | **Global occupiable PSD.** Calls the above, histograms, optional figure. Returns `psd_hist` or `(psd_hist, fig)`. |
| `compute_pore_sizes_pocketwise_rust` | ~1455 | Wraps Rust pocketwise sampler; returns `(sample_points, pore_sizes, region_labels)`. |
| `compute_psds_pocketwise_rust` | ~1490 | **Pocket-wise occupiable PSD.** Uses above + `compute_psds_from_points_pocketwise`; returns `pocket_psds` or `(pocket_psds, fig)`. |
| `compute_accessible_psd_from_results` | ~1544 | **Accessible PSD (core).** Uses `percolation_all_axes_winding` to get percolating regions at `probe_radius`, masks `regions_pbc`, then calls Rust pocketwise on masked regions and builds volume-weighted global + per-region. Returns dict: `pdf`, `cum`, `psd`, `bin_centers`, `per_region`, `accessible_region_ids`. |
| `compute_psd_global_accessible` | ~1721 | **Global accessible PSD.** Thin wrapper over `compute_accessible_psd_from_results`; returns same shape as `compute_psd_uniform_rust`: `psd_hist` or `(psd_hist, fig)`. |
| `compute_psds_pocketwise_accessible` | ~1813 | **Pocket-wise accessible PSD.** Wrapper over `compute_accessible_psd_from_results`; returns same dict shape as `compute_psds_pocketwise_rust`: `pocket_psds` or `(pocket_psds, fig)`. |

Accessible functions require a `results` dict from `do_periodic_segmentation_rust(..., return_dgrid=True)` (must include `dgrid` and segmentation/graph for percolation).

**Percolation:** `mofography/winding_percolation.py` — `percolation_all_axes_winding` (approx. line 542).

**Helper:** `mofography/helper.py` — `compute_psds_from_points_pocketwise` (approx. line 1075); builds per-pocket PSD dict from `(region_labels, pore_sizes, bins)`.

## Testing and validation

| What | Location |
|------|----------|
| Unit tests (import, return shapes, optional integration) | `mofography/tests/test_mofography.py` — e.g. `test_accessible_psd_functions_exported`, `test_compute_psd_global_accessible_return_format`, `test_compute_psds_pocketwise_accessible_return_format`. Run: `pytest mofography/tests/test_mofography.py -v`. |
| Standalone PSD validation script (NU-1000, occupiable + accessible at chosen probe radii, browser plots) | `docs/plot_accessible_psd_nu1000.py`. Run from repo root: `python docs/plot_accessible_psd_nu1000.py`. Uses cached segmentation in `docs/NU-1000_S1_segmentation.npz` if present. |

Quick test snippet (needs segmentation with `return_dgrid=True` and a CIF):

```python
import mofography as mgr
from ase.io import read

mat = read("mofography/examples/TestCIFs/P1-IRMOF-1.cif")
dgrid = mgr.dgrid_from_atoms_rust(mat, spacing=0.5)
results = mgr.do_periodic_segmentation_rust(dgrid, mat, spacing=0.5, return_dgrid=True, compute_connections=True)

# Occupiable
psd_hist = mgr.compute_psd_uniform_rust(dgrid, mat, num_samples_total=5000)
pocket_psds = mgr.compute_psds_pocketwise_rust(dgrid, results["regions_pbc"], mat, num_samples_per_region=500)

# Accessible at probe_radius
psd_hist_acc, fig = mgr.compute_psd_global_accessible(results, mat, probe_radius=1.0, return_fig=True)
pocket_acc = mgr.compute_psds_pocketwise_accessible(results, mat, probe_radius=1.0, num_samples_per_region=500)
```

## See also

- `docs/psd_detailed_notes.md` — Gelb–Gubbins method and implementation notes.
- `docs/pore_size_distributions_gelb_gubbins.md` — User-facing PSD workflow.
- `accessible_psd_analysis_prompt.md` (repo root) — Analysis that led to the accessible PSD API.

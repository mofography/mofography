---
name: mofography
description: >-
  Use the mofography Python module to segment the void space of crystalline nanoporous
  materials (MOFs, zeolites, polymers) into pockets, build pore graphs, and compute
  textural properties (LCD, PLD, LIS, accessible/total surface area and void volume, PSDs),
  pocket types, windows and percolation. Use when the user gives a CIF/ASE structure and
  asks for pore analysis in Python, or wants to write or debug code that imports mofography.
  Not for the `mgrcli` command-line tool or the Streamlit web app (separate repo).
---

# mofography (Python module)

Scope: the `mofography` package in this repo (Python API on a Rust core, `mgr_rust`).
It does **not** cover the `mgrcli` CLI or the web app (`mgr_cli_py`, a separate repo with its own README).

## Setup check

```bash
python -c "import mgr_rust, mofography as m; print(m.__version__)"
```

If the import fails, `mgr_rust` is not built. Build it in release mode (slow in debug):

```bash
python -m pip install maturin
cd mgr_rust && maturin develop --release && cd .. && python -m pip install -e .
```

After `git pull`, rebuild with `cd mgr_rust && maturin develop --release`.
Install scripts: `scripts/install_mofography.sh` (Linux/macOS), `scripts/install_mofography.ps1` (Windows).

## Core workflow

Input is an ASE `Atoms` object with a periodic cell. Everything starts from a distance grid.

```python
import numpy as np
import mofography as mgr
from ase.io import read

atoms = read("examples/TestCIFs/P1-IRMOF-1.cif")
spacing = 0.5  # Å; halving it gives 8x more voxels

# 1. distance grid (distance from each voxel to the nearest framework atom surface)
dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=spacing)

# 2. segmentation: maxima, watershed, periodic handling, windows, pore graph
seg = mgr.do_periodic_segmentation_rust(
    dgrid.astype(np.float32), atoms,
    mask_thickness=0.05, h=0.2, min_distance=1.5, threshold_abs=1.5,
    spacing=spacing,
    compute_connections=True,
    return_dgrid=True,            # keep dgrid in the result (needed by textural properties)
    disable_tqdm=True,
)
```

Useful keys of `seg`: `dgrid`, `regions`, `regions_pbc`, `regions_padded`, `all_maxima`, `maxima_pbc`, `pbc_groups`,
`unique_regions`, `connections` (windows between pockets), `close_or_inside_regions`, `maxima_info`.
`len(seg["maxima_pbc"])` is the pocket count; check it against the expected topology.

`mgr.quick_cluster(atoms)` is a one-call convenience that also builds figures. Prefer the explicit steps above in scripts.

## Textural properties

Pass the segmentation result `seg` (it carries the grid) plus `atoms`. Probe radius is in Å (1.2 ≈ helium-like point probe).

```python
# Diameters (Å): returns (LCD, PLD, LIS, dgrid). LIS = largest included sphere along percolating pathways.
lcd, pld, lis, _ = mgr.get_pore_dia_mgr_winding(results=seg, mat_atoms=None, include_dangling=True)

# Optional: winding percolation history (list, pld_locations, clusters); the paper pipeline stores it on seg
hist, pld_locations, clusters = mgr.compute_percolation_history_winding(seg, include_dangling=True)
seg["percolation_history"] = hist

# Surface area, Monte Carlo shell + winding percolation accessibility (m²/g)
np.random.seed(42)                         # MC is random; seed for reproducible numbers
sa = mgr.compute_surface_area_mc_winding(seg, atoms, probe_radius=1.2, n_samples=50000, shell_thickness=0.5)
sa["total_surface_area"], sa["accessible_surface_area"], sa["accessible_fraction"]
sa["by_pocket"], sa["by_region"]           # pocket- and region-level breakdowns

# Void volume, voxel counting + winding accessibility (cm³/g)
vol = mgr.compute_void_volume_from_results_winding(seg, atoms, probe_radius=1.2)
vol["total_void_volume"], vol["accessible_void_volume"], vol["by_pocket"]
```

Pore size distributions (mind radius vs diameter):

```python
# Global, uniform sampling. Returns (sample_points, RADII); diameter = 2 * radii.
# min_sphere_center_radius = probe clearance for the sphere center (1.657 Å = N2, matches PoreBlazer Total PSD).
_, radii = mgr.compute_pore_sizes_uniform_rust(seg["dgrid"], atoms, num_samples_total=50000,
                                               min_sphere_center_radius=1.657)
diameters = 2.0 * radii[radii > 1.657]

# Quick histogram only (counts over bins=linspace(0, 50, 200) by default; no probe clearance option)
counts = mgr.compute_psd_uniform_rust(seg["dgrid"], atoms, num_samples_total=10000)

# Pocket-wise (same sampling per pocket, so small pockets are not undersampled)
from mofography.pocket_typing import psd_bins
psds = mgr.compute_psds_pocketwise_rust(seg["dgrid"], seg["regions_pbc"], atoms,
                                        num_samples_per_region=1000,
                                        bins=psd_bins(seg["dgrid"]), probe_center_min_radius=1.2)
```

Meaning of the terms (keep them consistent in explanations):
- **Total** = occupiable by the probe center, in all segmented pore regions.
- **Accessible** = probe-center accessible, only regions that percolate for that probe radius.
- LCD = 2 x global max of the distance grid. PLD = largest probe diameter that still percolates. LIS = percolation-accessible analogue of LCD.
- Use **LIS** only. Do not use the older name DFSD. (The third return value of `get_pore_dia_mgr_winding` is LIS; it is called `dfsd` in older code.)
- Optional guardrails: `kwargs={"window_fp_for_lis": True}` on `get_pore_dia_mgr_winding` makes the LIS robust on small unit cells (window-fingerprint rescue of the largest cavity behind a narrow window); it is off by default. `window_fp_for_accessible` (on by default) is a separate switch for accessible surface area and volume; the two flags are independent. The old key `window_fp_for_dfsd` is still read but warns. `mgr.get_pore_dia_mgr_winding2(..., auto_rerun=True)` is the same but re-segments with `mask_thickness=1.5` when it detects deterministic-watershed spill (a pocket's self-loop radii spread over 30% of its radius); give it `mat_atoms`.
- Mesh surface area is also available: `mgr.compute_surface_area_mesh_winding(seg, atoms, probe_radius=...)`.
- Deprecated aliases you will see in older scripts: `compute_surface_area_from_results_mc` -> `compute_surface_area_mc_winding`, `compute_surface_area_from_results` -> `compute_surface_area_mesh_winding`. `compute_void_volume_from_results` is the underlying implementation that `compute_void_volume_from_results_winding` wraps. Prefer the `*_winding` names in new code.

Sanity check on `P1-IRMOF-1.cif` at spacing 0.5 (Bondi radii): LCD ≈ 15.1 Å, PLD ≈ 7.9 Å, 8 pockets, SA ≈ 3.9e3 m²/g, void volume ≈ 0.93 cm³/g. With CRC radii (below): SA ≈ 3.7e3 m²/g, void volume ≈ 0.94 cm³/g.

## Settings used for the paper (Zeo++/PoreBlazer comparison)

```python
from mofography.vdw_radii import mendeleev_vdw_radii_angstrom   # CRC handbook radii
radii = mendeleev_vdw_radii_angstrom(atoms)
dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=0.5, radii=radii).astype(np.float32)
seg = mgr.do_periodic_segmentation_rust(
    dgrid, mat_atoms=atoms, spacing=0.5, mask_thickness=0.05, h=0.2,
    threshold_abs=1.5, min_distance=1.5, extra_padding=5,
    watershed_backend="rust_det", return_dgrid=True, compute_connections=True,
    vdw_radii_angstrom=radii, vdw_radii_set="crc",
)
```

Probe 1.2 Å (He) for surface area and volume, 50,000 MC samples with shell thickness 0.5 Å, seed 42.
Run the unit cell directly: `dgrid_from_atoms_rust` handles non-orthogonal cells correctly, so no supercell is needed. Add `kwargs={"window_fp_for_lis": True}` to the LIS call for small cells so the LIS is not missed.

## Windows

`seg["connections"]` is a list; each entry is `[region_pair, valid, kind, window_positions, window_radii]`
(`kind` is `"periodic"` or `"internal"`). Window **radius** is the distance-grid value at the window voxel, so diameter = 2 x radius.

```python
window_diameters = []
for c in seg["connections"]:
    if not c[1]:                       # skip invalid entries
        continue
    window_diameters.extend(2.0 * float(r) for r in np.atleast_1d(c[-1]))
```

`mofography.helper.trim_connections_to_unique_regions(connections, seg["unique_regions"])` keeps only connections that touch at least one region inside the central cell (use it for unit-cell pore graphs and plots).
For IRMOF-1 the windows come out at 7.4 to 7.9 Å.

## Pocket types

Paper pipeline (`env_cluster` typing, default method):

```python
from mofography.pocket_typing import assign_pocket_types_from_segmentation, psd_bins

mgr.update_with_covalent_bonds(atoms)                       # bond info for the framework fingerprints
seg["pocket_psds"] = mgr.compute_psds_pocketwise_rust(
    seg["dgrid"], seg.get("regions_pbc", seg["regions"]), atoms,
    num_samples_per_region=1000, bins=psd_bins(seg["dgrid"]), probe_center_min_radius=1.2,
)
typed = assign_pocket_types_from_segmentation(seg, atoms, method="env_cluster", probe_radius_psd=1.2)
assert typed.ok, typed.error
typed.n_types, typed.n_pockets, typed.type_summary           # type_summary: list of dicts per type
```

Shortcut with the same result: `enriched = mgr.enrich_segmentation_with_pocket_psds(seg, atoms)` then
`typed = mgr.assign_pocket_types(enriched, atoms, method="env_cluster")` (a thin wrapper of the function above).
By default typing reads `max_r` bands from a second Bondi/CRC grid; pass `typing_radii="segmentation"` to use the segmentation grid instead.

## Visualization

```python
# Standalone interactive Three.js HTML (preferred backend)
mgr.show_mof_threejs(seg, atoms, typed_result=typed, framework_atoms=atoms, output_html="mof.html")

# Plotly figure for notebooks
fig = mgr.show_mof_plotly(seg, atoms, typed_result=typed, framework_atoms=atoms, stylize=True)
```

Three.js HTML is not stored in the repo; generate it on demand. Percolated-cluster scene (accessible pockets at the percolation probe, pore graph, PLD windows):

```python
from mofography.helper import trim_connections_to_unique_regions
from mofography.workflow import show_percolation_threejs

hist, pld_locations, _ = mgr.compute_percolation_history_winding(seg, step=0.05, expand_components=False)
probe_r = float(hist[-1][0])                       # probe radius at which percolation stops
analysis = {
    "percolation_history": hist, "pld_locations": pld_locations,
    "dgrid": seg["dgrid"], "regions": seg["regions"],
    "connections": trim_connections_to_unique_regions(seg["connections"], seg["unique_regions"]),
    "maxima": seg["close_or_inside_maxima"], "isovalue": 1.0,
}
show_percolation_threejs(analysis, atoms, probe_radius=probe_r, output_html="percolation.html")
```

Options: `layout="white"`, `camera_eye=(x, y, z)`, `camera_center=(x, y, z)` (Cartesian Å).

## Parameters to tune

| Parameter | Default | Effect |
|---|---|---|
| `spacing` | 0.5 Å | grid resolution; the main cost/accuracy knob |
| `mask_thickness` | 0.05 Å | void threshold; raise to merge nearby pockets |
| `h` | 0.2 Å | maxima persistence; raise if over-segmented, lower if pockets are missing |
| `min_distance` | 1.5 Å | minimum separation of maxima |
| `threshold_abs` | 1.5 Å | maxima must be at least this far from the framework |
| `watershed_backend` | `rust_det` | alternatives: `rust_bucket`, `rust_det_parallel`, `skimage` |

Large cells: `dgrid_from_atoms_rtree_rust` (less memory) or `dgrid_from_atoms_rtree_interp_rust` (coarse grid plus interpolation).

## Gotchas

- Pass a **float32** grid to `do_periodic_segmentation_rust` (`dgrid.astype(np.float32)`).
- Textural-property functions need `return_dgrid=True` on the segmentation call, or they cannot find `seg["dgrid"]`.
- `compute_pore_sizes_uniform_rust` returns radii; double them for diameters. `compute_psd_uniform_rust` returns histogram counts, not samples.
- Monte Carlo surface areas vary by about 1% between unseeded runs; seed NumPy when comparing.
- Log lines such as `[pocket_maxima_check] ... underestimates true region max` are diagnostics, not errors.
- Use the same `probe_radius` consistently when comparing surface area, volume and accessibility.
- Do not add the PBC/percolation filter by hand; the `*_winding` functions already apply it.

## Where to look

- Tutorials: `docs/` (`quick_start.qmd`, `compute_pore_diameters.ipynb`, `pld_and_percolation.ipynb`, `pockets_detection.qmd`, `connections_and_windows.qmd`, `hts_using_mofography.ipynb`)
- Surface area and accessibility theory: `docs/surface_area.qmd`, `docs/accessibility.qmd`
- Rust backend notes: `docs/rust_backend_documentation.md`
- Test structures: `examples/TestCIFs/`

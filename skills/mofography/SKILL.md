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
lcd, pld, lis, _ = mgr.get_pore_dia_mgr_winding(mat_atoms=atoms, results=seg)

# Surface area, Monte Carlo shell + winding percolation accessibility (m²/g)
sa = mgr.compute_surface_area_mc_winding(seg, atoms, probe_radius=1.2)
sa["total_surface_area"], sa["accessible_surface_area"], sa["accessible_fraction"]
sa["by_pocket"], sa["by_region"]          # pocket- and region-level breakdowns

# Void volume, voxel counting + winding accessibility (cm³/g)
vol = mgr.compute_void_volume_from_results_winding(seg, atoms, probe_radius=1.2)
vol["total_void_volume"], vol["accessible_void_volume"], vol["by_pocket"]

# Pore size distributions
psd = mgr.compute_psd_uniform_rust(seg["dgrid"], atoms, num_samples_total=10000)               # global
psds = mgr.compute_psds_pocketwise_rust(seg["dgrid"], seg["regions_pbc"], atoms,
                                        num_samples_per_region=1000)                          # per pocket
```

Meaning of the terms (keep them consistent in explanations):
- **Total** = occupiable by the probe center, in all segmented pore regions.
- **Accessible** = probe-center accessible, only regions that percolate for that probe radius.
- LCD = 2 x global max of the distance grid. PLD = largest probe diameter that still percolates. LIS = percolation-accessible analogue of LCD.
- Use **LIS** only. Do not use the older name DFSD.
- Mesh surface area is also available: `mgr.compute_surface_area_mesh_winding(seg, atoms, probe_radius=...)`.

Sanity check on `P1-IRMOF-1.cif` at spacing 0.5: LCD ≈ 15.1 Å, PLD ≈ 7.9 Å, 8 pockets, SA ≈ 3.9e3 m²/g, void volume ≈ 0.93 cm³/g.

## Pocket types

```python
enriched = mgr.enrich_segmentation_with_pocket_psds(seg, atoms)
typed = mgr.assign_pocket_types(enriched, atoms, method="env_cluster")
typed.n_types, typed.n_pockets, typed.type_summary   # type_summary: list of dicts per type
```

## Visualization

```python
# Standalone interactive Three.js HTML (preferred backend)
mgr.show_mof_threejs(seg, atoms, typed_result=typed, framework_atoms=atoms, output_html="mof.html")

# Plotly figure for notebooks
fig = mgr.show_mof_plotly(seg, atoms, typed_result=typed, framework_atoms=atoms, stylize=True)
```

Three.js HTML is not stored in the repo; generate it on demand. Percolation scenes: `show_percolation_threejs` (see `docs/pld_and_percolation.ipynb`).

## Parameters to tune

| Parameter | Default | Effect |
|---|---|---|
| `spacing` | 0.5 Å | grid resolution; the main cost/accuracy knob |
| `mask_thickness` | 0.05 Å | void threshold; raise to merge nearby pockets |
| `h` | 0.2 Å | maxima persistence; raise if over-segmented, lower if pockets are missing |
| `min_distance` | 1.5 Å | minimum separation of maxima |
| `threshold_abs` | 1.5 Å | maxima must be at least this far from the framework |
| `watershed_backend` | `rust_det` | alternatives: `rust_bucket`, `rust_det_parallel`, `skimage` |

Large cells: `dgrid_from_atoms_aabb_rust` (less memory) or `dgrid_from_atoms_aabb_interp_rust` (coarse grid plus interpolation).

## Gotchas

- Pass a **float32** grid to `do_periodic_segmentation_rust` (`dgrid.astype(np.float32)`).
- Textural-property functions need `return_dgrid=True` on the segmentation call, or they cannot find `seg["dgrid"]`.
- Log lines such as `[pocket_maxima_check] ... underestimates true region max` are diagnostics, not errors.
- Use the same `probe_radius` consistently when comparing surface area, volume and accessibility.
- Do not add the PBC/percolation filter by hand; the `*_winding` functions already apply it.

## Where to look

- Tutorials: `docs/` (`quick_cluster.qmd`, `compute_pore_diameters.ipynb`, `pld_and_percolation.ipynb`, `pockets_detection.qmd`, `connections_and_windows.qmd`, `hts_using_mofography.ipynb`)
- Surface area and accessibility theory: `docs/surface_area.qmd`, `docs/accessibility.qmd`
- Rust backend notes: `docs/rust_backend_documentation.md`
- Test structures: `examples/TestCIFs/`

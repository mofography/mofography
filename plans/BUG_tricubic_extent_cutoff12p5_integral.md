# BUG: Tricubic integral screening fails at `cutoff=12.5` (extent tiling incomplete)

**Date:** 2026-06-23  
**Status:** FIXED (`8a2df8a` — derivative scaling `extent/(N-1)` in Rust; remove incorrect extent folding in query)  
**Severity:** High for production screening if `cutoff=12.5` is used  
**Component:** mofography `TricubicField` + `egrid_query_extent_from_atoms`; mgr_sn integral v4  

---

## Summary

Increasing LJ `cutoff` from **6.0 Å** to **12.5 Å** should **not** materially change guest–framework energies in pore regions (LJ is negligible beyond ~8–10 Å). In practice, **xylene integral screening** on DONNAW_clean works at `cutoff=6.0` but returns **`n_kept=0`**, **`log_Z=-inf`**, **`winner_integral=None`** at `cutoff=12.5`.

This is **not a physics/chemistry issue** — it is a **grid build vs tricubic query coordinate** failure when `extent_a = 1/n_cells_a < 1` on multiple axes.

---

## Expected vs actual

### Expected

- Pocket integral on DONNAW (node 1, pocket type 0) at 298 K:
  - `winner_integral = 'pX'`
  - `margin_log_Z ≈ 18.5`
  - `log_Z` roughly: oX ~38, mX ~44, pX ~62
  - `n_poses_kept` ≫ 0 per isomer
- Energies at `cutoff=12.5` should match `cutoff=6.0` within small numerical tolerance in pore regions.

### Actual (`cutoff=12.5`)

```
winner_integral = None
margin_log_Z = None
oX/mX/pX: log_Z = -inf, n_kept = 0
```

### Actual (`cutoff=6.0`) — correct

```
winner_integral = pX
margin_log_Z ≈ 18.51
oX: log_Z ≈ 38.32, n_kept = 26
mX: log_Z ≈ 43.52, n_kept = 12
pX: log_Z ≈ 62.03, n_kept = 51
```

---

## Why cutoff changes more than “2×2×2”

`get_repetitions_to_meet_cutoff` (`mofography/compute.py`):

```python
nx_cells = int(np.ceil(2 * cutoff / lx_unit))
# same for y, z using orthogonal unit lengths lx_unit, ly_unit, lz_unit
```

**DONNAW_clean** (cell ≈ 17.5 × 6.7 × 11.9 Å):

| `cutoff` (Å) | `n_cells` `(nx, ny, nz)` | `extent` `(1/nx, 1/ny, 1/nz)` | Supercell images |
|--------------|--------------------------|-------------------------------|------------------|
| 6.0          | `(1, 2, 2)`              | `(1.0, 0.5, 0.5)`             | 4                |
| 12.5         | `(2, 4, 3)`              | `(0.5, 0.25, 0.333…)`         | 24               |

Grid array shape stays **unit-cell sized**; only `extent_*` shrinks per axis. More axes tile at 12.5 (including **x**), exposing query bugs.

---

## Root cause (hypothesis)

1. **Grid build** (`egrid_tricubic_from_atoms_rust`, `_egrid_unit_cell_grid_inputs`):
   - Supercell atoms supply MIC images.
   - Stored egrid covers **first tile only**: ASE span `extent_a` per axis.

2. **Query** (`TricubicField.interpolate` / `interpolate_batch`):
   - Must map ASE fractional coords → tile-local query `q ∈ [0,1)³`:
     ```text
     f_tile = f_ASE mod extent_a
     q_a    = f_tile / extent_a
     ```
   - Partial fix in `5fe820c`: `extent_x/y/z` on `TricubicField`, folding in `_interpolate_impl`.

3. **What still fails:**
   - **Grid-node parity** at NU-1000 with `cutoff=12.5`, `extent=(1,1,0.5)`: scalar egrid vs tricubic interpolator → **R²=1, RMSE=0** (`scripts/dev/parity_egrid_nu1000.py`).
   - **Arbitrary pocket poses** on DONNAW at `extent=(0.5, 0.25, 0.33)`: tricubic energies catastrophic → all poses rejected.
   - Suggests folding works at grid nodes but **not** for general ASE coords on multi-axis small extent, or a separate bug (e.g. build/query convention mismatch on non-orthogonal cells, derivative scaling, edge of tile).

4. **Not the cause:**
   - Accessibility (winding @ KD/2) — identical for both cutoffs.
   - `dgrid` geom mask — uses precomputed `results["dgrid"]` RGI on full cell (no extent); unchanged with cutoff.
   - Wrong pocket / material aggregation — same node 1, pocket type 0.

---

## How to reproduce

### Environment

```bash
conda activate mgr_env
pip install -e ~/repos/mofography   # commit >= 5fe820c (extent folding)
cd ~/repos/mgr_sn/mgr_org_scripts
```

**Pickles** (must exist locally):

- `gcmc_candidates_dict.pkl`
- `core2000_xylene_rags_env_cluster.pkl`
- `core2000_mgr_results_crc_radii_05spacing_005mask_02h_optim_dgrid_unitcells_with_perc_10June2026.pkl`

### Repro A — single pocket (user’s working cell)

```python
import pickle
from pathlib import Path
from xylene_integral_screening import (
    load_gcmc_candidates_pickle,
    rag_index_for_material,
    build_material_fields,
    dgrid_rgi_from_results,
    integral_xylene_isomers_at_node,
)

SCRIPT_DIR = Path("~/repos/mgr_sn/mgr_org_scripts").expanduser()
gcmc = load_gcmc_candidates_pickle(SCRIPT_DIR)
rags = pickle.load(open(SCRIPT_DIR / "core2000_xylene_rags_env_cluster.pkl", "rb"))
results_list = pickle.load(open(
    SCRIPT_DIR / "core2000_mgr_results_crc_radii_05spacing_005mask_02h_optim_dgrid_unitcells_with_perc_10June2026.pkl",
    "rb",
))

row = gcmc["pX"]["pockets"].query("material_name == 'DONNAW_clean'").iloc[0]
ridx = int(row["rag_index"])
rag = rags[ridx]
results = results_list[ridx]
node_id = int(row["region_label"])  # = 1 for DONNAW

for cutoff in (6.0, 12.5):
  fields, _ = build_material_fields(rag.graph["ase_atoms"], cutoff=cutoff)
  print("cutoff", cutoff, "extent", fields.extent)
  dgrid_rgi = dgrid_rgi_from_results(results)
  pocket = integral_xylene_isomers_at_node(
      rag, results, node_id, fields, dgrid_rgi, temperature_k=298.0,
  )
  print(pocket["log_Z"], pocket["winner"], pocket["margin_log_Z"])
  print({iso: pocket["per_isomer"][iso]["n_poses_kept"] for iso in ("oX", "mX", "pX")})
  print()
```

**Pass:** `cutoff=6.0` → `winner='pX'`, finite log_Z, `n_kept > 0`.  
**Fail:** `cutoff=12.5` → `winner=None`, `log_Z=-inf`, `n_kept=0`.

### Repro B — material screen

```python
from xylene_integral_screening import screen_material_xylene_integral

ridx = rag_index_for_material(gcmc, "DONNAW_clean")
for cutoff in (6.0, 12.5):
    one = screen_material_xylene_integral(
        rags[ridx], results_list[ridx], rag_index=ridx,
        temperature_k=298.0, cutoff=cutoff, n_rot=400, max_trans=1000,
    )
    m = one["material"]
    p = one["pockets"][0]
    print(f"cutoff={cutoff} winner={m['winner_integral']} margin={m.get('margin_log_Z')}")
    for iso in ("oX", "mX", "pX"):
        pi = p["per_isomer"][iso]
        print(f"  {iso} log_Z={pi['log_Z']} n_kept={pi['n_poses_kept']}")
```

### Repro C — extent metadata only

```python
from mofography.rust_fun import egrid_query_extent_from_atoms
from mofography.compute import get_repetitions_to_meet_cutoff

atoms = rag.graph["ase_atoms"]
for cutoff in (6.0, 12.5):
    print(cutoff, get_repetitions_to_meet_cutoff(atoms, cutoff=cutoff),
          egrid_query_extent_from_atoms(atoms, cutoff=cutoff))
```

### Repro D — mofography grid-node parity (should PASS even at 12.5)

```bash
python ~/repos/mofography/scripts/dev/parity_egrid_nu1000.py
# NU-1000_S1, cutoff=12.5, |E|<100 K → R²=1
```

### Automated smoke (mgr_org_scripts)

```bash
python test_screen_accessibility_parity.py   # accessibility OK (not this bug)
# Add/f extend: test_integral_cutoff_parity.py comparing 6.0 vs 12.5 on DONNAW
```

---

## Files involved

| Layer | File | Role |
|-------|------|------|
| Rust | `mgr_rust/src/tricubic.rs` | `TricubicField`, `ase_frac_to_query`, `_interpolate_impl` |
| Rust | `mgr_rust/src/egrid_tricubic.rs` | Grid build with `extent_x/y/z` |
| Rust | `mgr_rust/src/lib.rs` | Scalar `egrid_from_atoms` (same extent convention) |
| Python | `mofography/rust_fun.py` | `_egrid_unit_cell_grid_inputs`, `egrid_query_extent_from_atoms`, `tricubic_interpolator_from_grid` |
| Python | `mofography/helper.py` | `ase_frac_to_egrid_query` |
| Python | `mofography/compute.py` | `get_repetitions_to_meet_cutoff`, `repeat_to_meet_cutoff` |
| Local | `mgr_org_scripts/xylene_energy_shape_match_v4.py` | `TricubicProbeFields.build_from_atoms` passes extent |
| Local | `mgr_org_scripts/xylene_integral_screening.py` | `build_material_fields(cutoff=6.0 default)`, `screen_material_xylene_integral` |

---

## Suggested fix directions

1. **Verify** scalar `egrid_from_atoms_rust[i,j,k]` vs `TricubicField` at **pocket pose** fractional coords (not just grid nodes) for DONNAW @ 12.5 — isolate build vs query.

2. **Check non-orthogonal cells:** DONNAW is orthogonal; also test NU-1000 (γ=120°) integral poses at `cutoff=12.5`.

3. **Audit** whether `ase_frac_to_query` matches Python `ase_frac_to_egrid_query` for all `f_ASE ∈ [0,1)` and edge cases (`f=extent`, negative frac).

4. **Consider** building/querying in **centered** fractional coords (`-0.5 + …`) consistently end-to-end instead of ASE `[0,1)` + fold.

5. **Regression tests:**
   - DONNAW `integral_xylene_isomers_at_node` @ 6.0 and 12.5: winners agree, log_Z within tolerance.
   - Extent folding unit test with `extent=(0.5, 0.25, 0.333)` on synthetic linear field.

---

## Workaround (production)

Use **`cutoff=6.0`** in integral screening until fixed:

```python
INTEGRAL_SCREEN_KW = dict(
    temperature_k=298.0,
    n_rot=400,
    max_trans=1000,
    seed=42,
    spacing=0.5,
    cutoff=6.0,
)
```

`build_material_fields` default is already `6.0`. Do **not** assume mofography’s `egrid_from_atoms_rust` default `12.5` applies to tricubic pocket integrals.

---

## Related docs

- `plans/BUG_tricubic_grid_issues.md` §4 (extent tiling — marked fixed prematurely)
- `docs/tricubic_energy_interpolation.org` (NU-1000 grid-node parity example)
- `mgr_org_scripts/xylene_integral_screening_REVIEW_FINDINGS.md`

---

## Acceptance criteria

- [ ] DONNAW_clean @ 298 K, node 1: `cutoff=12.5` gives `winner_integral='pX'`, `n_kept > 0`, log_Z within ~1–2 K of `cutoff=6.0` (not exact — extra images may shift sub-kJ amounts).
- [ ] Same for at least one non-orthogonal material (e.g. NU-1000_S1) if in validation set.
- [ ] `cargo test tricubic` + new integration test in `scripts/dev/`.
- [ ] No regression at `cutoff=6.0` on DONNAW / validation quartet.

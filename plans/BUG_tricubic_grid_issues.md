# Tricubic grid issues

## 1. Missing unit-cell slice — returns full supercell

**File:** `mofography/rust_fun.py:1989-1994`

`egrid_tricubic_from_atoms_rust` computes `nx, ny, nz` from **supercell** lengths:

```python
cell = cellpar_to_cell(mat_atoms_repeat.cell.cellpar()).T   # supercell cell
cell_lengths = np.linalg.norm(cell, axis=0)                  # supercell lengths
nx = int(cell_lengths[0] / spacing) + 1                      # e.g. 61 for 3×10Å/0.5Å
egrid8_flat = mgr_rust.egrid_tricubic_from_atoms(frac, cell, eps, sig, nx, ny, nz, cutoff)
egrid8 = egrid8_flat.reshape(8, nx, ny, nz)
return egrid8    # ← returns full supercell, no slice
```

Compare `egrid_from_atoms_rust` (line 182-189) which has the same pattern but **does** slice:

```python
egrid = egrid_super[:nx_unit, :ny_unit, :nz_unit]  # ← this is missing in tricubic
```

**Impact:** `get_energy_interpolator_tricubic` builds `RegularGridInterpolator` over `[0,1]` on whatever shape it receives. With a supercell-sized grid, querying fractional `0.5` maps to the middle of the supercell (unit cell 1 of NX), not the middle of the first unit cell. Wrong energies for any coordinate.

**Fix:** slice `egrid8[:, :nx_unit, :ny_unit, :nz_unit]` before returning (same as `egrid_from_atoms_rust`), OR compute `nx` from unit-cell lengths and grid only the unit cell.

**Status (Jun 2026):** fixed — unit-cell grid with `extent = 1/n_cells` per axis when replicated; supercell atom frac (no wrap); `extent=1` + unit-cell atoms when `n_cells=1`.

### 1b. Duplicate atom inflation (latent in first unit-cell fix)

Wrapping supercell Cartesian positions into principal-cell fractional coords collapses all periodic images onto the same site → Rust sums LJ `n_cells³` times (e.g. 6× for rep `(3,2,1)`). Fixed by keeping supercell `get_scaled_positions()-0.5` when replicated.

## 2. Grid point formula `-0.5 + i/(N-1)` vs `find_coord_from_indices` `i/N`

**File:** `mgr_rust/src/lib.rs` (Rust kernel) and `mofography/helper.py:33`

The Rust tricubic kernel uses grid points at `frac = -0.5 + i/(N-1)`. But `find_coord_from_indices` (used everywhere for index→coordinate conversion) uses `frac = i/N`:

```python
# helper.py:33
return np.dot(cell, (indices/shape).T).T   # frac = i / N
```

For the tricubic grid with `N = 61`: max positional error = `60/(61×60) = 1/61 ≈ 0.016` fractional units = ~0.5 Å for a 30 Å supercell. Error scales as `i/(N·(N−1))`, largest at far face.

**Fix:** either:
- Change `find_coord_from_indices` to `indices/(shape-1)` (affects everything, risky), or
- Change the Rust kernel to use `i/N` instead of `i/(N-1)` (linspace-like, but then -0.5..0.5 range needs adjustment), or
- Accept the small systematic error (diminishes with finer grids — at N=100, max ~0.15 Å)

## 3. (Nice-to-have) RASPA scaling uses `ncells = N−1` — internally consistent

**File:** `mofography/rust_fun.py:1998-2001`

```python
ncells = np.array([nx - 1, ny - 1, nz - 1], dtype=np.float64)
egrid8[1] /= ncells[0]   # etc.
```

This uses `N−1` which matches the kernel's denominator. If the kernel switches to `i/N`, this scaling should switch to `N` accordingly.

## 4. Extent tiling vs ASE query coords — **fixed Jun 2026**

**Symptom:** When `get_repetitions_to_meet_cutoff` > 1 on an axis, the egrid covers only the first unit-cell tile (`extent = 1/n_cells`), but `TricubicField.interpolate` treated ASE `[0,1)³` as full-cell query coords. ~92% of DONNAW pocket translations with `y > 0.5` queried the wrong tile → catastrophic energies and flat `log_Z`.

**Fix:**
- `mgr_rust.TricubicField`: `extent_x/y/z` params; fold ASE frac → tile query in `_interpolate_impl`.
- `tricubic_interpolator_from_grid(..., extent_*)`, `egrid_query_extent_from_atoms`, `ase_frac_to_egrid_query`.
- `get_energy_interpolator` / `get_energy_interpolator_tricubic` accept optional extent.
- Local `TricubicProbeFields.build_from_atoms` passes extent from cutoff replication.

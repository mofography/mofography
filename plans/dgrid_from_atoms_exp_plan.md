# Plan: `dgrid_from_atoms_exp` — Optimized Exact PBC Distance Grid

## Goal

Create a new experimental Rust function `dgrid_from_atoms_exp` in `mgr_rust/src/lib.rs` with Python wrapper `dgrid_from_atoms_exp_rust` in `mofography/rust_fun.py`. The function must be:
1. **Exact** for all unit cell types (orthogonal through triclinic)
2. **Fast** — restore near-original speed for orthogonal cells, minimize overhead for non-orthogonal
3. **Benchmarked** against the existing `dgrid_from_atoms` (current 27-image fix) and `dgrid_from_atoms_aabb` (R-tree ground truth)

## Background: The MIC Bug and Theory

### What was broken

The original `dgrid_from_atoms` used per-component fractional wrapping for minimum image convention (MIC):

```rust
dx -= dx.round();   // per-component wrap into [-0.5, 0.5]
dy -= dy.round();
dz -= dz.round();
```

This is the **Wassenaar/Tuckerman fractional-NINT algorithm** — convert displacement to fractional coordinates, round each component to the nearest integer (NINT), subtract to wrap into [-0.5, 0.5], then convert back to Cartesian.

**This is the same algorithm used by LAMMPS (`Domain::minimum_image`), GROMACS (`pbc_dx`), and OpenMM.**

### When per-component wrapping fails

The fractional cube [-0.5, 0.5]³ is the **Voronoi cell of the origin** in fractional space. But the true Voronoi cell in Cartesian space is defined by the metric tensor G = H^T H, where H is the cell matrix. For non-orthogonal cells, the Cartesian Voronoi cell is a skewed polytope that does NOT coincide with the fractional cube.

The per-component wrapping minimizes ||Δf||_∞ (L∞ norm in fractional space) but we need to minimize ||H·Δf||₂ (L2 norm in Cartesian space). These are equivalent **only when H is diagonal** (orthogonal cells).

**Key result (Wassenaar 2006, PhD thesis, University of Groningen, Appendix C):**
The fractional-NINT wrapping IS exact for all displacements satisfying:

    ||dr_wrapped|| < d_max = ½ × min(||a||, ||b||, ||c||)

where ||a||, ||b||, ||c|| are the lattice vector lengths. For displacements exceeding d_max, up to 6 additional image checks (corresponding to the faces of the Voronoi cell) are needed.

**Wassenaar's full algorithm:**
1. Apply fractional-NINT wrapping (Algorithm 1)
2. If ||dr_wrapped|| < d_max → done, result is exact
3. Otherwise, check 6 additional face-adjacent images: ±v₁, ±v₂, ±v₃ of the reduced basis
4. Total: 1 + 6 = 7 checks maximum (not 27)

The 6-image fallback requires the basis to be **reduced** (Niggli/Minkowski reduction) so that |μ_ij| ≤ 0.5 for all Gram-Schmidt coefficients. Without reduction, up to 14 images may be needed.

### Which Bravais lattices are affected

| Lattice type | Cell angles | MIC bug possible? |
|---|---|---|
| Cubic (P/I/F) | 90°, 90°, 90° | **No** — H is diagonal |
| Tetragonal (P/I) | 90°, 90°, 90° | **No** |
| Orthorhombic (P/I/F/C) | 90°, 90°, 90° | **No** |
| Hexagonal | 90°, 90°, 120° | **Yes** — off-diagonal H terms |
| Rhombohedral (R) | α≠90°, all equal | **Yes** |
| Monoclinic (P/C) | 90°, β≠90°, 90° | **Yes** |
| Triclinic (P) | α≠β≠γ≠90° | **Yes** |

**5 out of 14 Bravais types** are affected. In MOF chemistry, triclinic and monoclinic are very common.

### The FOTNIN failure case (confirmed empirically)

FOTNIN has a rhombohedral cell (all angles 60°). Shortest lattice vector ≈ 15.4 Å, so d_max ≈ 7.7 Å. The pore center sits at distance ~8 Å from the nearest atom, just barely exceeding d_max. The fractional wrapping picks the wrong image (Cartesian dist 28.4 Å instead of 15.9 Å), inflating LCD from 33.6 to 37 Å.

### Current `dgrid_from_atoms_aabb` correctness

The AABB method is our **ground truth** but has a subtle caveat:
- It replicates all N atoms into 27 Cartesian images and builds an R-tree
- For each grid point, it queries the `n_neighbors=10` nearest R-tree entries
- It finds the 10 nearest atom **centers**, not the 10 nearest atom **surfaces**
- The true minimum surface distance could come from a farther center with a larger radius
- Max VDW radius spread in our table is 1.17 Å (Na=2.27 - H=1.10)
- With 10 neighbors covering ~3-5 Å of distance range, this is **practically safe** but **not mathematically guaranteed** for pathological cases
- To be strictly exact, we would need to query enough neighbors to cover the full radius spread, or use a different distance metric in the R-tree

**Verdict:** AABB is exact for all practical MOF structures. We use it as the benchmark reference.

## Implementation Plan

### Step 1: Create `dgrid_from_atoms_exp` in `mgr_rust/src/lib.rs`

New PyO3 function with the following structure:

```
#[pyfunction]
fn dgrid_from_atoms_exp(
    py: Python,
    fpoints: PyReadonlyArray2<f64>,   // atom fractional positions (shifted by -0.5)
    cell: PyReadonlyArray2<f64>,       // cell matrix (3×3)
    radii: PyReadonlyArray1<f64>,      // atom radii
    nx: usize, ny: usize, nz: usize,   // grid dimensions
) -> PyResult<Py<PyArray3<f32>>>
```

### Step 2: Three-tier MIC strategy

**Tier 0: Orthogonal cell fast path (zero overhead)**

At function entry, check if the cell matrix is diagonal:
```rust
let is_orthogonal = cell[[0,1]].abs() < 1e-12
    && cell[[0,2]].abs() < 1e-12
    && cell[[1,2]].abs() < 1e-12;
```

If orthogonal: use the original single-image `dx -= dx.round()` path. This is mathematically exact for orthogonal cells and has zero overhead. Restores original speed for cubic, tetragonal, and orthorhombic cells.

**Tier 1: Fractional-NINT fast path with d_max guard (Wassenaar Algorithm 1)**

For non-orthogonal cells:
1. Precompute d_max = 0.5 × min(||a||, ||b||, ||c||) once at setup
2. For each (grid_point, atom) pair:
   a. Compute fractional difference: `df = grid_frac - atom_frac`
   b. Wrap: `df_wrapped = df - df.round()`
   c. Convert to Cartesian: `dr = cell × df_wrapped`
   d. Compute squared Cartesian distance: `sq = dr·dr`
   e. If `sq < d_max²` → this IS the true minimum image, use it
   f. If `sq ≥ d_max²` → fall through to Tier 2

**Expected hit rate:** >99% of atom-gridpoint pairs will hit the fast path, since most atoms are close to the grid point relative to d_max.

**Tier 2: 27-image brute force fallback**

Only for the rare pairs that fail the d_max check. Uses the same 27-offset loop as the current fix, with squared-distance early exit.

### Step 3: Additional optimizations

**3a. Precompute cell inverse for fractional-NINT (avoid per-component wrapping on fractional diffs)**

Wait — our input is already in fractional coordinates, so the wrapping IS per-component `df -= df.round()`. The conversion to Cartesian is `dr = cell × df`. So we're already doing the right thing. No inverse needed.

**3b. Squared-distance early exit in the 27-image fallback**

Already implemented in the current fix. Skip `sqrt` when `sq_dist >= (min_dist + r_atom)²`.

**3c. Future: SIMD batching with `wide` crate (NOT in this PR)**

Process 4 atoms at once using `f64x4` SIMD. Requires SoA layout for atom positions. This is a separate optimization that can be layered on later. The `wide` crate works on stable Rust and provides `f64x4::round()`, `f64x4::mul_add()`, `f64x4::min()`.

### Step 4: Create Python wrapper in `mofography/rust_fun.py`

```python
def dgrid_from_atoms_exp_rust(ase_atoms, spacing, radii=None):
    """Experimental optimized distance grid with exact PBC handling.
    
    Uses three-tier MIC strategy:
    - Orthogonal cells: single-image fast path (exact)
    - Non-orthogonal: fractional-NINT with d_max guard (Wassenaar 2006)
    - Fallback: 27-image brute force for rare corner cases
    """
    import mgr_rust
    from ase.data import vdw_radii as _ase_vdw
    
    if radii is None:
        radii = _ase_vdw[ase_atoms.get_atomic_numbers()]
    
    fpoints = ase_atoms.get_scaled_positions() - 0.5
    cell = _cell_matrix_for_rust(ase_atoms)
    nx, ny, nz = (ase_atoms.cell.cellpar()[:3] / spacing).astype(int) + 1
    return mgr_rust.dgrid_from_atoms_exp(
        fpoints, cell, int(nx), int(ny), int(nz), radii
    )
```

### Step 5: Benchmark suite

Test the following materials at spacing=0.5 and 0.25:

| Material | Cell type | Expected behavior |
|---|---|---|
| P1-Cu-BTC | Cubic | Tier 0 fast path, matches original speed |
| TOVKOG_clean | Triclinic | Tier 1 fast path for most pairs |
| FOTNIN_clean | Rhombohedral (60°) | Some pairs hit Tier 2 fallback |
| WUPNIF_clean | Triclinic | Tier 1 fast path |
| NU-1000 | Hexagonal | Tier 1 fast path |

For each:
1. Compare LCD against AABB (ground truth)
2. Time against: (a) original buggy `dgrid_from_atoms`, (b) current 27-image fix, (c) new `dgrid_from_atoms_exp`

## Files to modify

1. `mgr_rust/src/lib.rs` — add `dgrid_from_atoms_exp` function
2. `mgr_rust/src/lib.rs` — register in `#[pymodule]` 
3. `mofography/rust_fun.py` — add `dgrid_from_atoms_exp_rust` wrapper

## What NOT to do

- Do NOT modify `dgrid_from_atoms` (current 27-image fix) or `dgrid_from_atoms_aabb` — these remain as-is for comparison
- Do NOT add SIMD (`wide` crate) in this PR — that's a separate optimization
- Do NOT modify any Python `compute.py` functions — those are separate fixes
- Do NOT change the `get_pore_dia_mgr_winding` or `winding_percolation.py` functions

## References

1. Wassenaar, T.A. (2006). "Molecular Dynamics of Sense and Sensibility in Signaling." PhD thesis, University of Groningen. Appendix C: Minimum Image Convention.
2. Deiters, U.K. (1998). "Efficient Coding of the Minimum Image Convention." Molecular Simulation 20(4), 239-244.
3. LAMMPS `Domain::minimum_image` — `src/domain.cpp` lines 2321-2400
4. GROMACS `pbc_dx` / `low_set_pbc` — `src/gromacs/pbcutil/pbc.cpp`
5. Bekker, H. et al. (1997). "Efficient algorithms in computational chemistry." J. Comp. Chem. 18:1930-1942.
6. scicomp.stackexchange.com/questions/3107 — MIC for triclinic unit cell

# Blended Priority Metric Watershed — Implementation Plan

> **Status**: In Progress  
> **Created**: 2025-06-09  
> **Goal**: Fix metric-aware watershed segmentation for non-orthogonal crystal cells  

---

## Problem Statement

The current `rust_det_metric` watershed uses `priority = dgrid * inv_w` where `inv_w` is
an inverse metric weight per 26-neighbor offset. This **multiplicatively mixes physical
clearance (Å) with per-edge geometry**, creating a direction-dependent effective landscape
that misplaces boundaries away from physical constrictions.

### Symptoms
- **VIGNOQ** (triclinic, angles ~95-105°): `rust_det` (no weights) blows past boundaries;
  `rust_det_metric` with CORRECT step matrix (`cell.T @ diag(1/n)`) also places boundaries
  incorrectly — not at constrictions.
- **NU-1000** (hexagonal, γ=120°): `rust_det` shows directional preference; 
  `rust_det_metric` with correct step matrix fixes directional bias but the `dgrid*inv_w`
  coupling is still fundamentally wrong.

### Root Cause (confirmed by reviewer + oracle)
- The step matrix `cell.T @ diag(1/n)` IS mathematically correct (ASE: rows = lattice vectors,
  `cart = frac @ cell.array`, so `step = cell.T @ diag(1/n)`).
- The coupling `dgrid * inv_w` is wrong because:
  1. At constrictions (low dgrid), the correction is weakest — exactly where it matters most
  2. It creates a direction-dependent effective height landscape
  3. Method A (`cell @ diag(1/n)`) only appeared to work on VIGNOQ because its WRONG metric
     produced stronger weight anisotropy (1.77 vs 1.46) that compensated for the coupling error
     — a "two wrongs make a right" coincidence.

---

## Algorithm: Blended Priority Watershed

### Core Idea
Separate two concerns into one priority value:

```
priority = dgrid_saddle × SCALE − accumulated_physical_path_length
```

| Component | Propagation | Purpose |
|-----------|-------------|---------|
| `dgrid_saddle` | **min** (standard watershed) | Places boundary at physical saddle/constriction |
| `accumulated_path` | **sum** (Dijkstra) | Equalizes flood rate in physical Å space |
| `SCALE` | tunable, default 200 | dgrid resolution = 1/SCALE Å per Å of path |

### Why This Works
- **dgrid component** (×SCALE) dominates → boundaries at physical saddles
- **path component** provides additive metric correction → uniform regardless of dgrid value
- For SCALE=200: a path difference of 10 Å overcomes a dgrid difference of 0.05 Å
- Boundary shift from physical saddle: < 0.05 Å (negligible)
- Flood rate equalized over paths of 10+ Å (sufficient for metric correction)

### Pseudocode (Rust)

```rust
fn heap_watershed_blended_impl(
    image: Array3<f32>,          // dgrid
    seeds: Array3<u32>,          // labeled maxima
    mask: Option<Array3<u8>>,    // void mask
    edge_lengths: Vec<f32>,      // 26 physical edge lengths in Å (from metric)
    scale: f32,                  // default 200.0
) -> Array3<u32> {
    // Per-voxel state
    let mut labels: Vec<u32> = vec![0; total];
    let mut finalized: Vec<bool> = vec![false; total];
    let mut dgrid_saddle: Vec<f32> = vec![f32::NEG_INFINITY; total]; // min-propagated
    let mut path_len: Vec<f32> = vec![0.0; total];                   // sum-propagated
    let mut best_prio: Vec<f32> = vec![f32::NEG_INFINITY; total];    // blended priority

    // Max-heap (highest priority = closest to seed, highest dgrid)
    let mut heap: BinaryHeap<(NotNan<f32>, usize)> = BinaryHeap::new();

    // Initialize seeds
    for each seed (i,j,k) with label L {
        let idx = flat(i,j,k);
        labels[idx] = L;
        dgrid_saddle[idx] = image[idx];
        path_len[idx] = 0.0;
        best_prio[idx] = image[idx] * scale;
        heap.push((NotNan::new(image[idx] * scale).unwrap(), idx));
    }

    // Flood
    while let Some((_, idx)) = heap.pop() {
        if finalized[idx] { continue; }
        finalized[idx] = true;

        for (offset_idx, &(dx,dy,dz)) in NEIGHBOR_OFFSETS_26.iter().enumerate() {
            let nidx = neighbor_flat(idx, dx, dy, dz);
            if finalized[nidx] || labels[nidx] != 0 || mask[nidx] == 0 { continue; }

            // Min-propagation: dgrid saddle
            let cand_dgrid = image[nidx].min(dgrid_saddle[idx]);
            // Sum-propagation: physical path length
            let cand_path = path_len[idx] + edge_lengths[offset_idx];
            // Blended priority
            let cand_prio = cand_dgrid * scale - cand_path;

            if cand_prio > best_prio[nidx] {
                best_prio[nidx] = cand_prio;
                dgrid_saddle[nidx] = cand_dgrid;
                path_len[nidx] = cand_path;
                labels[nidx] = labels[idx];
                heap.push((NotNan::new(cand_prio).unwrap(), nidx));
            }
        }
    }

    Array3::from_shape_vec((nx, ny, nz), labels).unwrap()
}
```

### Edge Length Computation (Python + Rust)

```python
cell = atoms.cell.array  # ASE: rows = lattice vectors
step = cell.T @ np.diag(1.0/nx, 1.0/ny, 1.0/nz)  # CORRECT step matrix
for (di, dj, dk) in 26_neighbors:
    n = np.array([di, dj, dk], dtype=float)
    edge_length = np.linalg.norm(step @ n)  # physical distance in Å
```

No clamping needed — raw physical edge lengths in Å.

---

## Implementation Plan

### Phase 1: Rust Sandbox (`watershed_blended/`)

Create a standalone Rust crate **separate from `mgr_rust/`** so we can iterate without
touching the production codebase.

#### File Structure
```
watershed_blended/
├── Cargo.toml          # PyO3 + maturin
├── src/
│   └── lib.rs          # All code in one file (~200 lines)
└── pyproject.toml      # maturin build config
```

#### `Cargo.toml`
```toml
[package]
name = "watershed_blended"
version = "0.1.0"
edition = "2021"

[lib]
name = "watershed_blended"
crate-type = ["cdylib"]

[dependencies]
pyo3 = { version = "0.22", features = ["extension-module"] }
numpy = "0.22"
ndarray = "0.16"
```

#### `src/lib.rs` — Functions to expose

1. **`compute_edge_lengths(cell_matrix, grid_shape) -> Vec<f32>`**
   - Input: 3×3 cell matrix (ASE `cell.array`), `(nx, ny, nz)` grid shape
   - Step matrix: `step = cell.T @ diag(1/n)`
   - For each of 26 offsets: `edge_length = || step @ [di,dj,dk] ||`
   - Returns: 26 `f32` edge lengths in Å

2. **`seeded_watershed_det_blended(image, seeds, edge_lengths, scale=200.0, mask=None)`**
   - Input: 3D float32 image, 3D uint32 seeds, 26 edge lengths, scale param, optional mask
   - Algorithm: blended priority watershed (pseudocode above)
   - Returns: 3D uint32 labeled regions

3. **`seeded_watershed_det(image, seeds, mask=None)`** [for A/B comparison]
   - Standard unweighted heap watershed (same as `mgr_rust::seeded_watershed_det`)
   - Needed so we can compare `rust_det` from this module too

#### Build
```bash
cd watershed_blended
maturin develop --release
```

### Phase 2: Python Wrapper (`test_blended_watershed.py`)

A copy of `test_metric_copy.py` (which is itself a copy of `do_periodic_segmentation_rust`)
with these changes:

1. **Import `watershed_blended` instead of `mgr_rust`** for the watershed call only
   - Keep `mgr_rust` for everything else (dgrid, maxima detection, window connections)
   - The ONLY functions called from `watershed_blended`:
     - `watershed_blended.compute_edge_lengths()`
     - `watershed_blended.seeded_watershed_det_blended()`
     - `watershed_blended.seeded_watershed_det()` (for comparison)

2. **Add `--backend` CLI flag** with choices:
   - `rust_det` — unweighted standard (via `mgr_rust`)
   - `rust_det_metric` — current metric (via `mgr_rust`, `dgrid * inv_w`)
   - `rust_det_blended` — new blended priority (via `watershed_blended`)
   - Default: `rust_det_blended`

3. **Backend dispatch** — replace the single `BACKEND = "rust_det_metric"` with:
   ```python
   if args.backend == "rust_det_blended":
       # Use watershed_blended module
       edge_lengths = watershed_blended.compute_edge_lengths(
           cell_matrix=_cell_matrix_for_metric_weights(mat_atoms),
           grid_shape=dgrid.shape,
       )
       # ... padded_regions = watershed_blended.seeded_watershed_det_blended(...)
   elif args.backend == "rust_det_metric":
       # existing path (unchanged)
   elif args.backend == "rust_det":
       # existing path (unchanged)
   ```

4. **HTML output** with backend name in filename:
   - `VIGNOQ_clean_rust_det_blended.html`
   - `NU-1000_S1_rust_det_blended.html`

### Phase 3: Validation

#### Test Matrix

| Material | Cell Type | `rust_det` | `rust_det_metric` | `rust_det_blended` |
|----------|-----------|------------|--------------------|--------------------|
| VIGNOQ_clean | Triclinic (95-105°) | Blows past boundaries | Wrong boundaries (Method B) | Should be correct |
| NU-1000_S1 | Hexagonal (120°) | Directional preference | Fixed direction, wrong coupling | Should be correct |
| LTA96 | Cubic | Baseline (correct) | Should match | Should match |

#### Validation Criteria

1. **Boundary placement**: boundaries sit at dgrid constrictions (narrowest points)
2. **Self-loop bottleneck**: matches PLD ± 0.15 Å
3. **No directional preference**: PBC groups have symmetric voxel counts
4. **Indexing**: 0 mismatches between maxima_info labels and regions
5. **Visual**: HTML shows clean segmentation with boundaries at pore throats

#### SCALE Tuning

Default `SCALE = 200`. If needed, try:
- SCALE = 100 (stronger path correction, coarser dgrid resolution)
- SCALE = 400 (weaker path correction, finer dgrid resolution)
- SCALE = ∞ (degrades to pure unweighted watershed = `rust_det`)

### Phase 4: Integration (after validation)

1. **Port `heap_watershed_blended_impl` into `mgr_rust/src/lib.rs`**
   - Add alongside existing `heap_watershed_impl` and `heap_watershed_metric_impl`
   - Same structure, ~50 lines of new code

2. **Add Python function `compute_edge_lengths` to `mgr_rust`**
   - Rename from `compute_metric_inv_weights` or add new function
   - Return raw edge lengths instead of inv_w

3. **Wire `rust_det_blended` backend into `_seeded_watershed_backend`**
   - Add case in the dispatch function
   - Accept `scale` parameter (default 200.0)

4. **Update `do_periodic_segmentation_rust`**
   - Add `watershed_backend="rust_det_blended"` as new option
   - Update docstring

5. **Clean up**
   - Remove `seeded_watershed_det_metric` and `compute_metric_inv_weights` (superseded)
   - Or keep for backward compatibility with deprecation warning
   - Delete `watershed_blended/` sandbox

---

## Key Files (current state)

| File | Role |
|------|------|
| `mofography/rust_fun.py:256-297` | `_compute_metric_inv_weights` (will be replaced) |
| `mofography/rust_fun.py:300-342` | `_seeded_watershed_backend` dispatch (add case) |
| `mofography/rust_fun.py:33-41` | `_cell_matrix_for_metric_weights` (keep as-is) |
| `mgr_rust/src/lib.rs:1250-1296` | `metric_inv_weights_from_cell` (will be replaced) |
| `mgr_rust/src/lib.rs:1629-1731` | `heap_watershed_metric_impl` (will be replaced) |
| `test_metric_copy.py` | Current test harness (template for Phase 2) |
| `test_rust_det_metric_vignoq.py` | VIGNOQ validation script |

## ASE Cell Convention (confirmed)

- `atoms.cell.array`: **rows = lattice vectors** (cell.py line 17)
- `cart = frac @ cell.array` (cell.py line 113)
- `step = cell.array.T @ diag(1/nx, 1/ny, 1/nz)` ← CORRECT step matrix
- `cellpar_to_cell()` returns same convention (rows = lattice vectors)

## Mathematical Correctness Chain

1. Grid voxel `(i,j,k)` → fractional position `(i/nx, j/ny, k/nz)`
2. Fractional displacement `(di/nx, dj/ny, dk/nz)` → Cartesian displacement:
   `d_cart = [di/nx, dj/ny, dk/nz] @ cell.array = cell.array.T @ [di/nx, dj/ny, dk/nz]`
3. Step matrix: `step = cell.T @ diag(1/n)` → `d_cart = step @ [di, dj, dk]`
4. Edge length: `||step @ n||` for offset `n = [di, dj, dk]`
5. Blended priority: `dgrid_saddle × SCALE − Σ(edge_lengths)` along path

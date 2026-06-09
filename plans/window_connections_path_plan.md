# Window Connections Path Plan (VIGNOQ/REHHOD)

**Date:** 2026-06-07  
**Status:** Ready for implementation handoff (`designer` subagent)  
**Goal:** Add a new window extraction backend that is less sensitive to seam-shifted watershed boundaries, while keeping all existing segmentation outputs and interfaces unchanged.

---

## 1) Problem and constraints

- Keep current two-stage segmentation pipeline in `do_periodic_segmentation_rust` unchanged (merged and unmerged labels both remain available).
- Do **not** modify existing Rust function `window_connections(...)`.
- Add a **new** Rust function (proposed name: `window_connections_path`) and temporarily wire it in Python behind a method switch.
- New function must return data in the **exact same schema** as `window_connections`.

Current issue on `VIGNOQ_clean.cif`:
- Region-pair boundaries like `1 <-> 19` can occur at PBC seams (cell wall) in padded labels.
- Existing `window_connections` uses max-distance points on that boundary cluster, which can select wall windows.
- This inflates winding PLD relative to face-merging / zeo++ reference behavior.

---

## 2) Exact API contract (Rust)

Create a new PyO3 function in `mgr_rust/src/lib.rs`:

```rust
#[pyfunction]
#[pyo3(text_signature = "(labels, dgrid, padding_thickness, cell, distance_threshold)")]
fn window_connections_path(
    py: Python<'_>,
    labels: PyReadonlyArray3<'_, u32>,
    dgrid: PyReadonlyArray3<'_, f32>,
    padding_thickness: usize,
    cell: PyReadonlyArray2<'_, f64>,
    distance_threshold: f64,
) -> PyResult<PyObject>
```

### Output format (must be identical to `window_connections`)

Return a Python list of tuples, one tuple per region pair:

```python
(
  regs_arr,      # np.ndarray shape (2,), dtype uint32
  connected,     # bool (always True for emitted edges)
  edge_type,     # str, default "internal" (Python side may relabel periodic)
  centers_arr,   # np.ndarray shape (n_windows, 3), dtype float64 (Cartesian)
  radii_arr      # np.ndarray shape (n_windows,), dtype float64
)
```

Same tuple ordering and types as current Rust `window_connections`.

### Module registration

- Add to `#[pymodule] fn mgr_rust(...)`:
  - `m.add_function(wrap_pyfunction!(window_connections_path, m)?)?;`
- Add `"window_connections_path"` to `__all__`.

---

## 3) Algorithm for `window_connections_path` (fast path-based window)

For each touching region pair `(reg1, reg2)`:

1. **Candidate boundary voxel collection**
   - Reuse current 6-neighbor scan used by `window_connections` to identify touching pairs and collect boundary voxels.
   - Keep same pair map and clustering behavior (`distance_threshold`) to preserve multiplicity semantics.

2. **For each boundary cluster, compute a maximin path bottleneck**
   - Build a local ROI around cluster voxels:
     - `bbox(cluster_voxels) + margin` (margin in voxels, start with 3-5).
     - Clip to padded grid bounds.
   - Allowed path voxels:
     - `labels != 0` and in `{reg1, reg2}` plus immediate interface voxels, OR
     - alternatively `dgrid > 0` if first variant is too restrictive.
   - Start set = cluster voxels that belong to `reg1`.
   - End set = cluster voxels that belong to `reg2`.
   - Run widest-path (maximin) search in ROI (heap-based):
     - path score = minimum `dgrid` along the path.
     - maximize that score between start/end sets.
   - Output for this cluster:
     - `radius = best_score` (float64)
     - `center = Cartesian coordinate of bottleneck voxel` (or midpoint of two bottleneck-adjacent voxels if preferred)

3. **Fallback behavior**
   - If path search fails for a cluster, fallback to current local behavior for that cluster:
     - choose max-distance points per region in cluster (same as existing function).
   - This guarantees no edge drops due to path solver corner cases.

4. **Coordinate conversion**
   - Use same `to_cart` mapping currently used by `window_connections`:
     - convert padded index -> fractional using central shape (`nx - 2*padding_thickness`, etc.)
     - multiply by `cell` matrix.

Performance note:
- Keep search local (cluster ROI), not full-grid.
- Use 6-neighbor in path search first for speed and determinism.
- Parallelize across region pairs/clusters with Rayon (same style as existing function).

---

## 4) Python wiring plan (`mofography/rust_fun.py`)

In `do_periodic_segmentation_rust(...)`:

1. Add kwarg:
   - `window_method="boundary"` (default keeps current behavior).
2. At connection call site:
   - if `window_method == "boundary"` -> call `mgr_rust.window_connections(...)`
   - if `window_method == "path"` -> call `mgr_rust.window_connections_path(...)`
   - else raise `ValueError`.
3. Keep all downstream processing exactly unchanged:
   - periodic/internal tagging logic
   - remap_dict2/reindexing
   - output keys and shapes

No changes to default behavior.

---

## 5) Validation checklist (must pass)

### A) API/schema compatibility

- `window_connections_path` output tuples are accepted by existing downstream code without modifications.
- `do_periodic_segmentation_rust(..., window_method="path")` returns same keys and data structures as baseline.

### B) Numerical sanity

On `examples/TestCIFs/VIGNOQ_clean.cif` (`spacing=0.25`):
- Path windows should reduce/avoid wall-dominated seam windows for `1<->19`-type pair.
- Winding PLD should move toward face-merging (`get_pore_dia_mgr2`) and zeo++ reference (~3.11 A), not remain ~5.9 A.

On `examples/TestCIFs/REHHOD_clean.cif`:
- No catastrophic regression in edge count / graph build.
- Compare winding PLD and face-merging PLD deltas vs baseline.

### C) Speed benchmark

Measure per-material:
- total `do_periodic_segmentation_rust` time
- connection phase time only
- path vs boundary method ratio

Target (initial): path method total runtime <= 1.5-2.0x boundary method on test CIFs.

---

## 6) Suggested benchmark script behavior

For each material in `{VIGNOQ_clean.cif, REHHOD_clean.cif}`:
- Build `dgrid` once.
- Run segmentation twice with same params except `window_method`.
- Record:
  - `n_connections`
  - top periodic/self-loop windows
  - winding PLD
  - face-merging PLD (`get_pore_dia_mgr2`)
  - timings
- Emit compact table:
  - material, method, pld_winding, pld_face, ratio, seg_time_s, conn_time_s.

---

## 7) Risks and mitigations

- **Risk:** Path solver too slow on large clusters.
  - **Mitigation:** ROI-based search, early-stop when no better bottleneck possible, 6-neighbor first.
- **Risk:** Path cannot be found for odd cluster geometry.
  - **Mitigation:** per-cluster fallback to existing boundary-max logic.
- **Risk:** subtle schema drift.
  - **Mitigation:** strict tuple/type parity test against existing function.

---

## 8) Implementation punch list (ordered)

1. Add `window_connections_path` in `mgr_rust/src/lib.rs` (copy skeleton from `window_connections`).
2. Reuse pair detection + clustering from `window_connections`.
3. Implement ROI maximin path solver and bottleneck extraction.
4. Preserve identical output tuple construction.
5. Export function in `#[pymodule]` + `__all__`.
6. Add `window_method` switch in `do_periodic_segmentation_rust`.
7. Add small benchmark utility script for VIGNOQ + REHHOD.
8. Run validation table and store results in a short markdown note.

---

## 9) Decision gate after prototype

If `window_method="path"`:
- fixes VIGNOQ directionally (PLD closer to ~3.0) **and**
- does not break REHHOD **and**
- runtime acceptable,

then keep as experimental option and expand test set.

Otherwise keep boundary default and iterate on path ROI / solver tuning.

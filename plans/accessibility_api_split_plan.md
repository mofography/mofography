# Accessibility API split: winding vs lattice

## Goals

1. **Keep winding and lattice entirely separate** — no cross-fallback between methods.
2. **`window_fp` only on winding** — requires pore graph; not exposed on lattice APIs.
3. **`window_fp=False` means winding without FP** — not lattice (Chao-style strict graph percolation).
4. **Results dict stores both methods side-by-side** when computed:
   - `percolation_history_winding`
   - `percolation_history_lattice` (+ `percolation_lattice_variant`: `"v1"` | `"v2"`)
   - `pore_dia_winding` — `{lcd, pld, dfsd, method, ...}`
   - `pore_dia_lattice` — `{lcd, pld, dfsd, method, variant, ...}`
5. Legacy key `percolation_history` remains readable via helpers but new writes go to method-specific keys.

## Layers (orthogonal)

| Layer | Options |
|-------|---------|
| Accessibility | `winding` (default) or `lattice` (explicit) |
| Winding add-on | `window_fp=True/False` only when accessibility = winding |
| Geometry | `mesh` SA, `mc` SA, `voxel` void volume |

## Public API (target)

### Winding (production default)

```python
compute_surface_area_mesh_winding(results, atoms, probe_radius, *, window_fp=True, ...)
compute_surface_area_mc_winding(results, atoms, probe_radius, *, window_fp=True, ...)
compute_void_volume_from_results_winding(results, atoms, probe_radius, *, window_fp=True, ...)

compute_percolation_history_winding(results, ...)  # → percolation_history_winding
get_pore_dia_mgr_winding2(results=results, ...)    # → pore_dia_winding
```

### Lattice (Zeo++-style, explicit)

```python
compute_surface_area_mesh_lattice(results, atoms, probe_radius, *, variant="v2", ...)
compute_surface_area_mc_lattice(results, atoms, probe_radius, *, variant="v2", ...)
compute_void_volume_from_results_lattice(results, atoms, probe_radius, *, variant="v2", ...)

compute_percolation_history_from_results_lattice(results, variant="v2")  # → percolation_history_lattice
get_pore_dia_lattice(results, variant="v2", ...)  # alias of get_pore_dia_mgr2 → pore_dia_lattice
```

### Compatibility shims (deprecated)

```python
compute_surface_area_from_results(...)  # → _winding; window_fp_for_accessible controls FP only
compute_surface_area_from_results_mc(...)
compute_void_volume_from_results(...)
```

## Lattice variant ranking

| Rank | Variant | Functions | Notes |
|------|---------|-----------|-------|
| **1** | **v2** | `find_pore_limiting_radius2`, `percolation_all_axes_periodic`, `compute_percolation_history2` | Per-axis PBC union-find; **default for lattice** |
| **2 (deprecated)** | v1 | `find_pore_limiting_radius`, `compute_percolation_history` | Full-grid wrap pad + single 3D label; legacy only |

## Results dict contract

```python
results["percolation_history_winding"]  # list[(threshold, {x|y|z: [[region_ids], ...]})]
results["percolation_history_lattice"]    # same shape; v1 may use 0,1,2 keys
results["percolation_lattice_variant"]    # "v1" | "v2" when lattice history stored

results["pore_dia_winding"]  # {lcd, pld, dfsd, method: "winding", plr, ...}
results["pore_dia_lattice"]   # {lcd, pld, dfsd, method: "lattice", variant: "v2", ...}

# Legacy (read-only fallback in helpers):
results["percolation_history"]
```

Both histories and both pore-dia blocks may coexist on the same `results` object.

## Implementation phases

### Phase 1 (this PR)

- [x] `mofography/percolation_results.py` — keys, get/set helpers, pore-dia storage
- [x] Winding history writer → `percolation_history_winding`
- [x] Lattice history writers → `percolation_history_lattice`
- [x] Split accessibility helpers; fix `window_fp=False` → winding without FP
- [x] Public `*_winding` / `*_lattice` SA, MC SA, void vol entry points
- [x] Shims on old names → winding
- [x] `tests/test_accessibility_api_split.py` — Chao MC regression

### Phase 2 (follow-up)

- Rename `get_pore_dia_mgr` → `get_pore_dia_lattice_v1` with deprecation
- `lattice_percolation.py` module move
- Update viz/workflow to read method-specific keys via helpers
- Remove legacy `percolation_history` writes

## Testing

- Import Chao `compute_surface_area_from_results_mc_chao` from Downloads
- Compare region sets: Chao == `mc_winding(window_fp=False)` 
- `sum(by_region) == accessible_surface_area` for winding and lattice paths
- Both histories storable on same results without overwrite

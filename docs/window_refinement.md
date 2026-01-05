# Window Refinement Prototype

This note documents a distance-field-driven refinement pass for re-estimating
window centers/radii without relying solely on watershed boundaries. The
prototype operates on the full output of `do_periodic_segmentation_rust`.

## Required Inputs

- `results`: dict returned by `do_periodic_segmentation_rust` with
  `return_dgrid=True` and `return_padded_regions=True`
- `results["connections"]`: original window list
- `results["regions_padded"]`: segmentation in the padded domain
- `results["input_params"]["padding_thickness"]`
- `results["unique_regions"]`: optional filter (skip padded-only windows)
- `cell_matrix`: `3x3` array whose **columns** are the lattice vectors (with ASE,
  use `mat_atoms.cell.array.T`)

## Workflow

1. Filter connections (default: keep only those touching a unit-cell label).
2. For each window center, grab a small padded sub-volume around the center
   using the Cartesian→fractional mapping and padding thickness.
3. Build a binary mask containing voxels from the two regions forming the
   connection. All computations stay inside this local block to keep the
   runtime manageable.
4. Produce new center/radius estimates with the morphological opening heuristic.
5. Return a new list of connection tuples matching the original format so
   downstream code can swap the array without extra bookkeeping.

## Morphological Opening via Distance Thresholds

Implemented as `refine_windows_via_morph_opening` in
`mofography/window_refinement.py`.

1. Perform a binary search on the distance threshold `t` (default 6 steps) to
   find the largest `t` where the voxels belonging to either region with
   `dgrid >= t` still remain connected.
2. When the set disconnects, gather the surviving voxels of each region and
   compute their centroids in global voxel coordinates.
3. Convert the centroids to Cartesian coordinates and average them to obtain
   the refined window center. Use the limiting threshold as the new radius.

This mimics a morphological opening driven purely by the scalar distance field
and typically pushes the window toward the tightest throat consistent with the
local geometry.

[What is a morphological opening?](https://en.wikipedia.org/wiki/Opening_(morphology))
In classical image processing, an opening is an erosion followed by a dilation
with the same structuring element—it trims protrusions while preserving the
overall mass of the object. Our distance-grid variant plays the same role:

- The threshold `t` is analogous to the structuring element size. We gradually
  erode the pore throats by demanding `dgrid >= t`, which peels away voxels near
  the walls.
- Once the two regions disconnect, the remaining voxels are the last “bridge”
  of free volume between them—the actual window.
- Averaging the surviving centroids from both sides effectively places the
  window center at the saddle of this bridge, while the threshold itself serves
  as the radius estimate.

This rigorous approach requires no extra boundary data and adapts to curved or
noisy interfaces, making it suitable for amorphous pores.

### Mini 2D Example

Consider a binary slice of a pore throat (1 = void, 0 = framework) where the
opening narrows in the middle. Our 3×3 structuring element (SE) contains ones.

```
Slice (t=0.0):            SE:
0 0 0 0 0 0 0 0           1 1 1
0 1 1 1 1 1 1 0           1 1 1
0 1 1 1 1 1 1 0
0 1 1 1 1 1 1 0
0 1 1 1 1 1 1 0
0 0 0 1 1 1 0 0  <-- throat narrows to 3 voxels
0 0 0 1 1 1 0 0
0 0 0 0 0 0 0 0
```

1. **Erode once (t≈1):**

```
0 0 0 0 0 0 0 0
0 0 1 1 1 1 0 0
0 0 1 1 1 1 0 0
0 0 1 1 1 1 0 0
0 0 0 1 1 1 0 0
0 0 0 1 1 1 0 0
0 0 0 0 0 0 0 0
```

2. **Erode twice (t≈2):**

```
0 0 0 0 0 0 0 0
0 0 0 1 1 0 0 0
0 0 0 1 1 0 0 0
0 0 0 1 1 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
```

3. **Erode thrice (t≈3):** the structure collapses entirely, meaning the previous
   step captured the last viable bridge.
4. **Dilate back (conceptually):** expanding the surviving column recreates the
   3-voxel-wide bridge—the window location.

In our 3D distance grid we mimic these steps by:

- Treating `dgrid` as a continuous analog of the binary mask.
- Thresholding by `t` (erosion) to keep only voxels far enough from the
  framework.
- Taking connectivity centers after the threshold crosses the critical point
  (dilation).

Thus the morph-opening effectively “slides” a spherical probe of radius `t`
through the interface and records the last location where it can pass between
the two regions—precisely the window we care about.

### Sub-voxel placement?

The prototype currently accepts the discrete centroid of the surviving voxels,
so the refined center lies on the grid. Interpolation can be layered on top
later (e.g., tri-linear interpolation within the local block or fitting a local
quadratic to `dgrid`) if sub-voxel accuracy becomes necessary, but the baseline
implementation stays grid-aligned for speed/simplicity.

The function returns a new connection array and the caller can optionally stash
it inside the segmentation dictionary as `results["connections_morph_opening"]`
for downstream plotting.

## Example Usage

```python
from mofography.window_refinement import refine_windows_via_morph_opening

results = mgr.do_periodic_segmentation_rust(
    dgrid,
    mat_atoms,
    return_dgrid=True,
    return_padded_regions=True,
    compute_connections=True,
)

cell_matrix = mat_atoms.cell.array.T  # ASE returns rows; transpose to columns

connections_morph = refine_windows_via_morph_opening(
    results,
    cell_matrix=cell_matrix,
)

results["connections_morph_opening"] = connections_morph
```

The refined list can replace `results["connections"]` before building graphs or
trimming edges. The code is a pure-Python prototype intended for experimentation
prior to a Rust port.

## Testing Outputs

`docs/test_window_refinement.py` automates regression checks on
`NU-1000_S1.cif`. The script regenerates histograms, parity plots, and Plotly
scenes for the baseline and refined windows:

- With the default parameters (`spacing=0.5`, `mask_thickness=0.5`, `h=0.2`,
  etc.), the morphological pass shifts just a small fraction of windows:
  mean center deviation ≈0.15 Å, max shift ≈113 Å (rare padded cases), and
  mean radius change ≈0.0003 Å.
- Generated artifacts are written to
  `docs/window_refinement_plots/`:
  - `center_shift_hist.html`, `radius_delta_hist.html`
  - `parity_morph.html`
  - `scene_baseline.html`, `scene_morph_opening.html`

Open those HTML files in a browser to inspect how the refinement nudges each
window while leaving most IRMOF/NU-1000 connections untouched.


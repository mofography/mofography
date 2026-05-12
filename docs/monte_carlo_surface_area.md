### Monte Carlo shell estimator for probe surface area

This note describes the experimental Monte Carlo (MC) method we use to estimate
probe–surface area from an existing distance grid. The goal is to reduce
sensitivity to global grid resolution and to avoid some of the artifacts of
marching–cubes–based surface meshes, while still reusing the full
segmentation/percolation pipeline.

#### Basic idea

Given a periodic unit cell and a distance–to–framework field \(f(x)\) (the
`dgrid`), the probe surface at radius \(r\) is the level set
\(f(x) = r\). Rather than explicitly meshing this level set, we:

1. Draw \(N\) random points \(x_i\) uniformly in the unit cell.
2. For each point, interpolate \(f(x_i)\) from the existing `dgrid`.
3. Define a thin shell of thickness \(\delta\) around the probe radius,
   i.e. points such that \(|f(x_i) - r| < \delta/2\).
4. Let \(p_{\text{shell}}\) be the fraction of samples in this shell.

For small \(\delta\), the surface area can be estimated by the identity
`SA ≈ V_cell * p_shell / δ`, where \(V_\mathrm{cell}\) is the unit–cell volume.
We convert from Å² to m²/g
using the same normalization as the existing code:

- Compute \(V_\mathrm{cell}\) in Å³ from `mat_atoms.get_volume()`.
- Compute the framework mass \(m\) in grams from the atomic masses.
- Apply the factor \(1 \text{ Å}^2 = 10^{-20} \text{ m}^2\) and divide by \(m\)
  to report SA in m²/g.

This gives a global **total surface area** estimate. To get **accessible**
surface area, we also use the per–region labels and percolation history.

#### Reusing existing data (`dgrid`, regions, percolation)

The MC method is designed to sit on top of the existing segmentation
infrastructure. For each structure we reuse:

- `results_dict["dgrid"]`:
  - 3D array of distances to the framework.
  - Used as a scalar field \(f(x)\); we never recompute distances.
- `results_dict["regions_pbc"]`:
  - 3D array of integer region labels on the same grid.
  - We interpolate labels at sample points with nearest–neighbor interpolation.
- `results_dict["percolation_history"]` (or we compute it if missing):
  - We call `select_percolation_threshold(percolation_history, probe_radius)`
    to obtain the percolating region IDs at the desired probe.
  - This is the same helper used in `compute_surface_area_from_results`.
- `mat_atoms`:
  - Provides the periodic cell and total mass, exactly as in the existing SA
    and void–volume routines.

Interpolations are done with `scipy.ndimage.map_coordinates`:

- `order=1` and `mode="wrap"` for `dgrid` (trilinear interpolation with
  periodic wrapping).
- `order=0` and `mode="wrap"` for `regions_pbc` (nearest–neighbor, categorical
  labels).

Accessible SA is obtained by restricting shell samples to those whose
interpolated region label is in the percolating–region set at the probe
radius. The estimator then becomes
`SA_acc ≈ V_cell * p_shell_acc / δ`, where \(p_{\text{shell,acc}}\) is the
fraction of random points in the shell
that also belong to percolating regions.

#### Relation to the marching–cubes method

The existing `compute_surface_area_from_results` function computes SA by
constructing a triangular mesh of the isosurface \(f(x) = r\) using marching
cubes, then measuring the mesh area with `trimesh`. It then splits this mesh
into per–region submeshes and derives per–region / per–pocket / per–cluster
SA from those pieces.

That approach is deterministic and fast, but it has two important sensitivities:

- **Grid resolution**: coarse `dgrid` spacing can lead to blocky surfaces and
  biased areas.
- **Region splitting at interfaces**: faces that straddle region boundaries
  can be dropped when constructing per–region submeshes, which can undercount
  surface area, especially for complex pore morphologies.

The MC estimator avoids building a mesh altogether. Instead it directly
estimates how much of the cell volume lies “near” the level set, which is
mathematically well–defined even when the discrete mesh is imperfect.

In practice we see better agreement with Zeo++ for SA because:

- The MC method is **less sensitive to missing triangles** at region
  interfaces and other meshing artifacts.
- Zeo++ itself uses a sampling–based approach for SA, so our estimator is
  conceptually closer to its measurement procedure than marching cubes on a
  coarse grid.

#### What the MC estimator returns

In its current form, the experimental MC routine focuses on global totals:

- `total_surface_area` (SA over all regions, MC estimate).
- `accessible_surface_area` (SA over percolating regions only, MC estimate).
- `sa_mc_diagnostics`:
  - `n_samples`
  - `shell_thickness`
  - `n_shell_hits_total`
  - `n_shell_hits_accessible`

The void–volume calculation is left unchanged and still uses the existing
grid–based method (`compute_void_volume_from_results`), since it already
agrees well with Zeo++.

The marching–cubes–based function remains available and continues to provide
rich breakdowns:

- `by_region` (per–region accessible SA),
- `by_pocket` (per–PBC–group SA),
- `by_percolated_cluster` (per–cluster SA),

which are useful for visualization and qualitative analysis. The MC estimator
is intended as a complementary, more robust global SA estimator that you can
use alongside the existing breakdowns and Zeo++ parity plots.

#### Implementation details

- **Location in codebase**
  - Core implementation: `mofography/compute.py`, function
    `compute_surface_area_from_results_mc(results_dict, mat_atoms, ...)`.
  - This function is part of the public `mofography.compute` API and is
    intended to be a drop-in alternative to
    `compute_surface_area_from_results` for SA.

- **Inputs**
  - `results_dict`:
    - Must contain at least:
      - `'dgrid'` or `'distance_grid'`: distance grid (Å).
      - `'regions_pbc'`: periodic region labels (integers).
      - `'pbc_groups'`: list or dict of PBC pockets (each group is a set of
        region IDs).
    - Optionally:
      - `'percolation_history'`: list of `(threshold, components)` tuples. If
        missing it is computed via `compute_percolation_history_from_results`.
  - `mat_atoms`:
    - ASE `Atoms` object for the (super)cell used in the distance grid.
    - Provides `get_volume()` and atomic masses for normalization to m²/g.
  - `probe_radius`:
    - Probe radius (Å) for which SA is estimated.
  - `n_samples`, `shell_thickness`:
    - Monte Carlo sample count and shell half-width control accuracy and cost.

- **Outputs**
  - Global SA:
    - `total_surface_area`: MC estimate over all regions (m²/g).
    - `accessible_surface_area`: MC estimate over percolating regions only
      (m²/g).
    - `accessible_fraction`: ratio of accessible to total SA.
  - Breakdowns (accessible SA only, mirroring the marching-cubes function):
    - `by_region`: dict `{region_id -> SA}` for percolating regions.
    - `by_pocket`: dict `{pocket_id -> SA}` from summing per-region SA over
      each PBC group.
    - `by_percolated_cluster`: dict `{cluster_id -> SA}` from summing
      per-region SA over each percolated cluster at the selected probe-threshold.
  - Diagnostics:
    - `sa_mc_diagnostics` with:
      - `n_samples`, `shell_thickness`,
      - `n_shell_hits_total`, `n_shell_hits_accessible`.



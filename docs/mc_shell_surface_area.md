# MC Shell Surface-Area Estimator

> *Quick reference for `compute_surface_area_mc_winding` in `mofography/compute.py`.*

## TL;DR

The Monte-Carlo shell estimator turns surface area into a *counting* problem:

1. Throw `N` random points inside the unit cell.
2. Keep the points that land in a thin shell of thickness `δ` centred on the probe surface
   (i.e. within `δ/2` of the isosurface `d(x) = r_probe`).
3. Compute `SA = V_cell · (n_shell / N) / δ`, then convert Å² → m²/g.

That's it. Everything else is bookkeeping around step 2 (percolation filtering, region
labelling, mass normalisation).

---

## The Core Idea

A 2-D surface *smeared* into a thin 3-D band has a volume proportional to its area:

```
V_shell  ≈  SA · δ          (in the limit of small δ)
```

So if you can measure `V_shell` and you know `δ`, you get `SA = V_shell / δ`.

The trick is measuring `V_shell` without explicitly constructing the surface. The
distance grid `d(x)` (already computed during segmentation) gives us the surface
implicitly: every point where `d(x) = r_probe` is *on* the probe surface, and every
point where `|d(x) - r_probe| ≤ δ/2` is in the shell.

Monte Carlo gives the shell volume:

```
V_shell ≈ V_cell · (n_shell / N)         (fractional volume of the shell)
```

Combining:

```
SA = V_cell · p_shell / δ       where   p_shell = n_shell / N
```

## What the Code Actually Does

Pseudocode of `compute_surface_area_mc` in `mofography/compute.py`:

```python
# 1. Random uniform points in fractional coords
u = np.random.rand(n_samples, 3)                       # (N, 3) in [0, 1)^3
coords = u * (shape - 1)                               # (N, 3) in voxel units

# 2. Interpolate the distance grid at each point
d_samples = map_coordinates(distance_grid, coords, order=1, mode="wrap")

# 3. Define the shell: points whose distance-to-wall is within ±δ/2 of r_probe
half_delta = shell_thickness / 2.0
shell_mask = np.abs(d_samples - probe_radius) <= half_delta
n_shell_total = shell_mask.sum()

# 4. Estimate SA
p_shell = n_shell_total / n_samples
total_surface_area = V_cell * p_shell / shell_thickness * 1e-20 / mass_g
```

The accessibility version additionally filters samples by whether their region
label is in `percolating_regions`, giving `accessible_surface_area` from the same
samples.

## Parameters

| Parameter | Default | What it does | Effect of changing |
|---|---|---|---|
| `n_samples` | 20 000 | Number of random points per structure | More samples → lower statistical noise (∝ 1/√N) |
| `shell_thickness` (`δ`) | 0.5 Å | Half-width of the band around the probe surface | Thicker → more hits per sample (lower variance) but more bias from averaging |
| `probe_radius` | He radius (1.2 Å) | Radius of the fictitious probe molecule | Larger probe → smaller accessible surface |

### Bias–Variance Tradeoff

The expected value of the estimator is essentially **independent of `δ`** (in the
small-δ limit), because `p_shell` scales linearly with `δ` and the `1/δ` factor
cancels. What *does* depend on `δ` is the **variance** and the **bias** from
smoothing:

- **Thin shell** (δ ≈ 0.2 Å): very little bias, but few samples per shell, so the
  estimate is noisy. Needs `n_samples` ≈ 50 000+ to reach the same precision.
- **Default** (δ = 0.5 Å): well-matched to a `SPACING = 0.5 Å` grid; a good
  default for screening.
- **Thick shell** (δ ≈ 1.0 Å): most samples per shell, but the SA is averaged over
  a band that is itself a sizable fraction of small pores. Only sensible when
  the distance field is locally flat (large pores).

The default `δ = 0.5 Å` is a good choice because it is a few times larger than
the grid spacing — so the shell contains a statistically meaningful number of
samples — and it is small compared to typical LCD values, so the bias is
negligible.

## Sanity Check

A quick way to confirm that `shell_thickness` isn't biasing your results is to
re-run the same MOF with two different `δ` values and check the agreement:

```python
for delta in (0.3, 0.5, 0.7):
    sa = mgr.compute_surface_area_mc_winding(
        results, atoms,
        probe_radius=1.2,
        n_samples=40_000,
        shell_thickness=delta,
    )
    print(f"δ = {delta:.1f} Å  →  SA = {sa['total_surface_area']:.1f} m²/g")
```

For typical MOFs the three numbers should agree to within a few percent — the
MC noise level for 40 000 samples.

## When to Use MC vs Mesh

| Method | Function | Pros | Cons |
|---|---|---|---|
| **Mesh** | `compute_surface_area_mesh_winding` | Deterministic; exact (modulo tessellation) | Slow for large grids; needs explicit surface triangulation |
| **MC shell** | `compute_surface_area_mc_winding` | Fast; trivially parallel over `n_samples`; tunable precision | Stochastic; needs enough samples for low noise |

The MC shell method is the natural choice for high-throughput screening, where
you can trade a few-percent statistical error for a 10×+ speedup over the mesh
method.

## References

- Implementation: `mofography/compute.py::compute_surface_area_mc` (≈ line 9524)
- Winding-accessible wrapper: `mofography/compute.py::compute_surface_area_mc_winding` (≈ line 10097)
- See also: `docs/accessibility_and_surface_area.org` for the broader SA/Vol API

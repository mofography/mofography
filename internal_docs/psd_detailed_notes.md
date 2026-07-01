# Pore Size Distribution (Gelb & Gubbins) - Detailed Notes

> **Technical documentation and advanced usage of Gelb & Gubbins PSD method**

---

## Overview

This document provides detailed technical information about the Gelb and Gubbins pore size distribution (PSD) method, its implementation in MOFography, and advanced usage patterns.

---

## Theoretical Background

### Gelb & Gubbins Method

**Original Paper**: Gelb, L.D., & Gubbins, K.E. (1998). "Characterization of porous solids and pore-size distributions using a mean-field method." *Langmuir*, 14(8), 1887-1902.

### Core Algorithm

For any point **P** in the void space of a porous material:

1. **Generate Distance Field**: Compute distance from framework atoms at every grid point
2. **Find Largest Inscribed Sphere**: For point **P**, find the largest sphere that:
   - Centered at some grid point **C** (not necessarily at P)
   - Sphere radius = distance from grid point **C** to nearest framework atom
   - **C** must be at distance ≥ sphere radius from **P** (sphere includes **P**)
3. **Repeat for Many Points**: Sample many points throughout the void space
4. **Create Histogram**: Bin all sphere diameters to create distribution

### Mathematical Definition

```
For each sample point p in void space:
    D(p) = max{ 2 × d(C, nearest_atom) | d(p, C) ≥ d(C, nearest_atom) }

where:
    d(C, nearest_atom) = distance from grid point C to framework atoms
    d(p, C) = Euclidean distance between sample point p and grid point C
```

**Key Properties:**
- **Pore diameter** = 2 × radius of largest inscribed sphere
- **Non-local**: Largest sphere center may be different from sample point
- **Geometrically exact**: Based on distances to atoms, not approximations

---

## Implementation in MOFography

### Rust Backend Functions

#### 1. `compute_pore_sizes_uniform()`

**Python Wrapper**: `compute_pore_sizes_uniform_rust(dgrid, ase_atoms, num_samples_total)`

**Purpose**: Uniform sampling across entire unit cell

**Algorithm**:
1. Generate `num_samples_total` random points uniformly in unit cell
2. For each sample point:
   - Sort all grid voxels by distance value (descending)
   - Iterate from largest distance to smallest
   - Find first voxel where sphere includes sample point
   - Return that sphere diameter

**Parameters**:
- `dgrid`: 3D distance grid (must be float32)
- `ase_atoms`: ASE Atoms object with cell and positions
- `num_samples_total`: Number of random points (default: 10000)

**Output**:
- `sample_points`: Array of shape (N, 3) - sample coordinates (fractional)
- `pore_sizes`: Array of shape (N,) - largest sphere radii (in Angstroms)

**Performance**:
- Parallelized using Rayon
- Sorted voxel iteration for efficiency
- Early termination when sphere found

---

#### 2. `compute_pore_sizes_for_psd_pocketwise()`

**Python Wrapper**: `compute_pore_sizes_pocketwise_rust(dgrid, regions_pbc, ase_atoms, samples_per_region)`

**Purpose**: Per-region sampling within each identified pocket

**Algorithm**:
1. Find unique region labels in `regions_pbc`
2. For each region:
   - Generate `samples_per_region` random points WITHIN that region
   - For each point:
     - Sort voxels by distance (descending)
     - Find first voxel where sphere includes point
     - Record sphere diameter
3. Return all points, pore sizes, and region labels

**Parameters**:
- `dgrid`: 3D distance grid (float32)
- `regions_pbc`: 3D array of region labels from periodic segmentation
- `ase_atoms`: ASE Atoms object
- `samples_per_region`: Number of points per region (default: 1000)

**Output**:
- `sample_points`: Array of shape (N_total, 3) - sample coordinates
- `pore_sizes`: Array of shape (N_total,) - largest sphere radii
- `region_labels`: Array of shape (N_total,) - corresponding region IDs

---

## Python Wrapper Functions

### 1. `compute_psd_uniform_rust()`

**Location**: `mofography/rust_fun.py:1377`

**Signature**:
```python
def compute_psd_uniform_rust(
    dgrid,
    ase_atoms,
    num_samples_total=10000,
    bins=None,
    density=False,
    return_fig=False,
    material="Material"
)
```

**Usage Example**:
```python
# Get histogram and optional figure
psd_hist, fig = mgr.compute_psd_uniform_rust(
    dgrid, mat_atoms,
    num_samples_total=10000,
    bins=np.linspace(0, 50, 200),
    density=False,
    material='NU-1000',
    return_fig=True
)

# Get raw data for custom analysis
sample_points, pore_sizes = mgr.compute_pore_sizes_uniform_rust(
    dgrid, mat_atoms,
    num_samples_total=10000
)

# pore_sizes are radii, multiply by 2 for diameters
pore_diameters = pore_sizes * 2
```

**Returns**:
- If `return_fig=True`: `(psd_hist, go.Figure)`
- If `return_fig=False`: `psd_hist` (histogram values only)

---

### 2. `compute_psds_pocketwise_rust()`

**Location**: `mofography/rust_fun.py:1483`

**Signature**:
```python
def compute_psds_pocketwise_rust(
    dgrid,
    regions_pbc,
    ase_atoms,
    num_samples_per_region=1000,
    bins=None,
    return_fig=False,
    material="Material"
)
```

**Usage Example**:
```python
# Get PSDs per pocket and optional figure
pocket_psds, fig = mgr.compute_psds_pocketwise_rust(
    dgrid, regions_pbc, mat_atoms,
    num_samples_per_region=1000,
    bins=np.linspace(0, 50, 200),
    material='NU-1000',
    return_fig=True
)

# pocket_psds is a dictionary:
# {
#     region_id_1: {
#         'hist': array of histogram values,
#         'bin_centers': array of bin centers,
#         'num_samples': number of samples in this region,
#         'mean_pore_size': mean pore diameter,
#         'std_pore_size': std dev of pore diameters
#     },
#     region_id_2: { ... },
#     ...
# }
```

---

## Helper Functions

### `compute_psds_from_points_pocketwise()`

**Location**: `mofography/helper.py:1075`

**Purpose**: Convert raw pore sizes and labels into per-pocket histograms

**Algorithm**:
1. Find unique region labels
2. Filter out background/invalid labels
3. For each region:
   - Extract pore sizes belonging to region
   - Filter out infinite/invalid values
   - Skip if insufficient samples (< 10)
   - Compute histogram
   - Calculate statistics (mean, std)
4. Return dictionary of pocket PSDs

**Parameters**:
- `sample_points_labels`: Array of region labels
- `pore_sizes`: Array of pore diameters
- `bins`: Bin edges for histogram

**Returns**:
```python
{
    region_id: {
        'hist': histogram_values,
        'bin_centers': bin_centers,
        'num_samples': count,
        'mean_pore_size': mean_diameter,
        'std_pore_size': std_diameter
    }
}
```

---

### `plot_pocket_psds()`

**Location**: `mofography/visualize.py:3873`

**Purpose**: Create comparison visualization of overall vs per-pocket PSDs

**Features**:
- Overall histogram (gray bars)
- Individual pocket curves (colored lines)
- Statistics in hover text
- Publication-ready formatting

---

## Parameter Guidelines

### Sampling Parameters

| Parameter | Global PSD | Pocket-wise PSD | Recommended Range |
|-----------|-------------|------------------|-------------------|
| `num_samples_total` | 10,000 | N/A (use `num_samples_per_region`) | 5,000 - 50,000 |
| `num_samples_per_region` | N/A | 1,000 | 500 - 5,000 |
| `bins` | 200 | 200 | 100 - 500 |
| `spacing` (dgrid) | 0.25 - 0.5 Å | 0.25 - 0.5 Å | Finer = more accurate but slower |

**Trade-offs**:
- **More samples**: Better statistics, but slower computation
- **More bins**: Higher resolution, but noisier histograms
- **Finer spacing**: More accurate pore sizes, but larger grid size

---

### Accuracy vs Performance

| Aspect | CSD (fast) | PSD (slow) | Recommendation |
|---------|-------------|-------------|-----------------|
| **Speed** | Very fast | Slow (10-100x) | Use CSD for screening |
| **Accuracy** | Approximate | Accurate | Use PSD for publication |
| **Memory** | Low | Moderate | Depends on grid size |
| **Use Case** | Pocket detection | Detailed analysis | Choose by goal |

---

## Advanced Usage

### Custom Binning

```python
# Use logarithmic bins for wide pore size ranges
import numpy as np
bins = np.logspace(np.log10(0.5), np.log10(50), 100)

psd_hist, fig = mgr.compute_psd_uniform_rust(
    dgrid, mat_atoms,
    num_samples_total=10000,
    bins=bins,
    density=False
)
```

### Analyzing Pocket Heterogeneity

```python
# Compare variance across pockets
pocket_psds = compute_psds_pocketwise_rust(...)

# Calculate coefficient of variation
means = [psd['mean_pore_size'] for psd in pocket_psds.values()]
stds = [psd['std_pore_size'] for psd in pocket_psds.values()]

# Heterogeneous: high coefficient of variation
if np.std(means) / np.mean(means) > 0.3:
    print("Material is highly heterogenous!")
    print("  Consider using pocket-wise analysis")
```

### Exporting for Analysis

```python
# Export to CSV for further analysis
import pandas as pd

# Global PSD
df_global = pd.DataFrame({
    'pore_diameter': pore_diameters
})
df_global.to_csv('global_psd.csv', index=False)

# Pocket-wise PSD
rows = []
for region_id, psd_data in pocket_psds.items():
    for diameter in range(len(psd_data['hist'])):
        rows.append({
            'region_id': region_id,
            'pore_diameter': psd_data['bin_centers'][diameter],
            'count': psd_data['hist'][diameter]
        })

df_pockets = pd.DataFrame(rows)
df_pockets.to_csv('pocketwise_psd.csv', index=False)
```

---

## Troubleshooting

### Common Issues

#### 1. Rust Backend Not Available

**Error**: `ModuleNotFoundError: No module named 'mgr_rust'`

**Solution**:
```bash
# Compile Rust backend
cd mgr_rust
maturin build --release

# Install
pip install --force-reinstall target/wheels/*.whl

# Reinstall mofography
cd ..
python setup.py install
```

#### 2. dtype Mismatch

**Error**: `ValueError: expected f32, got f64`

**Solution**:
```python
# Convert distance grid to float32
dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=0.25)
dgrid = dgrid.astype(np.float32)  # CRITICAL for Rust backend
```

#### 3. Memory Issues

**Error**: `MemoryError: Unable to allocate array`

**Solutions**:
1. Reduce `num_samples_total` or `num_samples_per_region`
2. Increase grid spacing (0.5 → 0.75 Å)
3. Use smaller material or unit cell
4. Close other applications before running

#### 4. No Pockets Found

**Error**: `Found 0 unique regions` or empty pocket_psds

**Causes**:
- Segmentation parameters too aggressive
- Material has no accessible void space
- Threshold too high (only background detected)

**Solutions**:
1. Lower `threshold_abs` in `do_periodic_segmentation_rust()`
2. Lower `mask_thickness` to include more void space
3. Increase `extra_padding` to capture boundary regions
4. Check structure with visualization first

#### 5. Unexpected Results

**Issue**: Pore sizes seem too small or distribution looks wrong

**Debugging Steps**:
```python
# 1. Visualize distance grid
import matplotlib.pyplot as plt
plt.imshow(dgrid[:,:,dgrid.shape[2]//2].T)
plt.colorbar(label='Distance to atoms (Å)')
plt.title('Distance Grid Check')
plt.show()

# 2. Check sample point locations
print("Sample point ranges:")
print(f"  X: [{sample_points_uniform[:,0].min():.2f}, {sample_points_uniform[:,0].max():.2f}]")
print(f"  Y: [{sample_points_uniform[:,1].min():.2f}, {sample_points_uniform[:,1].max():.2f}]")
print(f"  Z: [{sample_points_uniform[:,2].min():.2f}, {sample_points_uniform[:,2].max():.2f}]")

# 3. Check pore size distribution
print(f"Pore size statistics:")
print(f"  Min: {np.min(pore_sizes)*2:.2f} Å")
print(f"  Max: {np.max(pore_sizes)*2:.2f} Å")
print(f"  Mean: {np.mean(pore_sizes)*2:.2f} Å")
print(f"  Median: {np.median(pore_sizes)*2:.2f} Å")
```

---

## Performance Optimization

### Parallel Processing

Rust backend automatically uses Rayon for parallel processing:

```python
# Check number of threads being used
import os
n_threads = os.environ.get('RAYON_NUM_THREADS', 'auto-detected')
print(f"Using {n_threads} CPU threads")
```

### Memory Management

For large structures (> 10,000 atoms):

1. **Use Coarse Grid**: spacing = 0.5 Å instead of 0.25 Å
2. **Reduce Samples**: num_samples_total = 5,000 instead of 10,000
3. **Process in Batches**: Split analysis into smaller chunks
4. **Use Dask**: For very large datasets

### Benchmarking

```python
import time

# Time the operations
start = time.time()
sample_points, pore_sizes = mgr.compute_pore_sizes_uniform_rust(
    dgrid, mat_atoms, num_samples_total=10000
)
elapsed = time.time() - start

print(f"Computed {len(pore_sizes)} pore sizes in {elapsed:.2f} seconds")
print(f"Rate: {len(pore_sizes)/elapsed:.0f} samples/sec")
```

Typical performance (NU-1000, spacing=0.25 Å):
- Global PSD (10,000 samples): ~5-10 seconds
- Pocket-wise PSD (12 regions, 1000 samples each): ~30-60 seconds
- Speedup vs pure Python: 8-12x

---

## Comparison with Other Methods

### Gelb-Gubbins vs Voronoi-Based PSD

| Aspect | Gelb-Gubbins | Voronoi-Based | Note |
|---------|-------------|----------------|------|
| **Basis** | Distance grid | Voronoi tessellation | Different geometric models |
| **Accuracy** | High | Medium | G-G more accurate for non-spherical pores |
| **Speed** | Fast | Fast | Similar performance |
| **Implementation** | Rust | Python | Voronoi in pure Python |
| **PBC Handling** | Native | Complex | G-G simpler |

### PSD vs CSD for Pocket Detection

When choosing features for clustering pockets:

**Use CSD when**:
- Quick screening of many materials
- Computational resources limited
- Approximate geometric differences acceptable
- Similar pore sizes differ mainly by size

**Use PSD when**:
- Publication-quality analysis
- Materials have complex pore geometries
- Need to distinguish pockets by shape, not just size
- Sensitive to corner/edge effects
- Computational resources available

**Example**:
- IRMOF-1: CSD sufficient (2 distinct pocket types)
- NU-1000: PSD recommended (complex 1D channels)
- HKUST-1: PSD recommended (3 different cage types)

---

## Reference Implementation

### Rust Algorithm Details

**File**: `mgr_rust/src/lib.rs:5420` (find_largest_including_sphere_uniform)

```rust
fn find_largest_including_sphere_uniform(
    dgrid: &ndarray::ArrayView3<'_, f32>,
    cell: &ndarray::ArrayView2<'_, f64>,
    sample_point: &[f64; 3],
    sorted_voxels: &Vec<(f32, [usize; 3])>,
) -> f32
```

**Key Steps**:
1. Pre-compute cell components for PBC
2. Iterate through sorted voxels (highest distance first)
3. Convert grid coords to fractional coords
4. Calculate PBC-aware fractional distance
5. Convert to Cartesian distance
6. Check if sphere radius >= distance to nearest atom
7. Return first valid sphere diameter

**Complexity**: O(N × M) where N = samples, M = voxels checked

---

## Best Practices

### For Publication

1. **Use Adequate Sampling**
   - Global PSD: ≥ 50,000 samples for robust statistics
   - Pocket-wise PSD: ≥ 2,000 samples per region
   - Validate: Compare with different random seeds

2. **Report Complete Statistics**
   - Mean, median, mode, std, variance
   - Min, max, and range
   - 5th, 25th, 50th, 75th, 95th percentiles
   - Full histogram with bin specifications

3. **Document Parameters**
   - Grid spacing used
   - Number of samples
   - Bin edges
   - Random seed (for reproducibility)
   - Software versions

4. **Include Visualizations**
   - Full PSD histogram
   - Cumulative distribution function (CDF)
   - Log-log and semi-log plots
   - Comparison with other methods

### For Materials Screening

1. **Balance Speed vs Accuracy**
   - Use spacing = 0.5 Å (coarser grid)
   - Reduce samples to 5,000
   - Process multiple materials in parallel

2. **Standardize Workflow**
   - Use same parameters across all materials
   - Create automated pipeline
   - Store results in database

3. **Quality Control**
   - Check for convergence (repeat with different seeds)
   - Compare with experimental data if available
   - Validate physical reasonability of results

---

## Future Enhancements

### Planned Features

- [ ] Adaptive sampling (more samples in complex regions)
- [ ] Parallel region processing (for pocket-wise PSD)
- [ ] GPU acceleration for distance calculations
- [ ] Automated parameter optimization
- [ ] Integration with machine learning pipelines

---

## References

1. **Original Paper**
   Gelb, L.D., & Gubbins, K.E. (1998). "Characterization of porous solids and pore-size distributions using a mean-field method." *Langmuir*, 14(8), 1887-1902.

2. **Related Methods**
   - D. Nicholson & N. Parsonage (1982). Computer simulations and the structure of microporous carbons.
   - M. Thommes et al. (1995). Determination of pore size distributions.
   - S. Bhattacharya & D. Myers (2006). Statistical mechanics of pore networks.

3. **MOFography Documentation**
   - [Pore Size Distributions Tutorial](pore_size_distributions_gelb_gubbins.md)
   - [Rust Backend Documentation](rust_backend_documentation.md)
   - [Periodic Segmentation](outputs_from_periodic_segmentation.md)

---

**Last Updated**: December 2024
**Version**: 1.0
**Maintained by**: MOFography Development Team

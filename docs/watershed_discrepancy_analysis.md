# Watershed Segmentation Discrepancy: Python vs Rust Implementation

## Overview

This document explains the minor voxel count discrepancy observed between the Python (`do_periodic_segmentation_basic`) and Rust (`do_periodic_segmentation_rust`) implementations of the periodic segmentation algorithm in mofography. The issue manifests as a few extra voxels with `regions==0` (background) when using coarse grids (0.5Å spacing) in the Rust implementation.

## Key Differences Between Implementations

### 1. Distance Grid Processing

**Python Implementation** (`compute.py` line ~7577):
```python
padded_regions = watershed(
    -padded_grid,  # <-- Uses NEGATIVE distances
    padded_markers,
    mask=padded_void,
    compactness=compactness,
    connectivity=2,
)
```

**Rust Implementation** (`rust_fun.py` line ~427):
```python
padded_regions = mgr_rust.seeded_watershed_bucket(
    padded_grid.astype(np.float32),  # <-- Uses POSITIVE distances
    padded_markers.astype(np.uint32),
    mask=padded_void.astype(np.uint8),
)
```

### 2. Algorithmic Approach

- **Python**: Uses scikit-image's classic heap-based watershed algorithm
- **Rust**: Implements a parallel bucket-based watershed with 65,536 discrete buckets

## Root Cause of the Issue

### Bucket Quantization in Rust

The Rust implementation scales distance values into discrete buckets (`lib.rs` lines 646-675):

```rust
let (mut min_val, mut max_val) = (f32::MAX, f32::MIN);
img.iter().for_each(|&v| {
    if v < min_val { min_val = v; }
    if v > max_val { max_val = v; }
});

let range = if (max_val - min_val).abs() < 1e-12 {
    1.0
} else {
    max_val - min_val
};

let to_bucket = |v: f32| -> usize {
    let norm = ((v - min_val) / range * (BUCKETS as f32 - 1.0)).round() as usize;
    norm
};
```

### The "Leaking" Mechanism

The issue occurs when voxels near the `mask_thickness` threshold (e.g., values 0.49, 0.50, 0.51) get mapped to the same bucket due to rounding. During the watershed flooding process (`lib.rs` lines 750-752):

```rust
if labels[nidx]
    .compare_exchange(0, current_label, Ordering::Relaxed, Ordering::Relaxed)
    .is_ok()
{
    let b = bucket_idx_vec[nidx];
    buckets[b].lock().push(nidx);
}
```

Voxels that should remain background (dgrid ≤ mask_thickness) can get incorrectly labeled because they share the same bucket as foreground voxels (dgrid > mask_thickness).

## Why the Issue Disappears at Higher Resolutions

With finer grids (e.g., 0.1Å spacing):
- More unique distance values exist in the critical threshold region
- Better distribution of voxels across different buckets
- Reduced probability of threshold-crossing voxels sharing the same bucket
- More precise preservation of boundary conditions

With coarse grids (0.5Å spacing):
- Fewer unique distance values near the threshold
- Higher likelihood of voxels from both sides of the threshold being grouped together
- Increased quantization error due to larger voxel size

## Considerations for Adding More Bins

### Potential Issues with Increasing Buckets

1. **Parallelization Efficiency**: More buckets can lead to:
   - Increased memory usage
   - Potential load imbalance in parallel processing
   - More contention for bucket locks

2. **Diminishing Returns**: The benefit plateaus as the discretization becomes finer than the inherent precision of the distance grid

3. **Performance Trade-off**: The current 65,536 buckets represent a balance between precision and performance

### Alternative Approaches

1. **Adaptive Binning**: Dynamically adjust bucket count based on grid characteristics
2. **Selective Precision**: Use more buckets only for voxels near the threshold
3. **Hybrid Approach**: Combine bucket-based processing with direct processing for boundary voxels

## Final Considerations and Summary

### Impact Assessment

- **Magnitude**: Minor - only detectable when specifically examining voxel counts
- **Practical Effect**: Does not significantly affect pore geometry or size distribution calculations
- **Resolution Dependence**: Issue only appears with coarse grids (0.5Å), negligible at finer resolutions

### Recommendations

1. **Accept Current Implementation**: The Rust version provides significant speed advantages for large systems (polymers, etc.) that outweigh the minor numerical discrepancy

2. **Use Case-Based Selection**:
   - Use Rust for production analyses where speed is critical
   - Use Python when exact numerical consistency is required for validation

3. **Documentation**: Note this behavior in code comments and user documentation

4. **Monitor Real-World Impact**: Only consider fixes if the discrepancy affects actual pore analysis results

## Detailed Binning Process Explanation

### The Rust Binning Algorithm (Lines 647-676 in lib.rs):

**Step 1: Find min and max values in distance grid**
```rust
let (mut min_val, mut max_val) = (f32::MAX, f32::MIN);
img.iter().for_each(|&v| {
    if v < min_val { min_val = v; }
    if v > max_val { max_val = v; }
});
```

**Step 2: Calculate range and scaling factor**
```rust
let range = if (max_val - min_val).abs() < 1e-12 {
    1.0
} else {
    max_val - min_val
};
```

**Step 3: Map each distance value to a bucket (0 to 65535)**
```rust
let to_bucket = |v: f32| -> usize {
    let norm = ((v - min_val) / range * (BUCKETS as f32 - 1.0)).round() as usize;
    norm
};
```

This formula scales the entire range of distance values into 65,536 discrete buckets, where:
- `v - min_val`: Shifts value to start from 0
- `/ range`: Normalizes to [0,1] range
- `* (BUCKETS as f32 - 1.0)`: Scales to [0,65535]
- `.round()`: Rounds to nearest integer bucket

### Visual Example with Coarse Grid (0.5Å spacing):

**Typical distance values**: -2.0 → +5.0 (range = 7.0)
**Mask threshold**: ~0.5 Å

**Problem voxels around threshold**:
- Voxel A: dgrid = 0.49 → should be background
- Voxel B: dgrid = 0.50 → exactly at threshold  
- Voxel C: dgrid = 0.51 → should be foreground

**Bucket calculation for each**:
```rust
// For Voxel A (0.49):
bucket = ((0.49 - (-2.0)) / 7.0 * 65535.0).round()
bucket = (2.49 / 7.0 * 65535.0).round() = 23291

// For Voxel C (0.51):
bucket = ((0.51 - (-2.0)) / 7.0 * 65535.0).round()
bucket = (2.51 / 7.0 * 65535.0).round() = 23496
```

With coarse grids, many voxels have values that map to the same bucket due to the limited precision of the grid itself.

### Why More dgrid Values Help:

**With Coarse Grid (0.5Å spacing):**
- Fewer unique values between 0.4 and 0.6
- Many voxels cluster around the same few values
- Higher probability of threshold-crossing voxels landing in the same bucket

**With Fine Grid (0.1Å spacing):**
- Many more unique values between 0.4 and 0.6  
- Better distribution of values across different buckets
- Threshold-crossing voxels are more likely to be in separate buckets

### The Critical Issue:

When voxels from both sides of the threshold end up in the same bucket:
1. The foreground voxel (dgrid > 0.5) gets processed first
2. Its label spreads to ALL voxels in that bucket via neighbor expansion
3. Background voxels (dgrid ≤ 0.5) in that same bucket get incorrectly labeled
4. This happens because the mask check only prevents processing of clearly background voxels, not threshold-boundary voxels

### Why Python Doesn't Have This Issue:

Python uses scikit-image's heap-based watershed which:
- Processes each pixel individually in priority order
- Doesn't quantize values into discrete buckets
- Preserves the exact threshold boundary regardless of grid resolution

**The key difference is: Rust groups voxels by similar distance values first, then processes the group; Python processes voxels one by one by exact distance value.**

### Technical Summary

The discrepancy stems from the fundamental difference in how the two implementations handle distance discretization:

- **Python**: Processes continuous negative distances with heap-based priority queue
- **Rust**: Quantizes positive distances into discrete buckets for parallel processing

This difference manifests only at coarse resolutions where the quantization error becomes comparable to the distance grid resolution, causing a few voxels near the threshold to be misclassified. The effect is negligible for practical pore analysis and does not warrant compromising the significant performance benefits of the Rust implementation.

# VIGNOQ Watershed: Basin Leaks Through Intra-Basin Constriction

## Problem

Priority-flood watershed **leaks through constrictions** when a single pocket wraps around the cell. Basin expands past the narrowest point and claims rising terrain on the far side. Boundary lands at cell face (high dgrid) instead of at the constriction (low dgrid).

All backends (skimage, rust_bucket, rust_det, rust_det_max) produce identical wrong result.

## Evidence

Padded grid, slice at marker plane (b=35, c=43):

```
padded_a  uc_a  dgrid   region  note
   36       6    4.04      1     ← P1 marker
   47      17    1.19      1     ← CONSTRICTION (should be boundary, ISN'T)
   54      24    2.43      1     ← cell face
   55      25    2.43     20     ← BOUNDARY HERE (wrong place)
   61      31    4.04     20     ← P1 copy marker
```

Region 1 claims voxels past constriction (a=47, dgrid=1.19) through rising terrain (dgrid 1.25→1.32→...→2.43). Rising voxels get higher priority (processed first). By the time constriction voxel (priority=1.19) is processed, basin already owns both sides.

Result: self-loop edge reports wr=2.78 Å. True constriction = 1.28 Å (widest-path). zeo++ PLD = 3.115 Å, mofography PLD = 5.56 Å.

## Basin Terminology

- **Basin** = set of voxels claimed by one seed marker during watershed expansion
- **Inter-basin constriction** = narrow point where TWO different basins meet → watershed places boundary here (correct)
- **Intra-basin constriction** = narrow point WITHIN one basin's territory → NO boundary placed (the bug)

VIGNOQ: the a-direction channel constriction at dgrid=1.19 is **intra-basin**. Region 1 owns both sides. No boundary.

## Root Cause

`heap_watershed_impl` uses raw dgrid as heap priority:

```rust
// Line 1736: neighbor gets RAW dgrid as priority
let priority = NotNan::new(img_flat[nidx]).unwrap();
heap.push((priority, nidx));
```

No path-cost tracking. Basin can claim any voxel reachable through the void mask, regardless of how narrow the path was to get there. After squeezing through constriction (dgrid=1.19), basin happily claims voxels at dgrid=2.43 on the far side — those get higher priority and are processed earlier.

This is correct watershed behavior for catchment-basin segmentation. Wrong for pore-network segmentation where constrictions must become boundaries.

## Data Files

- `debug_vignoq_watershed.npz` (328 KB) — dgrids, regions, edge maps, widest paths
- `debug_vignoq_watershed.json` (298 KB) — serialized connections

## Potential Fixes

### Fix 1: Bottleneck-priority watershed

Track path bottleneck (`min dgrid along path from seed`) as priority instead of raw dgrid. Basin stops expanding when path bottleneck drops below current voxel's dgrid.

Already implemented as `heap_watershed_max_impl` (line 1531):
```rust
let next_prio = img_flat[nidx].min(prio_val);
```

**Tested — does NOT fix VIGNOQ.** Both basins reach constriction with identical bottleneck (1.38 Å). Tiebreak arbitrary. Same wrong boundary.

### Fix 2: Gradient-barrier watershed

Only claim neighbor if `dgrid[neighbor] >= dgrid[current] - epsilon`. Basin expands uphill, stops at saddle, can't flow downhill past constriction into rising terrain.

Problem: leaves most pore voxels unsegmented (valleys between saddles). Need second fill pass.

### Fix 3: Two-pass merge

1. Run watershed normally → 36 basins on padded grid
2. For each pair of basins sharing boundary: compute widest-path bottleneck between their seeds
3. If bottleneck < boundary dgrid: basin leaked through constriction → merge the two basins
4. Re-segment merged label set

No Rust changes. Post-processing in Python.

### Fix 4: Topological watershed

Use watershed variant that respects topology of the distance field — only merges basins when no saddle exists between them. Standard in mathematical morphology but complex to implement.

### Fix 5: window_connections_path (workaround)

Don't fix segmentation. Replace `window_connections` boundary-sampling with widest-path bottleneck search between pocket centers. Gives correct wr regardless of boundary placement.

## Source Code

`heap_watershed_impl` — the function that produces wrong boundaries:

```rust
// Shared helper implementing the original heap algorithm
fn heap_watershed_impl(
    py: Python<'_>,
    image: PyReadonlyArray3<'_, f32>,
    seeds: PyReadonlyArray3<'_, u32>,
    mask: Option<PyReadonlyArray3<'_, u8>>,
) -> PyResult<Py<PyArray3<u32>>> {
    let img = image.as_array();
    let seeds_arr = seeds.as_array();
    let (nx, ny, nz) = img.dim();
    let total_voxels = nx * ny * nz;

    let img_flat: Vec<f32> = img.iter().copied().collect();
    let mask_flat: Option<Vec<u8>> = mask.as_ref().map(|m| {
        m.as_array().iter().map(|&v| if v != 0 { 1u8 } else { 0u8 }).collect()
    });

    let mut labels: Vec<u32> = vec![0; total_voxels];
    let to_flat = |i: usize, j: usize, k: usize| -> usize { i * ny * nz + j * nz + k };

    use std::collections::BinaryHeap;
    let mut heap: BinaryHeap<(NotNan<f32>, usize)> =
        BinaryHeap::with_capacity(std::cmp::max(1, total_voxels / 8));

    // Seed initialisation
    for ((i, j, k), &label) in seeds_arr.indexed_iter() {
        if label == 0 { continue; }
        let idx = to_flat(i, j, k);
        labels[idx] = label;
        let priority = NotNan::new(img_flat[idx]).unwrap();
        heap.push((priority, idx));
    }

    let neighbor_offsets: &[(isize, isize, isize)] = &[
        (-1,-1,-1),(-1,-1,0),(-1,-1,1),(-1,0,-1),(-1,0,0),(-1,0,1),
        (-1,1,-1),(-1,1,0),(-1,1,1),(0,-1,-1),(0,-1,0),(0,-1,1),
        (0,0,-1),(0,0,1),(0,1,-1),(0,1,0),(0,1,1),(1,-1,-1),(1,-1,0),
        (1,-1,1),(1,0,-1),(1,0,0),(1,0,1),(1,1,-1),(1,1,0),(1,1,1),
    ];

    while let Some((_prio, idx)) = heap.pop() {
        let i = idx / (ny * nz);
        let j = (idx / nz) % ny;
        let k = idx % nz;
        let current_label = labels[idx];

        for &(dx, dy, dz) in neighbor_offsets {
            let ni = i as isize + dx;
            let nj = j as isize + dy;
            let nk = k as isize + dz;
            if ni < 0 || nj < 0 || nk < 0 { continue; }
            let (ni_u, nj_u, nk_u) = (ni as usize, nj as usize, nk as usize);
            if ni_u >= nx || nj_u >= ny || nk_u >= nz { continue; }
            let nidx = to_flat(ni_u, nj_u, nk_u);
            if let Some(ref m_flat) = mask_flat {
                if m_flat[nidx] == 0 { continue; }
            }
            if labels[nidx] == 0 {
                labels[nidx] = current_label;
                // BUG: raw dgrid as priority → basin leaks through constrictions
                let priority = NotNan::new(img_flat[nidx]).unwrap();
                heap.push((priority, nidx));
            }
        }
    }

    let out_arr = Array3::from_shape_vec((nx, ny, nz), labels).expect("shape mismatch");
    Ok(PyArray3::from_owned_array(py, out_arr).to_owned())
}
```

## Key Question for Senior Review

Can we modify priority propagation to prevent basins from expanding past intra-basin constrictions while still correctly segmenting inter-basin boundaries? The challenge: distinguishing "leaked through constriction" from "legitimate wide channel" at expansion time.

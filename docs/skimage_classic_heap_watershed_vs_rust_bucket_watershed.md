# skimage_classic_heap_watershed_vs_rust_bucket_watershed.md

## Watershed Segmentation: Classic Heap vs. Rust Bucket Implementation

---

## Introduction

Watershed segmentation is a standard method for dividing an image or 3D grid into contiguous regions using user-provided "seeds" or "markers". In our projects, it's common to run this on large, padded distance fields to segment structural features like cavities or channels. Here we explain:

- **Classic Heap-Based Watershed** (as in Python's `scikit-image`)
- **Rust Bucket-Based Watershed** (our high-performance parallel implementation)

We'll detail both algorithms, walk through the labeling process from initial state to flooding (with lots of step-by-step clarity!), and compare them for speed, accuracy, parallelization, and memory.

---

## Classic Heap-Based Watershed (scikit-image/Python)

### Overview
The Classic Heap-Based Watershed algorithm uses a priority queue (min-heap) to manage the order in which pixels (or voxels) are processed. This method is akin to filling a landscape with water, where the water flows from the lowest points (valleys) to the highest points (peaks).

### Detailed Explanation
1. **Initial Setup:**
   - **Input:** A 2D or 3D grid representing an image or distance field, a marker array where each known feature is assigned a unique label, and a mask to differentiate between foreground and background.
   - **Markers:** These are the starting points for the flooding process, like seeds from which regions will grow.

2. **Building the Min-Heap:**
   - Only the marker voxels are inserted into the heap initially. Each element in the heap contains the voxel's value, its position, and its label.
   - The heap ensures that the voxel with the smallest value (lowest point) is processed first.

3. **Flooding Process:**
   - The algorithm repeatedly pops the voxel with the lowest value from the heap.
   - For each voxel, its neighbors are examined. If a neighbor is within the mask and unlabeled, it is assigned the same label as the current voxel and added to the heap.
   - This process continues until the heap is empty, meaning all reachable voxels have been labeled.

4. **Boundary Formation:**
   - As the flooding progresses, boundaries naturally form where different regions meet. These boundaries are determined by the order in which voxels are processed.

### Analogy
Imagine a landscape where water flows from the lowest valleys to the highest peaks. The heap ensures that water always flows from the lowest point available, gradually filling up the landscape.

## Rust Bucket-Based Watershed

### Overview
The Rust Bucket-Based Watershed algorithm uses a bucket-based approach to process voxels in parallel. This method is like filling a landscape with water in layers, where each layer represents a range of values.

### Detailed Explanation
1. **Initial Setup:**
   - **Input:** Similar to the heap-based method, with a grid, markers, and a mask.
   - **Markers:** Serve as the starting points for the flooding process.

2. **Bucket Creation:**
   - The grid values are discretized into buckets, where each bucket represents a specific value or a small range of values.
   - Voxels are grouped into these buckets based on their values.

3. **Flooding Process:**
   - The algorithm processes buckets from the deepest (lowest value) to the shallowest (highest value).
   - Within each bucket, active voxels attempt to claim their unlabeled neighbors. Newly claimed neighbors are enqueued into the appropriate (shallower) bucket.
   - This process is highly parallelizable, as multiple voxels within the same bucket can be processed simultaneously.

4. **Atomic Claims:**
   - When multiple regions attempt to claim the same voxel, atomic operations ensure that only one label is assigned. This prevents conflicts and ensures consistency.

5. **Boundary Formation:**
   - Boundaries form as regions expand and meet. The use of atomic operations ensures that the process is deterministic and reproducible.

### Analogy
Imagine filling a landscape with water in layers, where each layer represents a specific elevation. All points at the same elevation are filled simultaneously, allowing for efficient and parallel processing.

## Key Differences and Implications
- **Flooding Order:** The heap-based method uses a global priority queue, while the bucket-based method processes voxels level-by-level.
- **Parallelization:** The heap-based method is challenging to parallelize due to the heap bottleneck, whereas the bucket-based method is highly parallelizable.
- **Memory Behavior:** The heap can grow large at boundaries, while buckets are more cache-friendly and predictable.

By explaining these algorithms in simple terms and using clear analogies, students can better understand the differences and advantages of each method. If you have any specific areas you'd like to focus on or further clarify, feel free to let me know!

## 1. Classic Heap-Based Watershed (scikit-image/Python)

### Initial State
- Input:
  - 3D or 2D grid (distance field or intensity image).
  - Marker array: At each known maximum or feature, assign a unique integer label. Everywhere else is 0 (unlabeled).
  - Mask: Foreground voxels (value 1/True), background voxels (0/False).
- Result: All voxels inside the mask but not at marker positions are initially unassigned (label 0).

### Detailed Step-by-Step: Heap Flooding Process

1. Build the Min-Heap (priority queue)
   - Insert only marker voxels into the heap.
   - Each heap element carries: (voxel value, voxel index, marker label).
   - The heap pops the smallest (lowest) value first.

2. Front Propagation (repeat until heap empty)
   - Pop the voxel (value, pos, label) with lowest value.
   - For each neighbor of pos (4/8 in 2D, 6/26 in 3D):
     - If the neighbor is inside the mask and unlabeled:
       - Assign it the same label.
       - Push (value_at_neighbor, neighbor_pos, label) into the heap.
   - Neighbors already labeled by another region are ignored, producing boundaries.

3. Stopping condition: When the heap is empty, all reachable foreground voxels are labeled.

### Worked 1D Example: Classic Heap-based Watershed

Suppose 7 voxels in a line (1D):

| Index | Value  | Notes                    |
|-------|--------|--------------------------|
| 1     | 0      | Wall                     |
| 2     | -0.5   | Well around marker 3     |
| 3     | -1     | Local maximum (marker)   |
| 4     | -0.5   | Well around marker 3     |
| 5     | 0      | Boundary / Wall          |
| 6     | -0.25  | Another maximum (marker) |
| 7     | 0      | Wall                     |

Markers: index 3 → label 10 (region A); index 6 → label 20 (region B)

Initial heap (min value first): [(-1, 3, 10), (-0.25, 6, 20)]

Step-by-step heap operations:

1. Pop (-1, 3, 10): neighbors 2 and 4 are unlabeled → label both 10, push (-0.5,2,10) and (-0.5,4,10).
2. Pop (-0.5, 2, 10): neighbor 1 unlabeled → label 1 as 10, push (0,1,10).
3. Pop (-0.5, 4, 10): neighbor 5 unlabeled → label 5 as 10, push (0,5,10).
4. Pop (-0.25, 6, 20): neighbor 5 already labeled 10 (ignored), neighbor 7 unlabeled → label 7 as 20, push (0,7,20).
5. Pop remaining zeros: (0,1,10), (0,5,10), (0,7,20) — no new labels.

Final labels: [10, 10, 10, 10, 10, 20, 20]

Explanation: Region A (marker at 3) expands first due to lower values and claims cells 1–5; region B claims 6 and 7.

---

## 2. Rust Bucket-Based Watershed

### Initial State
- Input: same as heap version (distance grid, marker array, mask).
- Only markers are initially labeled; others are 0.

### Bucket Flooding Overview
- Values are discretized into "buckets" (one bucket per distinct value or small ranges).
- Algorithm processes buckets from deepest (lowest value) upward. For each bucket, active cells attempt to claim their unlabeled neighbors and enqueue newly-claimed neighbors into the appropriate (shallower) bucket.
- Concurrency: multiple cells at the same bucket are processed together. Claims on neighbors use atomic compare-and-swap (or equivalent) so only one label wins per voxel.

### Worked 1D Example: Bucket-based Watershed (cleaned)

We use the same 7-voxel example to illustrate the bucket method. Values and markers:

| Index | Value  | Marker |
|-------|--------|--------|
| 1     | 0      |        |
| 2     | -0.5   |        |
| 3     | -1     | 10     |
| 4     | -0.5   |        |
| 5     | 0      |        |
| 6     | -0.25  | 20     |
| 7     | 0      |        |

Bucket grouping (by exact values here):
- bucket[-1]: [3]
- bucket[-0.5]: [2, 4]
- bucket[-0.25]: [6]
- bucket[0]: [1, 5, 7]

Initial labels: [0, 0, 10, 0, 0, 20, 0]

Processing proceeds from deepest to shallowest buckets. Newly-claimed neighbors are enqueued into shallower buckets.

Pass 1 — process bucket[-1]
- Active: [3(label=10)]
- Neighbors 2 and 4 are unlabeled → atomically claim both with 10, enqueue them into bucket[-0.5].
- Labels: [0, 10, 10, 10, 0, 20, 0]

Pass 2 — process bucket[-0.5]
- Active: [2(label=10), 4(label=10)]
- Cell 2 claims neighbor 1 → enqueue 1 into bucket[0] with label 10.
- Cell 4 claims neighbor 5 → enqueue 5 into bucket[0] with label 10.
- Labels after (claims applied atomically): [10, 10, 10, 10, 0, 20, 0]

Pass 3 — process bucket[-0.25]
- Active: [6(label=20)]
- Neighbor 5: in a single-threaded deterministic run, 5 is already claimed by 10 from Pass 2, so 6 cannot claim it.
- Cell 6 claims neighbor 7 → enqueue 7 into bucket[0] with label 20.
- Labels: [10, 10, 10, 10, 10, 20, 20]

Pass 4 — process bucket[0]
- Active: [1(label=10), 5(label=10), 7(label=20)]
- All neighbors already labeled; no further claims.

Final labels: [10, 10, 10, 10, 10, 20, 20]

Notes on atomicity and determinism:
- When multiple regions attempt to claim the same voxel, implementations must use atomic writes so only one succeeds. A deterministic tie-break rule (e.g., first-writer, lower label wins, or deterministic traversal order) ensures reproducible results across runs.
- With correct atomic semantics the bucket method produces the same watershed segmentation as the heap-based method for this example.

---

## Direct Comparison: Differences and Implications

| Feature              | Classic Heap (skimage)                 | Rust Buckets                                       |
|---------------------:|:---------------------------------------|:--------------------------------------------------|
| Flooding order       | Global priority queue (exact value order) | Level-by-level (buckets) with parallel processing |
| Parallelization      | Hard to parallelize (heap bottleneck)   | Highly parallel (process many cells per bucket)   |
| Memory behaviour     | Heap can grow large at boundaries       | Buckets are cache-friendly and predictable        |
| Boundary formation   | First pop wins (by value order)          | Atomic claim wins; deterministic tie-break needed |
| Output               | Watertight labels (with correct seeds)  | Same (with atomic claims & deterministic rules)   |

Practical tip: ensure at least one marker per connected foreground region to avoid unlabeled islands.

---

## Visual/Conceptual Analogy
- Heap-Classic: like a single-file line where the absolute closest person moves next.
- Bucket-Based: like whole steps (levels) being taken simultaneously by groups at the same elevation.

---

## Further Reading
- skimage Watershed Example: <https://scikit-image.org/docs/0.25.x/auto_examples/segmentation/plot_watershed.html>
- skimage Watershed API: <https://scikit-image.org/docs/0.25.x/api/skimage.segmentation.html>

*If you'd like, I can also run a quick consistency check script or add a small ASCII diagram of the propagation order.*


# Watershed Algorithms Explained: rust_det, rust_det_parallel, and skimage

This document explains in simple terms how the deterministic watershed implementations work (`rust_det`, `rust_det_parallel`), how scikit-image’s watershed works, and why our Rust heap implementation is faster than skimage. A short 1D example illustrates the flooding order.

---

## 1. The idea: seeded watershed by “flooding”

We have:

- **Image** (e.g. a 3D distance grid `dgrid`): each voxel has a value (e.g. distance to the nearest atom).
- **Seeds** (markers): some voxels already have a label (1, 2, 3, …). Labels are the “regions.”
- **Mask** (optional): only voxels where mask is true are considered; the rest are ignored.

**Goal:** Assign every masked voxel to exactly one seed’s region, so that “water” appears to flow from seeds and fill the volume in order of **priority**. For a **distance** grid we want “higher distance first” (pore centers before walls), so the flooding order is **by decreasing distance** (or equivalently by increasing `-dgrid`).

**Rule:** When we “process” a voxel, we look at its neighbors. Any **unlabeled** neighbor gets our label and is then processed later. The **order** in which we process voxels is what makes the result deterministic and physically meaningful: we always expand from the “highest priority” (e.g. largest distance) first.

---

## 2. Simple 1D example: how flooding order works

Take a 1D “image” of 7 cells with values (e.g. distance):

```text
index:    0    1    2    3    4    5    6
value:    1    3    5    4    2    6    2
```

Seeds: **label 1** at index 0, **label 2** at index 6. We flood in order of **decreasing value** (process the 6 first, then 5, then 4, …).

**Initial:** labels = `[1, 0, 0, 0, 0, 0, 2]`. Queue (heap or bucket frontier) contains the two seed indices with their values.

**Flooding steps (max-heap, pop largest first):**

1. Heap has `(1, 0)` and `(2, 6)`. Pop **max** → `(2, 6)` (value 2 at index 6). Process index 6: only neighbor 5 is unlabeled → assign label 2 to index 5, push `(6, 5)`.  
   Labels = `[1, 0, 0, 0, 0, 2, 2]`.

2. Pop `(6, 5)`. Process index 5 (value 6): neighbors 4 and 6; 4 unlabeled → assign 2, push `(4, 4)`.  
   Labels = `[1, 0, 0, 0, 2, 2, 2]`.

3. Pop `(4, 4)`. Process index 4 (value 4): neighbor 3 unlabeled → assign 2, push `(4, 3)`.  
   Labels = `[1, 0, 0, 2, 2, 2, 2]`.

4. Meanwhile we also pop `(1, 0)` and expand label 1: assign 1 to index 1 (value 3), then from index 1 assign 1 to index 2 (value 5). So index 2 gets label 1 before label 2 can reach it from the right, because we process by **value**: value 5 at index 2 is processed when it enters the heap from the left.

5. **Tie-breaking:** If two neighbors had the same value and both tried to claim the same cell, we’d break ties (e.g. by label or index) so the result is deterministic.

**Takeaway:** We store `(value, index)` and always pop the **largest** value. When we assign a neighbor, we push `(neighbor_value, neighbor_index)`. So we expand in order of decreasing value; **the first label to reach a cell wins**. In the bucket version we process the same order (high value first) using buckets instead of a heap; within a bucket we resolve conflicts with a deterministic sort. The 1D picture is the same in 3D: we still flood in priority order; only the number of neighbors (e.g. 26) and the tie-breaking details change.

---

## 3. rust_det: heap-based deterministic watershed

**Where it lives:** `mgr_rust.seeded_watershed_det(image, seeds, mask)` (Rust: heap-based, single-threaded).

**Idea:** Use a **max-heap** so we always process the voxel with the **largest** image value (highest distance) first. Same semantics as “flood from high to low” in a distance field.

**Steps:**

1. **Input:** `image` = distance grid (e.g. `dgrid`), `seeds` = label array (0 = no seed), optional `mask`.
2. **Initialize:** For every voxel where `seeds[i,j,k] != 0`, set `labels[i,j,k] = seeds[i,j,k]` and push `(image[i,j,k], index)` into a **binary max-heap** (priority = image value).
3. **Loop:** While the heap is not empty:
   - Pop the voxel `(priority, idx)` with the **largest** priority.
   - For each 26-neighbor of `idx` that is in the mask and has `labels[neighbor] == 0`:
     - Set `labels[neighbor] = labels[idx]`.
     - Push `(image[neighbor], neighbor)` onto the heap.
4. **Output:** `labels`.

**Determinism:** The heap breaks ties in an implementation-defined way (e.g. by index). So for a given image and seeds, the result is always the same. **First label to reach a voxel wins** (and that order is fixed by the priority order).

**Why “det”:** The result is deterministic and matches the classic “flood in order of decreasing distance” semantics used for region segmentation.

---

## 4. rust_det_parallel: bucket-based deterministic watershed

**Where it lives:** `mgr_rust.seeded_watershed_det_parallel(image, seeds, mask)` (Rust: bucket-based, parallel inside each level).

**Idea:** Get the **same** flooding order as the heap (high distance first) without a single global heap. We bucketize by **priority** and process bucket-by-bucket; inside each bucket we run in parallel, then resolve conflicts with a **deterministic** sort and sequential commit.

**Steps:**

1. **Bucketize:** Compute `f = -image` (so higher image value ⇒ smaller `f`). Map `f` into 65536 buckets in order (bucket 0 = smallest `f`, bucket 65535 = largest `f`). So **bucket index increases as image value increases**. We process buckets in **ascending** order (bucket 0 first, then 1, …), which is **from smallest f to largest f**, i.e. **from highest image value to lowest** — same as the heap.
2. **Initialize:** For each seed, set its label and add its index to the bucket corresponding to its image value.
3. **For each bucket level (0 to 65535):**
   - **Frontier:** Take the current list of voxel indices in this bucket (the “frontier”).
   - **Propose in parallel:** Each frontier voxel looks at its 26 neighbors. For every **unlabeled** neighbor it proposes `(neighbor_index, my_label, bucket_for_neighbor)`.
   - **Merge and sort:** Collect all proposals; sort by `(neighbor_index, label, bucket)` so the order is deterministic.
   - **Commit sequentially:** For each `(nidx, lab, b)` in order, if `labels[nidx] == 0` then set `labels[nidx] = lab` and append `nidx` to bucket `b`. So the “first” proposal in sort order wins — same idea as “first to reach wins.”
4. Repeat until no bucket has new work (inner loop per level until the current bucket’s list is empty), then move to the next bucket.

**Determinism:** The sort key `(index, label, bucket)` fixes a unique order. So the same input always gives the same output, and the semantics match the heap (same flooding order and tie-breaking philosophy).

**Why faster on big grids:** Within a bucket we only do parallel neighbor traversal and local buffers; the expensive part (heap) is replaced by cheap buckets and one sort per frontier. So we get parallelism without a shared heap bottleneck.

---

## 5. How skimage’s watershed works

**Where it lives:** `skimage.segmentation.watershed(image, markers, mask=..., connectivity=..., compactness=0)`.

**Convention:** skimage floods from **low** image values to **high**. So “valleys” get filled first. To get “high distance first” (pores before walls), we pass **negative** distance: `watershed(-dgrid, markers, mask=...)`.

**Mechanism (conceptually):**

- Uses a **priority queue** (min-heap): smaller image value = higher priority.
- Initially, only marker (seed) pixels are in the heap.
- Repeatedly pop the pixel with the **smallest** value; for each unlabeled neighbor in the mask, assign the same label and push the neighbor with **its** image value.
- So the flooding order is “increasing image value.” With `-dgrid`, that’s “decreasing dgrid” ⇒ high distance first, same as our rust_det.

So **algorithmically**, skimage’s seeded watershed is the same idea as rust_det: one priority queue, process in priority order, first label to reach a pixel wins. The only difference is the sign of the image (skimage: low first; we pass `-dgrid`. rust_det: high first on `dgrid`).

---

## 6. Comparison summary

| Aspect              | rust_det           | rust_det_parallel        | skimage                    |
|--------------------|--------------------|---------------------------|----------------------------|
| Structure          | Single max-heap    | Buckets + sort per level  | Single min-heap (on image) |
| Image sign         | Use +dgrid (high first) | Use +dgrid, bucketize -dgrid | Use -dgrid (low first)     |
| Threading          | Single-threaded    | Parallel per bucket       | Single-threaded            |
| Deterministic      | Yes                | Yes (sort + sequential commit) | Yes (heap order)           |
| Same semantics     | Yes                | Yes (by design)           | Yes (with -dgrid)          |

So all three can produce the same labeling for the same `dgrid`, seeds, and mask when used as intended (e.g. `skimage.watershed(-dgrid, markers, mask)` vs `rust_det(dgrid, seeds, mask)`).

---

## 7. Why is rust_det faster than skimage?

- **Language and runtime:** rust_det is compiled Rust; skimage is Python calling compiled C/Cython for the core. Even if the core of skimage uses a heap in C, there is still Python overhead (array handling, type checks, function call overhead) and possibly extra copies (e.g. converting to the format the C code expects).
- **Data layout and hot loop:** Rust uses flat `Vec`s and a single heap of `(priority, index)`; the inner loop is tight (neighbor offsets, index math, no Python interpreter). skimage may use more generic ndarray indexing and object handling in the hot path.
- **Single heap, cache-friendly:** One binary heap and contiguous label/image arrays are very cache-friendly in Rust. In Python, even with NumPy, the watershed implementation may do more indirection and allocation.
- **No GIL in the core:** The entire watershed loop in rust_det runs in native code without the Python GIL, so there’s no interpreter or GC in the middle of the loop.

So **rust_det is faster than skimage** mainly because it’s a tight, compiled, single-threaded implementation with minimal indirection and no Python overhead in the inner loop. **rust_det_parallel** can be faster still on large grids by spreading work across cores while preserving the same deterministic result as rust_det.

---

## 8. Summary

- **Flooding:** We assign labels by “flooding” from seeds in a fixed **priority order** (e.g. by decreasing distance). The 1D example shows that the structure that decides “who is processed next” (heap or buckets) controls who wins at boundaries.
- **rust_det:** Max-heap on image value; process highest first; first label to reach a voxel wins; deterministic; single-threaded; fast because it’s tight Rust.
- **rust_det_parallel:** Same order as the heap, but with buckets and per-level parallel propose + deterministic sort + sequential commit; same semantics, deterministic, parallel.
- **skimage:** Min-heap on image value (so use `-dgrid` for “high distance first”); same idea as rust_det; slower due to Python and implementation details.
- **Why rust_det beats skimage:** Compiled Rust, minimal overhead, cache-friendly hot loop, no Python in the core.

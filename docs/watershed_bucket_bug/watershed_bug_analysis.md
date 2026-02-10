# Watershed background-voxel bug: analysis and potential fixes

## Summary

Void voxels (dgrid > mask_thickness) sometimes remain labeled 0 after `seeded_watershed_bucket`, so they are never assigned to any region. The evidence is that many of these voxels have at least one labeled neighbor—i.e. they sit inside a region but were never flooded.

**Conclusion:** The bug is in the bucket-based watershed in `lib.rs` (`seeded_watershed_bucket`). Maxima detection and `do_periodic_segmentation_rust` in `rust_fun.py` are only relevant insofar as they feed the correct image, seeds, and mask into the watershed; they do not cause the missing labels.

---

## Simple explanation (plain-language summary)

### How does the bucket watershed work?

Instead of processing voxels one-by-one by value (like a heap), the algorithm groups voxels into **buckets** by their distance value. There are 65536 buckets. Processing goes from **highest value down to lowest** (like filling from the pore centers outward).

For each bucket level:

1. Take all voxels currently in that bucket (the "frontier").
2. In parallel, each frontier voxel tries to "claim" its unlabeled neighbors (26-connectivity). If a neighbor is in the mask and still has label 0, one region wins via atomic compare-and-swap.
3. Newly claimed neighbors get pushed into **their own** bucket (based on their distance value).
4. Move on to the next (lower) bucket. We never revisit a bucket.

**Relevant doc:** `docs/skimage_classic_heap_watershed_vs_rust_bucket_watershed.md` — "Rust Bucket-Based Watershed" (Overview, Bucket Creation, Flooding Process) and the worked 1D example (Pass 1–4).

### What is causing the error?

When we process bucket L, we drain it once. During the parallel loop we may **add** newly labeled voxels to bucket L (when a neighbor has the same distance value). Those voxels are added **after** we've already drained bucket L, so we never process them. They get a label but never get to propagate to **their** neighbors.

Any void voxel that can only be reached through such a "same-bucket" voxel stays 0.

**Relevant doc:** Same file — the bucket method assumes each bucket is processed once; the bug is that same-bucket voxels added mid-round are never processed.

### Will the error occur always or is it a special case?

It occurs whenever the pore geometry has **flat or near-flat regions** — i.e. many adjacent void voxels with the same (or very similar) distance value. They map to the same bucket. In curved, smoothly varying pores the effect is smaller. Pseudomaterials and MOFs with flat channels or plateaus are more likely to hit it.

### What can we do to fix it?

**Option 1 (recommended):** Re-process each bucket until quiescent. For each level L, keep draining and processing `buckets[L]` until a pass adds no new voxels to that bucket. Then move to L-1. Same-bucket chains get fully flooded.

**Option 2:** Use the heap-based watershed (`seeded_watershed` / `heap_watershed_impl` in lib.rs) instead — it processes every voxel when it’s popped, so no same-bucket omission. Slower but correct.

**Option 3:** Post-process: after the bucket watershed, assign any remaining void voxels with label 0 the label of a labeled neighbor. Band-aid, not a proper fix.

**Relevant doc:** `docs/skimage_classic_heap_watershed_vs_rust_bucket_watershed.md` — "Direct Comparison" table: heap uses exact value order; buckets use level-by-level. The fix aligns bucket behavior with that.

---

---

## 1. Data flow (rust_fun.py → lib.rs)

### 1.1 `do_periodic_segmentation_rust` (rust_fun.py)

- **Maxima:** `detect_maxima_selective_smoothing_parallel` is called on `dgrid` (no padding yet). It returns merged maxima in the unpadded grid.
- **Distance image:** `dt = dgrid.copy()` (no smoothing of the grid used for watershed; smoothing only affects maxima).
- **Mask:** `void_mask = dgrid > mask_thickness`, then `padded_void = np.pad(void_mask, ...)`.
- **Markers:** `markers` is built in the unpadded grid from merged maxima (labels 1..N), then `padded_markers = np.pad(markers, ...)`.
- **Watershed input:** `padded_grid = np.pad(dt, ...)` (distance transform), `padded_markers`, `padded_void`. So the Rust watershed sees a padded distance field, padded seeds, and a padded void mask. All are consistent.

So the Python side does not introduce inconsistent masks or seeds. The only way a void voxel stays 0 is if the watershed never assigns it a label.

### 1.2 Maxima / selective smoothing (lib.rs)

- `detect_maxima_selective_smoothing_parallel` builds a smoothed grid, runs h-maxima, finds maxima, merges them. It does not define the watershed mask or the watershed image.
- If a pore had no maximum, that whole pore would stay 0, but then we would not see “background voxels with labeled neighbors.” The fact that we do see interior background voxels means seeds exist and flooding reaches their neighbors; the bug is that some voxels that should be reached are not.

So the bug is not “missing seeds” or wrong smoothing; it is in the propagation logic of the watershed.

---

## 2. Root cause: `seeded_watershed_bucket` (lib.rs)

### 2.1 Algorithm

- Image values are mapped to `BUCKETS` (65536) levels via min–max normalization: `to_bucket(v)`.
- Seeds are stored with their label and their index is pushed to `buckets[b]` with `b = bucket_idx_vec[idx]`.
- Main loop: for `level` from `BUCKETS-1` down to 0:
  1. **Drain** `buckets[level]` into `frontier` (and leave that bucket empty).
  2. In parallel, for each index in `frontier`: for each neighbor in the mask with label 0, do `compare_exchange(0, current_label)`; on success, push the **neighbor’s** index to `buckets[bucket_idx_vec[neighbor]]`.

So propagation is “from high image value (high bucket) downward.” A voxel gets labeled when a neighbor that already has a label is processed and wins the compare_exchange for that voxel. The newly labeled voxel is then pushed to the bucket corresponding to **its own** image value.

### 2.2 Same-bucket bug

When we process level `L`:

- We only process indices that were in `buckets[L]` at the **start** of the iteration (the drained `frontier`).
- During the parallel loop we may **add** new indices to `buckets[L]`: when we label a neighbor whose image value also maps to bucket `L`, we push that neighbor’s index to `buckets[L]`.
- After the loop we move to `level = L-1` and never process `buckets[L]` again.

So any voxel that is **added** to `buckets[L]` in this iteration is **never** used as a frontier node. It has been assigned a label (correctly), but it never gets to propagate that label to **its** neighbors. Therefore:

- Any void voxel that can only be reached by propagation through such a “same-bucket” voxel will never be labeled and stays 0.
- This matches the observed bug: interior background voxels (void, label 0, with at least one labeled neighbor). They sit next to voxels that were labeled in the same bucket round and never got to propagate.

So the root cause is: **same-bucket propagation**. The algorithm assumes each bucket is processed once; but newly labeled voxels in the same bucket are added after the bucket has already been drained, so they never propagate.

### 2.3 Why it’s not the mask or seeds

- Void is defined as `dgrid > mask_thickness`; the same criterion is used for `void_mask` and then `padded_void`. So every void voxel is in the mask.
- Seeds are placed at maxima and padded; they are in the void. So the bug is not “mask excludes void” or “no seed in this pore.”
- The only way a masked void voxel stays 0 is that no labeled voxel ever won a compare_exchange for it. That happens when the only possible path goes through a same-bucket voxel that was labeled in the current round and never processed as frontier.

---

## 3. Potential fixes (no implementation)

### 3.1 Re-process each bucket until quiescent (recommended)

- For each level `L`, instead of a single drain-and-process:
  - Repeatedly: drain `buckets[L]` into `frontier`, process frontier (parallel propagate, pushing newly labeled voxels to their buckets as now), until a pass adds **no** new indices to `buckets[L]` (or until frontier is empty and nothing was added to L).
- Then advance to `L-1`.
- Effect: Every voxel that gets labeled at level L gets a chance to propagate within the same level, so same-bucket chains are fully flooded. Void voxels reachable from a seed should all get a label.

### 3.2 Inner “same-level” loop

- When processing level `L`, maintain a “current level” frontier.
- In an inner loop: process this frontier (parallel propagate); collect all indices that were written to `buckets[L]` in this pass (e.g. by having the propagation code record “added to L”). Use that set as the new frontier for level L. Repeat until no new indices are added to `buckets[L]`.
- Then drain any remaining contents of `buckets[L]` into the next level’s initial frontier if desired, or leave them for a final pass. Conceptually similar to 3.1 but with explicit tracking of “same-level” additions.

### 3.3 Heap-based watershed (reference behavior)

- Replace the bucket scheme with a priority queue (e.g. binary heap) keyed by image value (high first). Repeatedly pop the highest-value voxel, propagate to neighbors (push them with their value). No notion of “level” that we leave behind; every labeled voxel is eventually processed.
- This matches the classic seeded watershed and avoids the same-bucket bug entirely, at the cost of more allocations and possibly different performance. The existing `seeded_watershed` / `heap_watershed_impl` in lib.rs may already implement this; the bucket version was likely introduced for speed.

### 3.4 Post-process fill (mitigation only)

- After the current bucket watershed, do a pass over all void voxels with label 0 and assign them the label of an arbitrary labeled neighbor (or the majority label among neighbors). This would fix the symptom and leave the bucket algorithm as-is, but it is a band-aid and could introduce arbitrary choices at boundaries; not a substitute for fixing propagation.

---

## 4. What not to change

- **Maxima / selective smoothing:** They determine seed and image; they do not cause interior background voxels. Changing them would not fix this bug.
- **Mask or padding in rust_fun.py:** They are consistent (void = dgrid > mask_thickness, same for mask and for defining “void” when counting the bug). No change needed there for this bug.
- **compare_exchange / parallelism:** The race is “who gets to label an unlabeled voxel”; the first writer wins, which is correct. The bug is not a race on a single voxel but the **omission** of whole voxels from ever being in the frontier (same-bucket effect).

---

## 5. Verification after a fix

- Re-run the diagnostic script (`confirm_background_voxels.py`) with the same parameters (and optionally mask_thickness=0.5). After the fix, “Background voxels IN void” should be 0 (or the script could assert that).
- Adding a unit test that runs a small synthetic grid (e.g. a single pore with a flat plateau) through `seeded_watershed_bucket` and checks that all masked voxels receive a non-zero label would guard against regressions.

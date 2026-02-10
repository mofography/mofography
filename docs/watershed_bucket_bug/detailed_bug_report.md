# Watershed bucket bug: report

## Symptom

After `do_periodic_segmentation_rust`, some void voxels (`dgrid > mask_thickness`) remained with label 0 instead of receiving a region label. Those voxels are inside the pore (they can have labeled neighbors) and should have been assigned.

## Root cause

The bucket watershed processes levels **high to low**. When we are at a **low** level (e.g. neck), we can label a neighbor in a **higher** bucket (e.g. cavity) and push it to `buckets[b]` with `b > level`. Level `b` was already processed, so that voxel is never drained and never propagates. Void voxels only reachable “from below” (from a lower bucket) stayed unlabeled. This shows up in geometries with necks and cavities (common in MOFs, zeolites).

## Fix

Keep high-to-low order. When we push to `buckets[b]` with `b > level`, record `b` in a set `late_buckets`. After the main loop, **re-process** those buckets: repeatedly drain each and run the same propagation, adding any bucket that receives a push back to the set, until no bucket has pending work.

**Code:** `mgr_rust/src/lib.rs`, function `seeded_watershed_bucket`. Search for "late_buckets", "Re-process late buckets".

## Verify

From repo root:

- `python docs/watershed_bucket_bug/confirm_background_voxels.py` — should report "Background voxels IN void (BUG if > 0): 0".
- `pytest mofography/tests/test_mofography.py`

Test assets: `confirm_background_voxels.py`, `Case4-combined-revised.cif`, and `test_notes.org` (plotting) are in this folder.

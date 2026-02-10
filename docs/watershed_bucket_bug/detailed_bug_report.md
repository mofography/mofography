# Watershed bucket bug: detailed report for handoff

## Paths

- **Main repo (original source):** `/Users/arungopalan/repos/mofography`
- **Worktree (fix branch):** `/Users/arungopalan/repos/mofography-watershed-fix`
- **Rust watershed implementation:** `mgr_rust/src/lib.rs` — function `seeded_watershed_bucket` (around lines 902–1095 in the worktree).
- **Python caller:** `mofography/rust_fun.py` — `do_periodic_segmentation_rust` calls `mgr_rust.seeded_watershed_bucket` in several places (e.g. ~756, ~827, ~909).
- **Diagnostic script:** `docs/watershed_bucket_bug/confirm_background_voxels.py` — reproduces the bug with Case4 CIF and reports counts of void voxels with label 0.
- **Case4 CIF:** `docs/watershed_bucket_bug/Case4-combined-revised.cif`.
- **Analysis doc:** `docs/watershed_bucket_bug/watershed_bug_analysis.md` (initial analysis).  
- **Watershed algorithm notes:** `docs/skimage_classic_heap_watershed_vs_rust_bucket_watershed.md`.

---

## What the bug is (symptom)

After `do_periodic_segmentation_rust`, some **void voxels** (where `dgrid > mask_thickness`) remain with **label 0**. They should all receive a region label from the watershed. Observed:

- **~79,790** masked voxels still unlabeled on the **padded** grid (logged inside Rust).
- **518** background-in-void in the **unpadded** central box (diagnostic script with mask_thickness=0).
- Many of these have at least one labeled neighbor (“interior” to a region), so they are clearly inside the pore, not at the boundary.

---

## What was done so far

1. **Reproduce:** Confirmed with `confirm_background_voxels.py` and Case4 params from `docs/watershed_bucket_bug/watershed_bucket.md`.
2. **Hypothesis 1 (same-bucket):** The bucket algorithm processes each bucket once, high→low. If a voxel is labeled during the current level and pushed into the **same** bucket, it is never processed (bucket already drained). So we added an **inner loop** per level: drain `buckets[level]`, process frontier, repeat until the bucket is empty. Code change: in `seeded_watershed_bucket`, wrap the single drain+process in a `loop { ... }` so we keep re-draining and processing the same level until empty.
3. **Build/test:** `maturin build --release` in the worktree’s `mgr_rust`, then `pip install --force-reinstall target/wheels/mgr_rust-0.1.0-cp38-abi3-macosx_11_0_arm64.whl`. Running the diagnostic still showed 79,790 unlabeled and 518 in void — **no improvement**.
4. **Diagnostic logging in Rust:** Added in the worktree’s `lib.rs`:
   - On entry: `[seeded_watershed_bucket] using reprocess-same-bucket version (inner loop per level)`.
   - In the inner loop: when we push a newly labeled neighbor to `buckets[b]`, if `b == level` we increment `added_to_current_level` and log `[ws] level L pass P: frontier=N added_to_same_bucket=M`.
   - After finishing each level: count masked voxels with `bucket_idx_vec[idx]==level` and `label==0`; if > 0, log `[ws] level L DONE: N masked voxels in this bucket still unlabeled`.
   - At the end: total unprocessed masked count and summary line.
5. **Log results:** No `added_to_same_bucket` lines appeared. Many `level L DONE: N masked voxels in this bucket still unlabeled` lines (many levels with 12–48 unlabeled each).

---

## What we thought the issue was (and where)

- **Initial hypothesis:** Same-bucket propagation. When processing level L we add voxels to `buckets[L]` after draining L, so they are never processed. **Location:** single pass per level in `seeded_watershed_bucket`, main `for level in (0..BUCKETS).rev()` loop.
- **Fix applied:** Inner `loop { drain(level); if empty break; process(frontier); }` so same-bucket voxels get processed on the next pass. **Result:** Bug persists; logs show we never actually push to the current bucket (no `added_to_same_bucket`).

---

## What we think the issue is now

- **New hypothesis:** Voxels are being pushed into **already-processed (higher) buckets**. Algorithm processes levels **high → low** (65535 down to 0). When we are at a **low** level (e.g. narrow neck, low distance value), we can label a neighbor that lies in a **higher** bucket (e.g. cavity). We push that neighbor to `buckets[b]` with `b > level`. Level `b` was already processed earlier, so that voxel is never drained again and never propagates. So regions that are only reachable “from below” (from a lower bucket) stay unlabeled.
- **Evidence:** Logs show many “level L DONE: N masked voxels in this bucket still unlabeled” — those voxels were never in the frontier for level L. The only way they could get labeled later is when we process a **lower** level and label a voxel adjacent to them; that would push them to bucket L, but L is already done.
- **Proposed direction:** Either (1) process buckets **low → high** instead of high→low, or (2) keep high→low but **re-process any bucket that receives new voxels after we have passed it** (e.g. track “dirty” higher buckets and run a second pass or iterative reprocess until quiescent).

---

## What we are unsure of

- Whether **reversing** to low→high changes the intended watershed semantics (flood from peaks first) and whether that is acceptable or needs a different fix (e.g. late re-process of higher buckets only).
- Whether there are **disconnected** mask components (no seed) contributing to the counts, or if all 79,790 are explainable by the “added to already-passed bucket” effect.
- Exact bucket index convention: `to_bucket(v)` maps image value to `0..BUCKETS-1`; we iterate `(0..BUCKETS).rev()` so we process high bucket index first (high image value). Confirm that “higher bucket” = higher image value = pore center, “lower bucket” = lower value = nearer walls.

---

## Algorithm summary (for reference)

- **Constants:** `BUCKETS = 65536`. Image values scaled to bucket index via min–max.
- **Seeds:** Placed at maxima; each seed gets a label and is pushed to `buckets[bucket(seed)]`.
- **Main loop:** `for level in (0..BUCKETS).rev()` — for each level, inner `loop`: drain `buckets[level]` into `frontier`; if empty, break; else in parallel for each idx in frontier, for each 26-neighbor n in mask with label 0, `compare_exchange(0, label)`; on success push `n` to `buckets[bucket(image[n])]`.
- **Output:** Copy labels to output array.

---

## How to recompile and test

1. **Go to worktree:**  
   `cd /Users/arungopalan/repos/mofography-watershed-fix/mgr_rust`

2. **Build release wheel:**  
   `maturin build --release`  
   Wheel is written to `target/wheels/mgr_rust-0.1.0-cp38-abi3-macosx_11_0_arm64.whl` (or equivalent for your platform).

3. **Install (force reinstall so the new backend is used):**  
   `pip install --force-reinstall target/wheels/mgr_rust-0.1.0-cp38-abi3-macosx_11_0_arm64.whl`  
   (Use the exact filename printed by maturin.)

4. **Run diagnostic from main repo:**  
   `cd /Users/arungopalan/repos/mofography`  
   `python docs/watershed_bucket_bug/confirm_background_voxels.py`  
   Stderr will show the `[seeded_watershed_bucket]` and `[ws]` logs. The script prints “Background voxels IN void (BUG if > 0): N”; with a correct fix this should be 0.

5. **Optional — run tests:**  
   `pytest mofography/tests/test_mofography.py`  
   (Uses a different CIF; ensures segmentation still completes and return shapes are correct.)

---

## Fix applied (re-process late buckets)

**Root cause (confirmed):** Voxels were being pushed into **already-processed (higher) buckets**. The main loop processes levels **high to low** (`(0..BUCKETS).rev()`). When at a **low** level (e.g. neck), we can label a neighbor in a **higher** bucket (e.g. cavity) and push it to `buckets[b]` with `b > level`. Level `b` was already drained and never revisited, so that voxel never propagated; masked voxels only reachable "from below" stayed unlabeled.

**Physical/geometric conditions (not rare):** The bug is *not* the "same-bucket" case (flat surfaces, not enough propagation within one bucket). That was an earlier hypothesis; an inner loop per level already ensured same-bucket voxels are drained until empty, and the bug persisted. The real trigger is **cavity voxels that are first reached only from a neck (lower distance value)**. In that situation we process the cavity's bucket (high level) before any voxel in that cavity is in the frontier; later, when we process the neck's bucket (low level), we label cavity voxels adjacent to the neck and push them to the cavity bucket—which we have already passed—so they never propagate. This geometry (narrow windows/necks opening into larger cages) is **common** in porous materials: many MOFs, zeolites, and slit-like or channel structures with constrictions exhibit it. So the physical aspects leading to the bug are **not rare**; they are typical of cage-and-window and neck-and-cavity pore space. Case4 (and similar test cases) hit the bug because of that geometry.

**Fix:** Keep high-to-low order. Record any bucket index `b` when we push to `buckets[b]` and `b > level`. After the main loop, **re-process** those buckets: repeatedly drain and run the same propagation for each such bucket, and add any bucket that receives a push during this phase back to the "to process" set, until no bucket has pending work.

**Code references (worktree `mofography-watershed-fix`, `mgr_rust/src/lib.rs`):**

- **Late-bucket set:** Lines 1013-1015 — `late_buckets: Arc<Mutex<HashSet<usize>>>` created before the main loop.
- **Recording late pushes (main loop):** Lines 1054-1058 — inside the propagation closure, after `let b = bucket_idx_vec[nidx]`, `if b > level { late_buckets.lock().insert(b); }` then `buckets[b].lock().push(nidx)`.
- **Re-process loop:** Lines 1065-1115 — after the main `for level in (0..BUCKETS).rev()`, a `loop` drains `late_buckets` into `to_process`; for each `b` in `to_process`, an inner `loop` drains `buckets[b]` and runs the same `frontier.par_iter().for_each(...)` propagation; on each push to `b'` we do `late_buckets.lock().insert(b_prime)` so that bucket is processed again.
- **Performance:** The per-level full scan over all voxels (previously after each level) and the hot-path atomic/`eprintln!` diagnostics were removed so the 65k-level loop only does drain + propagate. A single post-loop unprocessed count is kept under `#[cfg(debug_assertions)]` (lines 1117-1129).

**Verification:** Run `docs/watershed_bucket_bug/confirm_background_voxels.py`; "Background voxels IN void (BUG if > 0)" should be 0. `pytest mofography/tests/test_mofography.py` should pass.

---

## Files to inspect in the worktree

- `mgr_rust/src/lib.rs`: function `seeded_watershed_bucket` (search for "Main flooding loop", "late_buckets", "Re-process late buckets"). The main loop (lines ~1019-1063), late-bucket recording (1055-1057), and re-process loop (1066-1115) implement the fix.

# Watershed bucket bug: detailed report for handoff

## Paths

- **Repo:** `mofography` (fix is on `master`; the former worktree has been removed.)
- **Rust watershed implementation:** `mgr_rust/src/lib.rs` — function `seeded_watershed_bucket` (see "Fix applied" section for line ranges).
- **Python caller:** `mofography/rust_fun.py` — `do_periodic_segmentation_rust` calls `mgr_rust.seeded_watershed_bucket` in several places (e.g. ~756, ~827, ~909).
- **Diagnostic script:** `docs/watershed_bucket_bug/confirm_background_voxels.py` — reproduces the bug with Case4 CIF and reports counts of void voxels with label 0.
- **Case4 CIF:** `docs/watershed_bucket_bug/Case4-combined-revised.cif`.
- **Analysis doc:** `docs/watershed_bucket_bug/watershed_bug_analysis.md` (initial analysis).
- **Watershed algorithm notes:** `docs/skimage_classic_heap_watershed_vs_rust_bucket_watershed.md` (if present in repo).

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
3. **Build/test:** `maturin build --release` in `mgr_rust`, then `pip install --force-reinstall target/wheels/...`. Running the diagnostic still showed 79,790 unlabeled and 518 in void — **no improvement**.
4. **Diagnostic logging in Rust (temporary, now removed):** Added in `lib.rs`:
   (Since removed in the fix.) It showed we never pushed to the current bucket but many levels had unlabeled voxels.
5. **Log results:** No same-bucket pushes; many levels still had unlabeled voxels, which led to the "late push to higher bucket" hypothesis.

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

## Resolved (no longer unsure)

- **Reversing** to low→high would break semantics (seeds are in high buckets; we would never expand into lower buckets). The correct fix is re-process late buckets.
- Bucket convention: `to_bucket(v)` maps image value to `0..BUCKETS-1`; we iterate `(0..BUCKETS).rev()` so high bucket index = higher image value = pore center, lower = nearer walls.
---

## Algorithm summary (for reference)

- **Constants:** `BUCKETS = 65536`. Image values scaled to bucket index via min–max.
- **Seeds:** Placed at maxima; each seed gets a label and is pushed to `buckets[bucket(seed)]`.
- **Main loop:** `for level in (0..BUCKETS).rev()` — for each level, inner `loop`: drain `buckets[level]` into `frontier`; if empty, break; else in parallel for each idx in frontier, for each 26-neighbor n in mask with label 0, `compare_exchange(0, label)`; on success push `n` to `buckets[bucket(image[n])]`. When `bucket(n) > level`, record that bucket in `late_buckets`.
- **Re-process phase:** After the main loop, repeatedly drain `late_buckets` and process each listed bucket (same propagation as above) until no bucket has pending work.
- **Output:** Copy labels to output array.

---

## How to recompile and test

1. **From repo root:**  
   `cd mgr_rust`  
   (or `cd /path/to/mofography/mgr_rust`)

2. **Build release wheel:**  
   `maturin build --release`  
   Wheel is written to `target/wheels/mgr_rust-0.1.0-cp38-abi3-<platform>.whl`.

3. **Install (force reinstall so the new backend is used):**  
   `pip install --force-reinstall target/wheels/mgr_rust-0.1.0-cp38-abi3-<platform>.whl`  
   (Use the exact filename printed by maturin.)

4. **Run diagnostic from repo root:**  
   `python docs/watershed_bucket_bug/confirm_background_voxels.py`  
   With the fix, the script should report "Background voxels IN void (BUG if > 0): 0".

5. **Optional — run tests:**  
   `pytest mofography/tests/test_mofography.py`

---

## Fix applied (re-process late buckets)

**Root cause (confirmed):** Voxels were being pushed into **already-processed (higher) buckets**. The main loop processes levels **high to low** (`(0..BUCKETS).rev()`). When at a **low** level (e.g. neck), we can label a neighbor in a **higher** bucket (e.g. cavity) and push it to `buckets[b]` with `b > level`. Level `b` was already drained and never revisited, so that voxel never propagated; masked voxels only reachable "from below" stayed unlabeled.

**Physical/geometric conditions (not rare):** The bug is *not* the "same-bucket" case (flat surfaces, not enough propagation within one bucket). That was an earlier hypothesis; an inner loop per level already ensured same-bucket voxels are drained until empty, and the bug persisted. The real trigger is **cavity voxels that are first reached only from a neck (lower distance value)**. In that situation we process the cavity's bucket (high level) before any voxel in that cavity is in the frontier; later, when we process the neck's bucket (low level), we label cavity voxels adjacent to the neck and push them to the cavity bucket—which we have already passed—so they never propagate. This geometry (narrow windows/necks opening into larger cages) is **common** in porous materials: many MOFs, zeolites, and slit-like or channel structures with constrictions exhibit it. So the physical aspects leading to the bug are **not rare**; they are typical of cage-and-window and neck-and-cavity pore space. Case4 (and similar test cases) hit the bug because of that geometry.

**Fix:** Keep high-to-low order. Record any bucket index `b` when we push to `buckets[b]` and `b > level`. After the main loop, **re-process** those buckets: repeatedly drain and run the same propagation for each such bucket, and add any bucket that receives a push during this phase back to the "to process" set, until no bucket has pending work.

**Code references (main repo, `mgr_rust/src/lib.rs`):**

- **Late-bucket set:** Lines 1013-1015 — `late_buckets: Arc<Mutex<HashSet<usize>>>` created before the main loop.
- **Recording late pushes (main loop):** Lines 1054-1058 — inside the propagation closure, after `let b = bucket_idx_vec[nidx]`, `if b > level { late_buckets.lock().insert(b); }` then `buckets[b].lock().push(nidx)`.
- **Re-process loop:** Lines 1065-1115 — after the main `for level in (0..BUCKETS).rev()`, a `loop` drains `late_buckets` into `to_process`; for each `b` in `to_process`, an inner `loop` drains `buckets[b]` and runs the same `frontier.par_iter().for_each(...)` propagation; on each push to `b'` we do `late_buckets.lock().insert(b_prime)` so that bucket is processed again.
- **Performance:** The per-level full scan over all voxels (previously after each level) and the hot-path atomic/`eprintln!` diagnostics were removed so the 65k-level loop only does drain + propagate. A single post-loop unprocessed count is kept under `#[cfg(debug_assertions)]` (lines 1117-1129).

**Verification:** Run `docs/watershed_bucket_bug/confirm_background_voxels.py`; "Background voxels IN void (BUG if > 0)" should be 0. `pytest mofography/tests/test_mofography.py` should pass.

---

## Files to inspect

- **`mgr_rust/src/lib.rs`:** function `seeded_watershed_bucket`. Search for "Main flooding loop", "late_buckets", "Re-process late buckets". The main loop, late-bucket recording (`if b > level`), and the re-process loop after the main loop implement the fix.

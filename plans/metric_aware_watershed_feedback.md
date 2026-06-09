# Review: `metric_aware_watershed_plan.org`

## Verdict

Plan direction strong. Root-cause framing mostly right: non-orthogonal cell makes 26-NN graph physically anisotropic, current flood treats neighbors too uniformly.  
But one key algorithm detail needs correction before merge.

## High-severity findings (fix before implementation)

1. **Priority update rule in proposal is not path-consistent yet.**  
   Proposed inner-loop `adjusted = raw / w` then push can make same voxel receive different priorities depending on last hop only.  
   For watershed/flood semantics, metric penalty should be applied **inside bottleneck propagation**, not as pure local replacement.

   Recommended update:
   - Compute `cand = img_flat[nidx] * inv_w[offset]` (or `/ w`).
   - Then propagate bottleneck: `next = min(parent_prio, cand)`.
   - Push `next`.

   This keeps monotone flood property and avoids path-order artifacts.

2. **Doc says “orthogonal cells unchanged” for all 26 neighbors; mathematically false if normalized by mean face length only.**  
   In orthogonal grid, corner/edge lengths are not equal to face length.  
   If `w = L / mean(face)` then orthogonal diagonals get penalized (`sqrt(2), sqrt(3)`), which changes behavior.

   Two clean options:
   - **Option A (recommended):** define baseline `L0(offset)` from orthogonal geometry and set `w = L / L0(offset)`. Then orthogonal case stays exactly unchanged.
   - **Option B:** explicitly accept behavior change and update claim/tests accordingly.

3. **“Expected overhead <1%” likely optimistic without microbenchmark.**  
   Heap flood neighbor loop hot path. One extra op can still be noticeable at scale.  
   Use reciprocal multiply (`inv_w`) and benchmark before claiming.

## Medium-severity findings

4. **Need deterministic tie behavior.**  
   Weighted priorities increase ties/near-ties. Ensure stable tie-break (existing index order or secondary key) preserved.

5. **Need explicit guardrails for extreme skew.**  
   Clamp weights (`w_min`, `w_max`) or alpha schedule to prevent over-penalizing one direction in pathological cells.

6. **Validation target should include segmentation topology metrics, not only PLD.**  
   Add:
   - component count
   - percolation class
   - boundary/saddle dgrid histogram
   - IoU vs reference mask (if available)

## Will this likely fix issue?

**Likely yes for shortcut leakage class**, if implemented with path-consistent bottleneck update and orthogonal-invariant normalization.  
Confidence: **medium-high** for VIGNOQ/TOVKOG type failures, **medium** for all PLD mismatches (some may still involve maxima/merge heuristics downstream).

## Faster path into Rust backend

### Fastest safe rollout (2-stage)

1. **Stage 1 (quick, low-risk): anisotropy gate + feature flag**
   - Compute anisotropy ratio from 26 lengths once (`max/min` or corner spread).
   - If ratio < threshold (e.g. `1.08`), keep current fast path (all weights = 1).
   - Else use weighted path.
   - Behind flag `metric_aware=true` default-on only for skewed cells.

2. **Stage 2 (performance polish)**
   - Precompute fixed arrays:
     - `NEIGHBOR_OFFSETS: [(i8,i8,i8); 26]`
     - `INV_W: [f32; 26]`
   - Replace division with multiply.
   - Keep arrays contiguous, avoid 3D index lookup in hot loop.
   - Reuse computed metric per call; no alloc inside loop.

### Suggested Rust core formula

```rust
let cand = img_flat[nidx] * inv_w[k];   // metric-aware local capacity
let next = cand.min(parent_prio);       // bottleneck propagation
heap.push((NotNan::new(next).unwrap(), nidx));
```

## Suggested edits to plan

1. In formulas section, replace `Priority = dgrid / w[offset]` with  
   `Candidate = dgrid / w[offset]; next_prio = min(parent_prio, Candidate)`.

2. Replace “orthogonal cells unchanged” statement with explicit invariant definition using `L0(offset)`.

3. Add benchmark acceptance criteria:
   - runtime overhead target (e.g. `<5%` on representative grids)
   - no regression on orthogonal corpus
   - deterministic output bitwise-equal when anisotropy gate not triggered

4. Add one ablation:
   - compare `alpha = 1.0` vs `0.5` in `cand = dgrid / w^alpha`.

## Bottom line

Best direction among listed options.  
Not ready exactly as written.  
Small math/algorithm corrections make it robust and still fast for Rust backend.

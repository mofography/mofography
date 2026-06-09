# Watershed Over-Segmentation: Extra Maxima in Channels

**Date:** 2026-06-07  
**Status:** Draft — needs senior review  
**Affected structures:** REHHOD_clean.cif (Cd-MOF), likely any MOF with large-VDW atoms in channel walls  
**Impact:** PLD can be underestimated by up to **40%** (4.30 Å vs 7.28 Å for REHHOD)

---

## 1. Problem Statement

When using Mendeleev VDW radii (Cd = 2.18 Å), mofography's watershed segmentation creates **extra pocket nodes in narrow channels** where there should be **windows (saddle points)**. This causes the winding percolation algorithm to compute an artificially low PLD.

### Observed Behavior (REHHOD_clean.cif, sp=0.25 Å)

| Method | PLD (Å) | Pore Graph Nodes |
|---|---|---|
| zeo++ (Mendeleev) | **7.28** | N/A (Voronoi network) |
| mofography + CRC radii | **7.20** | 4 (all large, r≈4.6 Å) |
| mofography + Mendeleev radii | **4.30** | 6 (2 large at r≈4.7, 4 small at r≈2.9) |

The largest inter-pocket window is 3.60 Å in **both** cases (CRC and Mendeleev). But with Mendeleev, the percolating path is forced through narrow windows (~2.15 Å) between small channel-fragment nodes.

---

## 2. Root Cause Analysis

### 2.1 Extra Maxima in Channels

With large VDW radii (Cd = 2.18 Å in Mendeleev), the distance grid develops **local maxima inside narrow channels** between cavities. The watershed's `peak_local_max` step detects these as separate maxima, creating extra pocket regions.

These channel maxima should be **saddle points** (windows), not pocket centers. The consequence:

- **CRC (Cd = 1.58 Å):** Channels are wide → distance grid has 4 true cavity maxima → 4 pockets → multiple independent wide paths → non-contractible loop at wr = 3.60 → **PLD = 7.20**
- **Mendeleev (Cd = 2.18 Å):** Channels are narrow → distance grid has 2 cavity maxima + 4 channel maxima → 6 pockets (2 large + 4 small) → only one wide path (1→4 at wr=3.60) → no non-contractible loop at high threshold → **PLD = 4.30**

### 2.2 Why Winding Percolation Fails with Over-Segmented Graphs

Winding percolation requires a **non-contractible loop** (closed path with non-zero lattice offset). The BFS traverses edges in both directions and detects when the same node is reached at two different accumulated offsets.

**CRC graph (4 nodes):**
```
Path 1: 1 →2(1,0,0) →4(0,1,0)  = node 4 at offset (1,1,0)
Path 2: 1 →3(0,-1,0)→4(0,0,-1) = node 4 at offset (0,-1,-1)
→ winding = (1,2,1) → PERCOLATION at bottleneck 3.60 Å → PLD = 7.20
```

**Mendeleev graph (6 nodes):**
```
At wr ≥ 3.60: only one edge survives (1→4 at offset (1,0,0))
BFS: 1 → 4 at (1,0,0) → reverse traversal back to 1 at (0,0,0) → winding = (0,0,0)
→ NO PERCOLATION. Must lower to wr = 2.15 for loop through small nodes.
```

With only 2 large nodes, any out-and-back traversal has **zero winding**. The 4 small nodes don't provide alternative wide paths — they only connect through narrow windows.

### 2.3 Why `get_pore_dia_mgr_winding2` Doesn't Catch This

The spill guardrail in `winding2` checks `_has_suspicious_self_loop_asymmetry`: it flags nodes whose self-loop radii spread exceeds 30% of pocket radius. For REHHOD, self-loop spread/pocket_r = 0.032 (well below 0.30). The guardrail didn't trigger.

The guardrail is designed for a different pathology (watershed spill across thin boundaries by 1-2 voxels). Our issue is **extra maxima in channels**, which is a distinct segmentation artifact.

### 2.4 Parameter Variations Don't Help

All tested parameter combinations produce the same PLD = 4.30:

| Parameter | Values Tested | PLD |
|---|---|---|
| `watershed_backend` | `rust_det`, `rust_bucket` | 4.30 |
| `mask_thickness` | 0.05, 1.5 | 4.30 |
| `h` (watershed depth) | 0.2, 0.5, 1.0, 2.0, 4.0 | 4.30 |
| `min_distance` | 1.5, 2.0, 3.0, 4.0, 5.0 | 4.30 |

The maxima are present regardless of flooding parameters. Higher `h` values reduce to 2 nodes (both large), but 2 nodes can never form a non-trivial loop.

---

## 3. Proposed Diagnostics (Detectable from Pore Graph)

Three checks that can flag over-segmentation:

### Check 1: Bimodal Pocket Size Distribution
```python
pocket_radii = [d['pocket_radius'] for _, d in G.nodes(data=True)]
ratio = max(pocket_radii) / min(pocket_radii)
if ratio > 1.5:
    flag("Possible over-segmentation: pocket size ratio = {:.2f}".format(ratio))
```
- CRC: ratio = 1.00 (healthy)
- Mendeleev: ratio = 1.61 (flagged)

### Check 2: PLD vs Max-Window Gap
```python
max_window = max(d['window_radius_min'] for u, v, d in G.edges(data=True) if u != v)
pld_ratio = pld / (2 * max_window)
if pld_ratio < 0.85:
    flag("PLD/max_window gap: ratio = {:.2f} (expected ~1.0)".format(pld_ratio))
```
- CRC: ratio = 1.00 (healthy)
- Mendeleev: ratio = 0.60 (flagged)

### Check 3: Path Redundancy at Max-Window Threshold
```python
# Filter graph to edges with wr >= max_window - 0.1
# Check: does any node have >= 2 independent wide paths?
# If not, the max window cannot participate in a percolating loop
```
- CRC at wr ≥ 3.50: 4 edges survive, 2 independent paths to node 4 → healthy
- Mendeleev at wr ≥ 3.50: 1 edge survives → flagged

---

## 4. Proposed Fixes (Three Options, Ranked by Effort)

### Option A: Simple-Path Percolation Fallback (Recommended — Quick Win)

When the pore graph fails the winding percolation (PLD comes out suspiciously low), fall back to a **simple-path percolation** check: does any inter-node edge with non-zero periodic offset exist at this threshold?

```python
if n_unique_cavity_nodes <= 3 or pld_ratio < 0.85:
    # Simple-path fallback: percolation = existence of any inter-node edge
    # with non-zero offset at this threshold
    max_wr_nonzero = max(
        d['window_radius_min']
        for u, v, d in G.edges(data=True)
        if u != v and any(o != 0 for o in d['relative_offset'])
    )
    pld_fallback = 2 * max_wr_nonzero
    if pld_fallback > pld:
        pld = pld_fallback  # Use the better estimate
```

**Result:** PLD = 7.20 Å (matches zeo++ within 1.1%)  
**Risk:** May overestimate PLD for structures where winding percolation's loop requirement is physically meaningful (i.e., where a single wide path doesn't imply true percolation).

### Option B: Maxima Suppression in Channels (Medium Effort)

Modify the segmentation pipeline to detect and suppress channel maxima before watershed:

1. After `peak_local_max`, identify pairs of maxima that are "channel maxima" (local max in a narrow passage between two larger maxima)
2. Suppress channel maxima by:
   - Increasing `min_distance` adaptively based on local grid curvature
   - Or: merge maxima whose distance grid values are below a fraction of the global maximum
   - Or: use the Hessian eigenvalues at each maximum to distinguish cavities (one negative eigenvalue) from channels (two negative eigenvalues)

**Risk:** May suppress genuine small cavities. Needs testing on a diverse structure set.

### Option C: Voronoi-Based PLD (High Effort)

Implement a Voronoi network-based PLD calculation as an alternative to watershed-based approach, similar to zeo++:

1. Build Voronoi network from atomic positions + VDW radii
2. Compute bottleneck radius along each Voronoi edge
3. Find widest path through the network using Dijkstra-like algorithm
4. PLD = 2 × minimum bottleneck radius along the widest percolating path

**Risk:** Significant implementation effort. Would essentially be reimplementing zeo++ in Python.

---

## 5. Current Mitigation

The default VDW radii in mofography have been switched from Mendeleev to CRC/Bondi (1964) in `mofography/vdw_radii.py`. This eliminates the over-segmentation for REHHOD and similar Cd-MOFs. However, the issue may still affect structures with intrinsically narrow channels regardless of radii choice.

---

## 6. Recommended Next Steps

1. **Review this plan** with senior team
2. **Decide on fix strategy:** Option A (fallback) is low-risk and immediately useful
3. **Test on diverse structures:** Run the diagnostic checks on the full CCR dataset (2000 CIFs) to assess how common this pathology is
4. **Implement chosen fix** and validate against zeo++ on a representative subset
5. **Add diagnostic warnings** to `get_pore_dia_mgr_winding` output when over-segmentation is detected

---

## Appendix: Key Files

- `mofography/vdw_radii.py` — CRC/Bondi radii (new default) and Mendeleev (legacy)
- `mofography/winding_percolation.py:888` — `get_pore_dia_mgr_winding`
- `mofography/winding_percolation.py` — `get_pore_dia_mgr_winding2` (spill guardrail)
- `mofography/winding_percolation.py` — `_has_suspicious_self_loop_asymmetry` (current spill check)
- `mofography/winding_percolation.py` — `percolation_winding_multidigraph` (BFS core)
- `mofography/winding_percolation.py` — `_iter_neighbor_edges` (bidirectional traversal)
- `mofography/rust_fun.py:33` — `dgrid_from_atoms_rust`
- `mofography/rust_fun.py:258` — `do_periodic_segmentation_rust`

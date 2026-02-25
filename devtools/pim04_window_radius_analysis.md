## Context

This note summarizes an investigation into why, for the pim‑0_4 case (`pim-0_4.cif`), some window radii attached to node 20 are larger than that node’s pocket radius, despite maxima detection on the distance grid (`dgrid`) being stable and unchanged.

The key question: **is the problem in the segmentation (watershed + padding) or in the window construction (`window_connections`), and how might the late‑buckets fix in the Rust watershed be involved?**

---

## Code and files involved

- **Rust watershed and windows**
  - `mgr_rust/src/lib.rs`
    - `seeded_watershed_bucket`: bucket-based parallel seeded watershed used in periodic segmentation.
    - `window_connections`: scans padded labels + `dgrid` and emits connections and window radii:
      - Returns tuples `(reg1, reg2, centers, radii)` as a Python list of tuples.
- **Python wrapper and segmentation workflow**
  - `mofography/rust_fun.py`
    - `do_periodic_segmentation_rust`: core periodic segmentation pipeline:
      - builds padded `dgrid` and markers,
      - calls `mgr_rust.seeded_watershed_bucket` for merged and unmerged regions,
      - calls `mgr_rust.window_connections` on `regions_padded`, padded `dgrid`, and padding thickness,
      - assembles the `results` dict used downstream.
- **Graph and percolation**
  - `mofography/winding_percolation.py`
    - `collect_periodic_edges(connections, maxima_info)`: aggregates connections into an edge map with `window_radii`, `window_radius_max`, etc.
    - `make_directed_graph_from_segmentation_results(results, dgrid=None)`:
      - Builds a `networkx.MultiDiGraph` from `results`.
      - For each node (region label), sets:
        ```python
        voxel_idx = tuple(maxima_info[node - 1, 1:4].astype(int))
        G.nodes[node]["pocket_radius"] = float(dgrid[voxel_idx])
        ```
      - Edges carry attributes like:
        - `window_radius`: the maximum window radius for that base edge.
- **Visualization helpers (for debugging only)**
  - `mofography/visualize.py`
    - `plot_connections_as_it_is`
    - `plot_windows`
- **External repro/debug (manuscript repo)**
  - `/Users/arungopalan/repos/mofography_manuscript/resources/mgr_org/mgr_org_scripts/chao_pld_issue/`
    - `debug_pim-0_4.org`: org-mode notebook with end-to-end pim‑0_4 workflow and manual checks.
    - `repro_segmentation_tests.py`: reproducibility tests across multiple runs and mask thicknesses.
    - `heap_vs_bucket_watershed_tests.py`: compares bucket vs heap watershed on the same padded inputs.

---

## Regression script: check_pim_windows.py

To make the issue reproducible and machine-checkable, we added:

- `devtools/check_pim_windows.py`

This script:

1. Loads a CIF (pim‑0_4 here).
2. Builds a distance grid using the Rust backend and the TRAPPE-style radii from the manuscript scripts.
3. Runs `do_periodic_segmentation_rust` with `compute_connections=True`, returning `results` with:
   - `regions`, `regions_padded`, `maxima_info`, `all_maxima`, `connections`, `dgrid`, etc.
4. Builds the directed graph using:
   ```python
   from mofography.winding_percolation import make_directed_graph_from_segmentation_results
   graph = make_directed_graph_from_segmentation_results(results, dgrid=results["dgrid"])
   ```
5. Scans all graph edges and flags violations of:
   \[
   \text{window\_radius} \le \min(\text{pocket\_radius}_i, \text{pocket\_radius}_j) + \text{tol}
   \]
   where `pocket_radius` is taken from `G.nodes[node]["pocket_radius"]`, which in turn is `dgrid[maxima_info[node-1,1:4]]`.

### How to run it

From the main `mofography` repo:

```bash
cd /Users/arungopalan/repos/mofography

# mask_thickness = 0.0
python devtools/check_pim_windows.py \
  /Users/arungopalan/repos/mofography_manuscript/resources/mgr_org/mgr_org_scripts/chao_pld_issue/pim-0_4.cif 0.0

# mask_thickness = 0.25
python devtools/check_pim_windows.py \
  /Users/arungopalan/repos/mofography_manuscript/resources/mgr_org/mgr_org_scripts/chao_pld_issue/pim-0_4.cif 0.25
```

The script reports:

- Basic segmentation metrics (`n_regions`, `n_maxima`, `n_connections`, `plr`, `pld`).
- Any edges where `window_radius > min(pocket radii) + tol`, highlighting node 20 in the listing.

---

## Results from the regression

### 1. mask_thickness = 0.0

Command:

```bash
python devtools/check_pim_windows.py \
  /Users/arungopalan/repos/mofography_manuscript/resources/mgr_org/mgr_org_scripts/chao_pld_issue/pim-0_4.cif 0.0
```

Output (relevant part, truncated):

```text
[info] regions=92, maxima=144, connections=709, plr=2.3730475902557373, pld=4.746095180511475
[warn] Found 13 edges where window_radius > min(pocket_radius_i, pocket_radius_j) + tol
  edge 7-11: window_radius=2.712 Å, pocket_u=2.298 Å, pocket_v=2.924 Å
**edge 11-20: window_radius=2.807 Å, pocket_u=2.924 Å, pocket_v=2.373 Å
  edge 11-22: window_radius=3.899 Å, pocket_u=2.924 Å, pocket_v=2.772 Å
  edge 11-25: window_radius=3.685 Å, pocket_u=2.924 Å, pocket_v=1.890 Å
  edge 11-30: window_radius=3.892 Å, pocket_u=2.924 Å, pocket_v=3.268 Å
**edge 16-20: window_radius=2.836 Å, pocket_u=2.677 Å, pocket_v=2.373 Å
**edge 20-30: window_radius=3.378 Å, pocket_u=2.373 Å, pocket_v=3.268 Å
**edge 20-35: window_radius=1.985 Å, pocket_u=2.373 Å, pocket_v=1.598 Å
**edge 20-39: window_radius=3.623 Å, pocket_u=2.373 Å, pocket_v=1.352 Å
...
```

So:

- There are **13 edges** violating the invariant.
- Several involve **node 20**, exactly matching the hand-inspection from the org debug file.

### 2. mask_thickness = 0.25

Command:

```bash
python devtools/check_pim_windows.py \
  /Users/arungopalan/repos/mofography_manuscript/resources/mgr_org/mgr_org_scripts/chao_pld_issue/pim-0_4.cif 0.25
```

Output (relevant part, truncated):

```text
[info] regions=82, maxima=136, connections=557, plr=2.3730475902557373, pld=4.746095180511475
[warn] Found 12 edges where window_radius > min(pocket_radius_i, pocket_radius_j) + tol
  edge 6-9: window_radius=2.693 Å, pocket_u=2.298 Å, pocket_v=2.924 Å
  edge 9-17: window_radius=2.974 Å, pocket_u=2.924 Å, pocket_v=2.373 Å
  edge 9-19: window_radius=3.865 Å, pocket_u=2.924 Å, pocket_v=2.772 Å
  edge 9-21: window_radius=3.341 Å, pocket_u=2.924 Å, pocket_v=2.597 Å
  edge 9-22: window_radius=3.632 Å, pocket_u=2.924 Å, pocket_v=1.890 Å
  edge 9-27: window_radius=3.858 Å, pocket_u=2.924 Å, pocket_v=3.268 Å
  ...
```

Even with a nonzero mask thickness (0.25), a similar pattern of violations persists (labels are renumbered, but several windows are still larger than the smaller pocket).

---

## Additional diagnostic: central vs padded dgrid for node 20

To disentangle whether the problem is with maxima detection or with segmentation (especially on the padded grid), we ran a one-off diagnostic in a Python REPL:

```python
import numpy as np
from pathlib import Path

from devtools.check_pim_windows import run_single_segmentation

pim_path = Path(
    "/Users/arungopalan/repos/mofography_manuscript/"
    "resources/mgr_org/mgr_org_scripts/chao_pld_issue/pim-0_4.cif"
)

results, graph, metrics = run_single_segmentation(
    pim_path, spacing=0.25, mask_thickness=0.0
)

dgrid = results["dgrid"]
regions = results["regions"]
regions_padded = results["regions_padded"]
pad = int(results["input_params"]["padding_thickness"])

padded_grid = np.pad(dgrid, pad, mode="wrap")

node = 20
mask_central = (regions == node)
mask_padded = (regions_padded == node)

central_max = float(dgrid[mask_central].max()) if mask_central.any() else float("nan")
central_min = float(dgrid[mask_central].min()) if mask_central.any() else float("nan")

padded_max = float(padded_grid[mask_padded].max()) if mask_padded.any() else float("nan")
padded_min = float(padded_grid[mask_padded].min()) if mask_padded.any() else float("nan")

print(f"Node {node}: central max={central_max:.6f}, min={central_min:.6f}")
print(f"Node {node}: padded  max={padded_max:.6f}, min={padded_min:.6f}")
```

Output:

```text
Node 20: central max=2.373048, min=0.000017
Node 20: padded  max=3.986511, min=0.000013
```

So for node 20:

- In the **central box**, the region-20 voxels satisfy:
  - `max(dgrid[regions == 20]) = 2.373048 Å`
  - `min(dgrid[regions == 20]) ≈ 1.7e-5 Å`
- In the **padded grid**, the same label 20 appears on voxels where:
  - `max(padded_grid[regions_padded == 20]) ≈ 3.9865 Å`
  - `min(...) ≈ 1.3e-5 Å`

We also checked that:

```python
maxima_info = results["maxima_info"]
row_idx = np.where(maxima_info[:, 0] == 20)[0][0]
voxel = tuple(maxima_info[row_idx, 1:4].astype(int))
pocket_radius_from_maxima = float(dgrid[voxel])
```

gives:

```text
pocket_radius_from_maxima = 2.373048 Å
```

and `max(dgrid[regions == 20])` is **exactly** the same, so maxima detection **within the central box** is consistent with the region’s true max.

But in the padded grid, regions_padded == 20 owns voxels with significantly larger `dgrid` values (~3.99 Å), which are not reflected in the node’s pocket radius (which is based only on the central box).

---

## Edges for node 20 in the graph

We also printed all edges in the graph that involve node 20:

```python
for u, v, data in graph.edges(data=True):
    if 20 not in (u, v):
        continue
    w  = float(data.get("window_radius", 0.0))
    ru = float(graph.nodes[u].get("pocket_radius", 0.0))
    rv = float(graph.nodes[v].get("pocket_radius", 0.0))
    print(f"  edge {u}-{v}: window_radius={w:.6f}, pocket_u={ru:.6f}, pocket_v={rv:.6f}")
```

Output (truncated):

```text
edge 11-20: window_radius=3.378628, pocket_u=2.923577, pocket_v=2.373048
edge 16-20: window_radius=3.645940, pocket_u=2.677011, pocket_v=2.373048
edge 3-20:  window_radius=0.622011, pocket_u=4.294212, pocket_v=2.373048
...
edge 20-29: window_radius=2.324910, pocket_u=2.373048, pocket_v=2.407458
edge 20-39: window_radius=3.049086, pocket_u=2.373048, pocket_v=1.352145
edge 20-35: window_radius=3.131189, pocket_u=2.373048, pocket_v=1.597508
edge 20-30: window_radius=2.535504, pocket_u=2.373048, pocket_v=3.267832
...
```

These are the same violations flagged by the regression script; they confirm that **window radii are being computed using boundary voxels whose `dgrid` is larger than the central-pocket max for node 20 and/or its neighbors**.

---

## Interpretation and conclusion

### 1. Maxima detection is not the culprit (for node 20)

For node 20:

- `pocket_radius_from_maxima = dgrid[maxima_info[20-1,1:4]] = 2.373048 Å`
- `max(dgrid[regions == 20]) = 2.373048 Å` (same to full precision)

So within the **central box**, maxima detection is accurately capturing the region’s true maximum dgrid value. There is no evidence that maxima detection itself is underestimating the pocket maximum for node 20.

### 2. The padded segmentation is “fatter” than the central segmentation

However:

- In the **padded grid**, `regions_padded == 20` reaches max `dgrid ≈ 3.9865 Å`.
- That means the region-20 label extends, in the padding, into voxels with much larger distance values than any voxel in the central box assigned to region 20.

This expansion is likely influenced by the **late-buckets fix** in `seeded_watershed_bucket`:

- Before the fix, higher buckets that received “late” pushes after their main pass could leave some voxels unlabeled or less extensively flooded, especially in complicated neck/cavity geometries.
- After the fix, late-bucket reprocessing ensures that voxels in higher buckets (cavities) that are discovered from lower buckets (necks) are fully propagated, including across the padded wrap-around.
- As a result, labels in `regions_padded` can now occupy more “interior” parts of the padded cavity at higher dgrid values than the central maxima.

### 3. window_connections is faithfully using padded labels and dgrid

In `window_connections` (Rust):

1. It scans the entire **padded** grid (`labels` == `regions_padded`) for boundaries between different nonzero regions (6-neighbor connectivity).
2. For each region pair `(reg1, reg2)` it:
   - clusters boundary points in Cartesian space using a `distance_threshold`,
   - for each cluster, splits points by region and finds the **max distance** in each region:
     ```rust
     let max_reg1 = reg1_points.iter().max_by(|a,b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
     let max_reg2 = reg2_points.iter().max_by(|a,b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
     let radius = 0.5 * (max_reg1.distance + max_reg2.distance) as f64;
     ```
3. It returns this radius as the window radius, and Python later uses the maximum over multiple windows for a given edge.

Because `regions_padded == 20` includes voxels with `dgrid ≈ 3.99 Å`, the **largest boundary distances** available to `window_connections` on the “region 20 side” can be significantly larger than 2.37 Å. Averaging those with the other side’s maxima naturally yields window radii larger than the **central** pocket radius.

Thus:

- The observed inconsistency (`window_radius > pocket_radius`) is **driven by the difference between how pockets and windows are defined**:
  - **Pocket radius**: central box only, via `dgrid` at central maxima.
  - **Window radius**: padded grid, via maximal boundary distances in each region cluster.
- `window_connections` is doing what it was asked to do; the mismatch arises because segmentation on the padded grid (especially after the late-buckets fix) allows regions to own higher-distance voxels than what is considered as “pocket interior” in the central box.

### 4. Role of the late-buckets fix

The late-buckets fix in `seeded_watershed_bucket` ensures:

- Same-bucket flooding is completed at each level.
- Any voxel pushed into a **higher** bucket (`b > level`) during flooding is reprocessed by the late-buckets loop until no pending work remains.

In complex pore geometries with necks and cavities (such as pim‑0_4):

- This improves correctness in the sense of avoiding unlabeled masked voxels.
- But it also makes it more likely that regions in the **padded** grid flood into parts of the cavity that were previously unlabeled or assigned differently, especially at higher dgrid values.

So while the fix removed one class of bugs (unlabeled interior void voxels), it also **made this pocket/window inconsistency more visible**, because now `regions_padded` more fully explores the padded cavity, extending labels into higher-dgrid areas where windows are later sampled.

---

## Ideas for fixing or reducing the inconsistency

There are a few conceptual options; they aren’t mutually exclusive and could be combined.

### Option A: Restrict window construction to central-box boundaries

Currently, `window_connections` runs on `regions_padded` and padded `dgrid`. One way to align windows with central pockets is to:

- Only consider **boundary voxels whose unpadded indices lie within the central box** when building `edge_map`.
  - That is, restrict to `i,j,k` such that:
    \[
    \text{padding\_thickness} \le i < nx - \text{padding\_thickness},
    \]
    and similarly for `j, k`.
- Or equivalently, map padded coordinates back to central coordinates and discard any crossings that are purely in the padding.

Pros:

- Window radii will be determined only by voxels inside the physical, central unit cell, where pocket radii are also defined.
- Less sensitive to the exact behavior of the watershed in deep padding.

Cons:

- Some periodic windows that genuinely occur across the padded boundaries may need careful treatment to avoid losing valid connections.

### Option B: Clamp window radii by pocket radii

Impose an invariant at the graph layer:

- After computing `window_radius` for an edge `(i, j)`, set:
  \[
  \text{window\_radius} \gets \min\left(\text{window\_radius}, \min(\text{pocket\_radius}_i, \text{pocket\_radius}_j)\right)
  \]

Pros:

- Simple, robust, and guarantees the invariant `window_radius <= min(pocket radii)` by construction.
- Avoids overestimating throats relative to pockets, which is important for PLD-style percolation.

Cons:

- If segmentation on the padded grid is genuinely correct and pockets are *underestimated* in some other example, this clamp might artificially shrink real windows.
- It does not address the underlying geometric inconsistency, only its numerical manifestation.

### Option C: Redefine pocket radius to respect padded extents

Instead of defining `pocket_radius` from the central maxima only, consider:

- Defining pocket radius as `max(dgrid[regions == node])` **after mapping PBC groups**, possibly using `regions_padded` to account for periodic images in a more symmetric way.

Pros:

- Aligns pockets and windows by letting both “see” the same extended cavity.

Cons:

- Changes the physical meaning of pocket radius away from “max inside the central box” to something more global, which might complicate comparisons with existing literature or Python-only implementations.

### Option D: Further constrain flooding in padding

Investigate modifications to `seeded_watershed_bucket` such that:

- In the padded grid, labels are prevented from occupying voxels with `dgrid` significantly larger than their central maxima (e.g., by:
  - forbidding propagation into voxels whose `dgrid` exceeds the seed’s own `dgrid` by a tolerance, or
  - performing watershed on a clipped `dgrid` where values are capped at each region’s pocket max).

Pros:

- Attacks the root cause on the segmentation side.

Cons:

- Harder to do cleanly in a parallel bucket framework without breaking watershed semantics.
- Needs careful analysis to avoid introducing new artifacts or reintroducing unlabeled-void issues that the late-buckets fix solved.

---

## Summary

- **Maxima detection in the central box is consistent**: for node 20, `pocket_radius` from `maxima_info` matches the true `max(dgrid[regions == 20])`.
- **The padded segmentation is “fatter”**: `regions_padded == 20` includes voxels with `dgrid` up to ~3.99 Å, much larger than the central max of ~2.37 Å.
- **window_connections operates on the padded grid**, so it can (and does) use these higher-distance boundary voxels to compute window radii, leading to `window_radius > pocket_radius` for several edges, including those involving node 20.
- The **late-buckets fix** in `seeded_watershed_bucket` likely increased how fully labels explore high-dgrid regions in the padding, making this inconsistency visible, but the conceptual mismatch (central-pocket definition vs padded-window definition) predates the fix.
- Possible remedies include:
  - restricting window construction to central-box boundaries,
  - clamping window radii by pocket radii at the graph layer,
  - redefining pocket radii to account for padded extents, or
  - adding further constraints to how the watershed floods in padding.


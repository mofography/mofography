# Mofography fork sync changelog

**Prepared:** 2026-06-11 (updated after push to `master`)  
**Author:** Arun Gopalan (`arung-northwestern`)  
**Upstream fork:** [Shi-Research-Group/mofography](https://github.com/Shi-Research-Group/mofography) (`master`)  
**Source repo:** [arung-northwestern/mofography](https://github.com/arung-northwestern/mofography) — branch **`master`** (same tip as `clean-master`)

This document summarizes what **`arung-northwestern/mofography` @ `2d63da9`** contains that **Shi-Research-Group/mofography** (`master`, last pushed **2026-03-25**) does not. It is intended for manual fork updates on GitHub.

> **Ready for GitHub UI:** As of **2026-06-11**, `master` and `clean-master` on `arung-northwestern/mofography` both point to **`2d63da9`** and are pushed. You can open a cross-fork PR **now** (see [Fork update workflow](#fork-update-workflow-github-ui) below).

> **Note on org name:** the Shi group repo is under **`Shi-Research-Group`**, not `shi-group-github`.

---

## At a glance

| Metric | Value |
|--------|------:|
| Common ancestor | `1444c2e` — *consolidate the overlap leakage stats* (2026-03-25) |
| Commits on Northwestern `master` since ancestor | **107** |
| Shi `master` tip | `1cce57f` (merge from `arung-northwestern:master`, 2026-03-25 era) |
| Northwestern tip (`master` **=** `clean-master`) | **`2d63da9`** (2026-06-11) |
| Ahead / behind vs Shi `master` | **107 ahead**, **40 behind** (Shi merge commits only) |

**Branch naming (updated):** Northwestern development lived on `clean-master` during the spring; **`master` and `clean-master` are now identical** on `arung-northwestern/mofography`. For Shi fork updates, use **`arung-northwestern/mofography` → `master`** (simplest — matches Shi’s default branch name).

**Shi-only history:** Shi `master` contains ~40 merge commits (`Merge branch 'arung-northwestern:master' into master`) not on Northwestern `master`. These are integration noise, not separate feature work. The substantive delta is the **107 commits** since `1444c2e`.

### Commits after the first changelog draft (`6c36b13` → `2d63da9`)

| Commit | Summary |
|--------|---------|
| `b487af5` | Drop failed metric/blended watershed prototypes; default **`rust_det`**; add this changelog. |
| `840ee9a` | Untrack agent state (`graphify-out/`, Letta, Serena, …); move root dev scripts to `scripts/dev/`; repo hygiene. |
| `c1e38db` | Update this changelog for synced `master`/`clean-master` and GitHub UI fork instructions. |
| `2d63da9` | Changelog tip/hash alignment (doc-only). |

---

## Highlights (new capabilities)

These are the largest functional additions since the March 2025 fork point.

### 1. Interactive Three.js / Figure3D visualization (`mofography/threejs_viz.py`, `threejs_plotting.py`, `viz_backend.py`)

- Plotly-like Python API for standalone Three.js HTML export (static scenes + movies).
- Rich viewer HUD: orthographic toggle, axis views, exploded cell view, mesh opacity, stick/spacefill, GIF/WebM recording, lazy sidecar loading for large movies.
- Instanced framework representation; aromatic ring rendering; periodic bond fixes for flexible movies.
- **Key commits:** `62adc898` (Figure3D API), `b6547f6` (instanced framework), `60acaccd` (movie HTML module), `f5f1c7b` / `9186ef8` (movies + percolation scenes).

### 2. Deterministic watershed segmentation (Rust)

- **`rust_det`** is now the **default** production backend in `do_periodic_segmentation_rust` (was `rust_bucket`).
- **`rust_det_parallel`** and **`rust_bucket`** remain available for benchmarks.
- Experimental **metric / blended** watershed backends (`rust_det_metric`, `rust_det_blended`, `rust_det_max`, `watershed_blended/` crate) were **prototyped and removed** — they did not improve parity on skewed cells and are not in the current tree.

### 3. Optimized distance grid (`dgrid`) for non-orthogonal cells (`mgr_rust/src/dgrid_exp.rs`)

- Tiered exact MIC (Wassenaar 7-image + 27-image paths), SoA + SIMD orthogonal fast path.
- Promoted optimized path as default `dgrid_from_atoms`.
- Mendeleev default VDW radii in Rust wrappers; CRC/Bondi radii module in Python (`mofography/vdw_radii.py`).
- **Key commits:** `305a13a`, `6773f638`, `96e5c10c`, `606c7c48`, `b10999a`.

### 4. Pocket typing pipeline (`mofography/pocket_typing.py`, `pocket_framework_env.py`)

- Ported/upstreamed from `mgr_cli_py`: PSD + annulus fingerprints + connectivity clustering.
- Benchmarked end-to-end (7/7 pass after fixes).
- **Key commits:** `dd19eaf`, `81364a4`, `328bf22`, `542edd1`.

### 5. Window fingerprinting for accessible textural properties (`mofography/winding_percolation.py`, `mofography/compute.py`)

- Shell-vector fingerprints on pore-graph windows; BFS rescue for sub-probe grid undershoot (Mechanism B).
- Default **on** for accessible SA / void volume via `window_fp_for_accessible=True`.
- One-sided `wr_ref − wr_cand` gate blocks annex over-bridging (XALDOE class; Mechanism D).
- Documentation: `docs/window_fingerprinting.org`; prototype `examples/window_types_lta96_report.py`.
- **Key commits:** `b10999a`, `ff06dbc`, `f723f00`, `6c36b13`.

### 6. CIF I/O, bonds, aromatics in Rust (`mgr_rust/src/cif.rs`, `bonds.rs`, `aromatic.rs`)

- Faster CIF parsing, covalent bond detection, aromatic cycle detection used by viz and pocket typing.
- **Key commits:** `dd19eaf`, `7589dcc`.

### 7. High-level plotting & workflow helpers

- `show_mof_plotly()` / `show_mof_threejs()` convenience APIs (`02e7966`, `51c356c`).
- `mofography/movies.py` — probe-smoothing and screen-door HTML movies.
- **`docs/plotting_guide.org`** — consolidated Plotly + Three.js / Figure3D guide (supersedes the older `docs/tutorial_figure3d.org` for new work; the tutorial org file remains as historical reference only).
- **`docs/upstream_features_from_mgr_cli_py.org`** — runnable walkthrough of features ported from `mgr_cli_py` (see caveat below).

---

## Not yet downstreamed to `mgr_cli_py` / web app

Many fixes in this fork live only in **`arung-northwestern/mofography` (`clean-master`)**. The Shi-group web app and CLI at **`~/repos/mgr_cli_py`** still track an older integration point and **do not** include:

- Optimized MIC dgrid (`dgrid_exp.rs`) and the `.T` cell-matrix fix in `dgrid_from_atoms_rust`
- Default `rust_det` watershed, skip-endchannel merge, window fingerprinting for accessible SA/vol
- One-sided window-radius FP gate (XALDOE class)
- Pocket-typing pipeline fixes from `81364a4`

`docs/upstream_features_from_mgr_cli_py.org` documents what was **ported from** `mgr_cli_py` **into** mofography, not what has been pushed back. Treat the public web app as **still buggy** for non-orthogonal cells and accessible textural properties until those paths are rebased onto `clean-master`.

---

## Bug fixes (correctness & embarrassing issues)

These are worth calling out explicitly — several affected **non-orthogonal cells**, **accessible property reporting**, or **viewer blank screens**.

### Distance grid: cell-matrix transpose + incorrect minimum image (MIC)

Two separate bugs affected the legacy brute-force dgrid path. Both are fixed in **`dgrid_from_atoms_rust`** (`mgr_rust/src/dgrid_exp.rs`, Python wrapper `mofography/rust_fun.py::dgrid_from_atoms_rust`). Full theory: `docs/dgrid_rust_mic_optimization.org`.

#### Bug A — wrong cell matrix passed to Rust (transpose)

ASE uses row-vector fractional coordinates: $\mathbf{r}_\mathrm{cart} = \mathbf{f}\,\mathbf{M}$ where $\mathbf{M} = \texttt{cellpar\_to\_cell}(\ldots)$ has **rows** = lattice vectors. The Rust kernel evaluates $\mathbf{r}_\mathrm{cart} = \mathbf{C}\,\Delta\mathbf{f}$ (column form), so we need $\mathbf{C} = \mathbf{M}^\top$.

**Wrong (pre-`4ec70b3`):** Python passed $\mathbf{M}$ instead of $\mathbf{M}^\top$, i.e. Rust computed $\mathbf{M}\,\Delta\mathbf{f}$ instead of $\mathbf{M}^\top\,\Delta\mathbf{f}$.

**Correct:** `_cell_matrix_for_rust()` returns `cellpar_to_cell(cellpar).T`.

For orthogonal cells $\mathbf{M}$ is diagonal so $\mathbf{M} = \mathbf{M}^\top$ and the bug is invisible. For **hexagonal** (e.g. NU-1000, $\gamma = 120°$), **monoclinic**, and **triclinic** MOFs, Cartesian distances were wrong by $\mathcal{O}(10\%)$ and displacement **directions** were wrong, breaking MIC image selection.

| Function | Transpose bug? | MIC algorithm | Status on `clean-master` |
|----------|----------------|---------------|--------------------------|
| `dgrid_from_atoms_rust` | **Fixed** (`_cell_matrix_for_rust`) | Wassenaar 7-image + orthogonal SIMD | **Use this** |
| `dgrid_from_atoms_aabb_rust` | **Fixed** | 27-replica AABB neighbor search (exact MIC via replication) | Safe |
| `dgrid_from_atoms_aabb_interp_rust` | **Fixed** | Same AABB family | Safe |
| `dgrid_from_atoms_cpu` / `dgrid_from_atoms` | N/A (Freud `AABBQuery` + `Box.from_matrix(cell.T)`) | Freud periodic box MIC | Generally safe |
| `dgrid_from_atoms_cpu_no_aabb` | Cell uses `complete_cell().T` ✓ | **Per-axis NINT only** — see Bug B | **Do not use** on skewed cells |
| `dgrid_from_atoms_gpu` | Cell uses `complete_cell().T` ✓ | **Per-axis NINT only** — see Bug B | **Do not use** on skewed cells |

#### Bug B — per-axis fractional NINT is not true MIC on skewed cells

The legacy brute-force path (Python `gpd_np` / CuPy `gpd_cp`, and pre-`dgrid_exp` Rust) used **per-component** wrapping:

$$
\Delta f_k = f_{\mathrm{grid},k} - f_{\mathrm{atom},k}, \qquad
\Delta f'_k = \Delta f_k - \mathrm{NINT}(\Delta f_k)
$$

$$
d = \left\| \mathbf{C}\,\Delta\mathbf{f}' \right\|_2 - R_i
$$

This minimizes $\|\Delta\mathbf{f}\|_\infty$ in fractional space, **not** the periodic Cartesian distance $\|\mathbf{C}\,\Delta\mathbf{f}\|_2$ under the metric $\mathbf{G} = \mathbf{C}^\top\mathbf{C}$. It is exact only when $\mathbf{C}$ is diagonal (cubic / tetragonal / orthorhombic).

**Correct (current `dgrid_exp.rs`):** After NINT, apply the Wassenaar certificate. If

$$
\left\| \mathbf{C}\,\Delta\mathbf{f}' \right\|_2 < d_\mathrm{max}
= \tfrac{1}{2}\min\bigl(\|\mathbf{a}\|, \|\mathbf{b}\|, \|\mathbf{c}\|\bigr)
$$

then NINT is provably exact; otherwise check the six face-adjacent images $\Delta\mathbf{f}' \pm \mathbf{e}_k$:

$$
d_\mathrm{MIC} = \min_{\mathbf{n} \in \{0, \pm\mathbf{e}_1, \pm\mathbf{e}_2, \pm\mathbf{e}_3\}}
\left\| \mathbf{C}\,(\Delta\mathbf{f}' + \mathbf{n}) \right\|_2 - R_i
$$

(Seven images total; matches unconditional 27-image search on benchmark MOFs.)

**Empirical failure:** rhombohedral MOF FOTNIN (60° angles) — LCD inflated from $\sim 33.6$ Å to $\sim 37$ Å because a pore voxel picked a wrong periodic image (28.4 Å vs 15.9 Å true MIC).

**Affected systems:** any **non-orthogonal** Bravais lattice — hexagonal ($\gamma \neq 90°$), monoclinic, triclinic, rhombohedral. Orthogonal cells are unaffected by Bug B (NINT is exact).

**Downstream impact:** wrong dgrid shifts LCD/PLD and local maxima placement, but **segmentation and pore graphs are less sensitive** than raw pore metrics: the watershed floods connected void domains and does not depend on exact Ångström-level dgrid values at every voxel. Expect larger changes in **LCD/PLD/PSD** than in **region count / connectivity topology** for typical MOFs.

**Still broken in this repo (verified in `mofography/compute.py`):**

- `dgrid_from_atoms_cpu_no_aabb` — lines 328–331: `(diff - np.around(diff))` then `np.dot(cell, diff.T)`
- `dgrid_from_atoms_gpu` — lines 957–959: same pattern with CuPy

Do **not** use these for production or parity checks on skewed cells; call `dgrid_from_atoms_rust` instead.

| Commit | Summary |
|--------|---------|
| `4ec70b3` | **Cell-matrix `.T` fix** for Rust dgrid wrappers (`dgrid_from_atoms_rust`, `dgrid_from_atoms_aabb_rust`, …). |
| `606c7c4` | Interim **27-image** brute-force MIC in legacy `lib.rs` `dgrid_from_atoms` (superseded by `dgrid_exp`). |
| `305a13a`–`96e5c10` | **`dgrid_exp.rs`**: tiered Wassenaar MIC + SoA/SIMD; promoted as default `mgr_rust.dgrid_from_atoms`. |

### Segmentation / grid / watershed (other)

| Commit | Summary |
|--------|---------|
| `31ef5d3` | Default **skip endchannel merge** (bad merges on some structures); remove experimental `rust_det_metric` from `mgr_rust`. |
| `4518ac8` | **by_pocket SA/V reporting** — partial PBC pockets dropped incorrectly. |
| `ee49495` | Spill guardrail for deterministic watershed backends. |

### Pocket typing & connections

| Commit | Summary |
|--------|---------|
| `81364a4` | Pocket typing pipeline fixes (7/7 benchmark pass). |
| `328bf22` | End-to-end snippet: 3 bugs + `pocket_regions_cluster` helper. |
| `f13e46e` | Import paths, `mask_thickness`, Three.js backend parity gaps. |
| `542edd1` | Reuse segmentation pocket PSDs; surface error handling. |
| `c3d052d` | Unit-cell connection typing and window plotting defaults. |
| `89efd32` | Window fractional coords; far-outside reclassification. |

### Accessible volume / window FP

| Commit | Summary |
|--------|---------|
| `6c36b13` | **FP over-bridge to annex apertures** (XALDOE) — cosine-only gate let sub-probe annex throats match open mouths; fixed with one-sided window-radius gate. |

### Three.js viewer (many small but user-visible bugs)

| Commit | Summary |
|--------|---------|
| `bd48ff8`, `1726cd3` | **Blank scene** — JS init order / duplicate declarations. |
| `b6547f6` | Blank-canvas TDZ bug with instanced framework. |
| `e0cbab8` | `_explodeOrig` before `updateFrameworkVisibility()`. |
| `d430670` | Fresnel/MatCap shaders broke mesh rendering (later reverted). |
| `94fa395` | Guest atoms rendered pitch-black. |
| `4be607d` | Periodic boundary bonds in flexible framework movies. |
| `5a477f0` | `trimesh_from_pyvista` mixed face-format meshes. |
| `9ea16dd` | Sidecar streaming cache when HTTP Range unsupported. |
| `adb0741` | Batched bond cylinders when bonds align with ±Z. |
| `057824f` | Wrong `step_size` default. |
| `51c356c` | Wrong Taubin smoothing defaults. |

---

## New & modified Python modules (`mofography/`)

| Path | Change | Notes |
|------|--------|-------|
| `threejs_viz.py` | **Added** | Standalone Figure3JS HTML export |
| `threejs_plotting.py` | **Added** | Lower-level Three.js trace builders |
| `viz_backend.py` | **Added** | Backend selection (Plotly vs Three.js) |
| `movies.py` | **Added** | HTML movie generation |
| `pocket_typing.py` | **Added** | Pocket-type assignment |
| `pocket_framework_env.py` | **Added** | Framework shell annulus for typing |
| `cif_io.py` | **Added** | CIF helpers |
| `vdw_radii.py` | **Added** | CRC/Bondi/Mendeleev radii tables |
| `export_mp4.py` | **Added** | MP4 export helper |
| `winding_percolation.py` | **Modified** | Window FP, canonical edges, DFSD helpers |
| `compute.py` | **Modified** | Accessible SA/vol + window FP defaults |
| `workflow.py` | **Modified** | Segmentation backends, radii, skip endchannel |
| `rust_fun.py` | **Modified** | New Rust entry points (dgrid_exp, watershed) |
| `visualize.py` | **Modified** | Plotly/Three.js integration |
| `analyze.py` | **Modified** | Bonds, aromatics hooks |

---

## Rust changes (`mgr_rust/`)

| Path | Change |
|------|--------|
| `mgr_rust/src/dgrid_exp.rs` | **Added** — optimized MIC distance grid |
| `mgr_rust/src/cif.rs` | **Added** |
| `mgr_rust/src/bonds.rs` | **Added** |
| `mgr_rust/src/aromatic.rs` | **Added** |
| `mgr_rust/src/lib.rs` | **Modified** — exports for new modules; default watershed is `rust_det` |

---

## Documentation & examples added

Not exhaustive; largest user-facing additions:

- `docs/window_fingerprinting.org` — accessible SA/vol mechanisms A–D
- `docs/plotting_guide.org` — **primary** Plotly + Three.js / Figure3D guide (replaces `tutorial_figure3d.org` for new work)
- `docs/upstream_features_from_mgr_cli_py.org` — features ported *from* `mgr_cli_py` (not downstream status)
- `docs/dgrid_rust_mic_optimization.org`, `docs/PLAN_fix_smoothing_and_cell_matrix.md`
- `docs/connections_and_windows.md`, `docs/pore_graph.md`, `docs/region_segmentation.md`
- `examples/window_types_lta96_report.py` — window-type clustering prototype

---

## Recommended fixes for legacy `compute.py` dgrid paths

The per-axis NINT bug in `dgrid_from_atoms_cpu_no_aabb` and `dgrid_from_atoms_gpu` should be addressed in one of these ways (ordered by preference):

1. **Delegate to Rust (simplest, recommended).** Replace the body of both functions with a call to `dgrid_from_atoms_rust(ase_atoms, spacing, radii=...)`. Keeps one MIC implementation; removes ~200 lines of duplicated Mendeleev tables.

2. **Deprecate and redirect.** Emit `warnings.warn(..., DeprecationWarning)` and call `dgrid_from_atoms_rust`. Leave stubs for API compatibility.

3. **Port Wassenaar MIC into Python/CuPy** (only if GPU path must stay on device). Extract the seven-image loop from `dgrid_exp.rs::mic_min_dist_skewed` into a shared helper used by `gpd_np` and the CuPy kernel:
   - After NINT wrap, compute $d_\mathrm{max} = \tfrac{1}{2}\min(\|\mathbf{a}\|,\|\mathbf{b}\|,\|\mathbf{c}\|)$ from column norms of $\mathbf{C}$
   - If $\|\mathbf{C}\Delta\mathbf{f}'\|_2 \ge d_\mathrm{max}$, evaluate the six face offsets and take the minimum
   - Add regression tests against `dgrid_from_atoms_rust` on NU-1000 and a triclinic cell

4. **Leave `dgrid_from_atoms_cpu` / `dgrid_from_atoms` on Freud AABB** — these use `freud.box.Box` periodic queries and are not subject to the per-axis NINT bug; optionally switch their default backend to Rust AABB for consistency.

5. **Unify cell helpers.** Replace scattered `complete_cell(...).T` with `_cell_matrix_for_rust()` from `rust_fun.py` (or move that helper to `mofography/helper.py`) so transpose convention is documented in one place.

---

## Full commit list (`1444c2e..2d63da9`, 107 commits)

Newest first. Use `git show <hash>` on `master` for patches.

```
2d63da9 docs(changelog): set handoff tip to c1e38db (106 commits)
c1e38db docs(changelog): update fork handoff for master @ 840ee9a
840ee9a chore: untrack agent state, graphify cache, and root clutter
b487af5 fix(watershed): drop failed metric/blended backends; default rust_det
6c36b13 fix(winding): one-sided window FP gate; add docs and metric v1 watershed
f723f00 feat(compute): default window FP for accessible SA and void volume
ff06dbc feat(winding): add window-fingerprint accessible region expansion
31ef5d3 fix(segmentation): default skip endchannel merge; drop rust_det_metric
2e9d8f8 feat(watershed): add blended-priority metric watershed (rust_det_blended)
5c37652 fix(watershed): wire rust_det_metric straight through main pipeline
c317e32 feat(watershed): metric merged flood for rust_det_metric regions_pbc
36539d0 feat(watershed): add rust_det_metric backend for skewed cells
b10999a feat: CRC/Bondi VDW radii, canonical edges, DFSD threshold & window fingerprint
2b928fa docs: add org-roam note on optimized Rust dgrid MIC
96e5c10 refactor(rust): promote optimized MIC dgrid as default dgrid_from_atoms
6773f63 perf(rust): SoA + SIMD + Wassenaar 7-image MIC for dgrid_from_atoms_exp
305a13a feat(rust): add dgrid_from_atoms_exp with tiered exact MIC
e415f03 feat(rust): default Mendeleev VDW radii for dgrid Rust wrappers
606c7c4 fix(rust): use 27-image MIC in dgrid_from_atoms for non-orthogonal cells
e0b82e2 docs(plan): add dgrid_from_atoms_exp MIC optimization plan
ee49495 feat(winding): add spill guardrail for deterministic watershed backends
2a2ebfb feat(watershed): add deterministic rust_det_max backend
256a14e docs: plotting guide, CSD reference, and upstream feature notes
f5f1c7b feat(movies): expand probe-smoothing and screen-door HTML movie generation
9186ef8 feat(viz): Figure3JS helpers, percolation scene, and window plot API
542edd1 fix(pocket-typing): reuse segmentation pocket PSDs and surface errors
c3d052d fix: unit-cell connection typing and window plotting defaults
8dd08ce docs(dev): cell/smoothing fix plan, AABB benchmark notes, and viz examples
94f3668 docs: add Figure3D tutorial and upstream mgr_cli feature notes
424e8f5 feat(viz): add threejs_viz module for standalone Figure3JS HTML export
7589dcc perf(analyze): faster covalent bonds and aromatic cycle detection
89efd32 fix(viz): window fractional coords and far-outside reclassification
aa4a221 refactor(rag): separate CSD centroid histograms from pocket PSD node keys
4ec70b3 fix(rust): correct dgrid cell matrix and physical smoothing for non-orthogonal cells
b6547f6 feat: instanced framework representation + fix blank-canvas TDZ bug
e0cbab8 fix: move _explodeOrig before updateFrameworkVisibility() init call
5088718 port crystallographic explode plane controls from mgr_cli_py
057824f fix step_size default + port merged framework opacity slider
12496e7 docs: update §7 — show() usage, corrected taubin defaults, quality notes
51c356c fix taubin defaults + add show() + add taubin to plot_surface_dgrid
e03af5f docs: replace §7 with show_mof_plotly/threejs snippets
02e7966 add show_mof_plotly/threejs + fix 5 supporting bugs
328bf22 fix pocket typing end-to-end snippet: 3 bugs, add pocket_regions_cluster helper
2fe9506 docs: update upstream_features_from_mgr_cli_py.org for 81364a4 fixes
81364a4 fix pocket typing pipeline: 7/7 benchmark pass
f13e46e fix upstream gaps: pocket typing import, mask_thickness, Three.js backend parity
18f7d5b docs: add upstream_features_from_mgr_cli_py.org — org-mode runnable walkthrough
dd19eaf upstream from mgr_cli_py: Rust CIF+bonds+aromatic, pocket typing, window improvements
28eacbc checkpoint: pre-upstream-from-mgr_cli_py baseline
adb0741 Fix batched bond cylinders when bonds align with ±Z.
aeeb16c Add aromatic ring rendering to Three.js viewer
3e21721 Ignore Playwright browsers, quick-start HTML output, and huge tutorial blend
8ddd582 Update Figure3D tutorial with movie builders, sidecar streaming, and Plotly comparison
9ea16dd Fix sidecar streaming to cache full response text once when Range unsupported
2a6eaa4 Implement per-frame HTTP Range streaming for sidecar movies
22c6eb7 Add binary sidecar packing for numerical frame data (2x smaller sidecars)
6a73d3c Implement lazy sidecar loading for large movie HTMLs
94fa395 Fix guest atom pitch-black coloring in Three.js viewer
5a477f0 Fix trimesh_from_pyvista for mixed face-format meshes
4be607d Fix periodic boundary bonds in flexible framework movies
95411a4 Add pore radius (dgrid) and region label to maxima hover tooltip
984a121 Add built-in GIF recording for movie HTMLs
6f3e1b1 Add video recording (WebM) to interactive viewer
0640d0e Fix hover on local maxima (THREE.Points) traces
55cc8ad Fix cell_pbc to only follow its own explode-group traces
5ca6326 Fix cell-copy visibility to track actual trace visibility
0052dfa Fix exploded cell-copy visibility for all group toggles
d25f82f Exploded view: group-specific cell copies tied to legend checkboxes
b1f437e Exploded view: include stick and spacefill variants in framework shift
5bddfa9 Fix exploded view surface trace IDs
c58b444 Multi-directional exploded view: surfaces split in Y
063d499 Exploded view: add stationary cell_center copy at origin
ba11ce9 Exploded view: bidirectional separation
9d5c0c5 Fix hover on stacked surfaces: pick topmost rendered mesh
6843c4f Set camera from cell dims (1.5x la,lb,lc); add surface stacking note
bbdfdf0 Regenerate demo HTMLs with latest viewer template
62bb5d6 Revert depthWrite changes; enhance axis views and camera distance
ff8bb61 Enhance interactive viewer: ortho toggle, axis views, camera distance, explode range, depthWrite fix
1cad230 Fix viewer: unified framework opacity, single exploded cell, vdw radii
bd48ff8 Fix blank scene: remove duplicate cellObj declaration, init framework visibility
1726cd3 Fix blank scene: move _explodeGroups before cell copy creation
19df805 Add space-filling representation, per-mesh opacity, and exploded cell copies
5db2e43 Add stick representation toggle for framework mesh
5ef8ee6 Fix controls: mesh opacity targets all meshes, redesigned exploded view, default camera (1.5,1.5,1.5), remove atom labels and fw scale
c90d9fd Add dynamic HUD controls: mesh opacity, framework scale, point size, speed, lights, bg, PNG res, measure tool, atom labels, reset cam, auto-rotate, exploded view + A/B/C axis labels
5368769 Fix orthographic rendering: point visibility, lighting contrast, CSS filter
64554ad Revert ortho material mods, minimal light boost only — still pale in ortho
f345a4c Remove shader system (fresnel/matcap), restore clean default rendering
d430670 Fix shaders: fresnel/matcap were breaking mesh rendering
d3f632f Snappier rotation (dampingFactor 0.08->0.25)
981599c Rich hover, cell axis labels, light popup, MatCap/Fresnel shaders
9b6674d Right-side panel, title, PNG export, axis labels, white mode
5b873cb HUD bottom-centered horizontal layout + mobile responsive + Plotly trace API
9089f98 Surface region hover + ortho scatter depth fix
fcc0d33 Shiny surfaces: clearcoat + emissive for Plotly-like appearance
b73ffb2 Add appearance/colorscale APIs + Figure3D tutorial
b8862a7 Tone down ortho lighting for Pastel surfaces
f66ee05 Surface colors: HSV wheel -> px.colors.qualitative.Pastel
49fde99 Periodic pore graph hidden by default via default_off flag
7e51dc9 Use px.colors.qualitative.Pastel palette, 4x ortho point boost
68929f9 Fix ortho lighting/coloring and scatter point sizes
62adc89 Figure3D: Plotly-like Three.js figure API for static + movie scenes
e1e0b7d feat: add watershed_connectivity param to do_periodic_segmentation_rust for skimage backend
4518ac8 Fix by_pocket SA/V reporting: keep partial PBC pockets; cluster from threshold components
60acacc Add movie HTML module with flexible framework and docs
5c6f9d4 Add Three.js viewer script with framework hover matching Plotly traces.
df1a1b3 feat(devtools): add camera and background controls to animation script
```

---

## Fork update workflow (GitHub UI)

**Status: ready now.** Northwestern `master` and `clean-master` are pushed and identical at `2d63da9`. No need to wait for further commits before opening the PR.

### Option A — Cross-fork pull request (recommended for Shi group)

1. Open: [Compare across forks](https://github.com/Shi-Research-Group/mofography/compare/master...arung-northwestern:mofography:master)  
   - **base repository:** `Shi-Research-Group/mofography`, branch **`master`**  
   - **head repository:** `arung-northwestern/mofography`, branch **`master`**  
   (If the link 404s, use **New pull request** → **compare across forks** and set those branches manually.)
2. Title suggestion: *Sync from arung-northwestern/mofography master (107 commits, through 2d63da9)*  
3. Attach or link this changelog: `docs/CHANGELOG_arung_vs_Shi-Research-Group_2026-06-11.md`
4. Review conflicts — expect few on library code; possible noise on Shi-side merge commits. Resolve by keeping Northwestern feature commits for `mofography/`, `mgr_rust/`, `docs/`.
5. **Merge** when CI/review passes (merge commit or squash per Shi group preference).

### Option B — Northwestern fork of Shi (if you maintain a GitHub fork)

On **your** fork of `Shi-Research-Group/mofography`, use **Sync fork** → **Update branch** only if GitHub offers it from `arung-northwestern`. Often cross-fork PR (Option A) is clearer when both repos are org/user remotes you control.

### After merge — verify locally

```bash
pytest mofography/tests/
cd mgr_rust && maturin develop --release
# Spot-check dgrid on a hexagonal MOF: dgrid_from_atoms_rust (not dgrid_from_atoms_cpu_no_aabb)
```

- Skim `docs/window_fingerprinting.org` — accessible SA/vol defaults (`window_fp_for_accessible=True`).
- Regenerate graphify locally: `graphify update .` (`graphify-out/` is no longer in git).
- **Do not assume `mgr_cli_py` / web app** pick up these fixes until rebased — see [Not yet downstreamed](#not-yet-downstreamed-to-mgr_cli_py--web-app).

---

## Reproduce this report locally

```bash
git fetch origin master
git fetch https://github.com/Shi-Research-Group/mofography.git master:shi/master

git rev-list --left-right --count shi/master...origin/master
# → 40  107  (Shi merge commits vs Northwestern feature commits)

git log shi/master..origin/master --oneline --no-merges | wc -l
git diff --stat shi/master..origin/master -- mofography mgr_rust docs examples scripts
```

---

*Handoff to Shi Research Group. Northwestern: `master` = `clean-master` @ `2d63da9` (2026-06-11).*

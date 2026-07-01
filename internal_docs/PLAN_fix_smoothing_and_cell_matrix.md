¶docs/PLAN_fix_smoothing_and_cell_matrix.md#47A7
1:# Plan: Fix Smoothing Anisotropy + Cell Matrix Transposition Bugs
2:
3:**Date**: 2026-05-31
4:**Status**: Pending Review (REVISED after senior review caught inverted Bug 2 analysis)
5:**Scope**: `mgr_rust/src/lib.rs`, `mofography/rust_fun.py`, `mgr_rust/test_compute_pore_sizes.py`
6:
7:---
8:
9:## Critical Convention Note
10:
11:ASE uses **row-vector** convention for fractional coordinates:
12:```
13:cart = frac @ M    (where frac is a 1×3 row vector)
14:```
15:This is equivalent to `cart = M.T @ frac` (column vector form).
16:
17:The Rust code computes:
18:```rust
19:cart_x = c00*dx + c01*dy + c02*dz   // = C[0,:] · frac = (C @ frac)[0]
20:```
21:This is `C @ frac`. So for correctness: `C @ frac` must equal `frac @ M = M.T @ frac`, therefore **C = M.T**.
22:
23:**Conclusion**: Passing `cellpar_to_cell(...).T` to Rust is **CORRECT**. Locations WITHOUT `.T` are the bugs.
24:
25:---
26:
27:## Bug 1: Smoothing Radius Unit Mismatch + No Anisotropy
28:
29:### Status: CONFIRMED CORRECT by senior review
30:
31:### What's wrong
32:
33:The `detect_maxima_selective_smoothing*` functions in `lib.rs` apply a Gaussian smoothing kernel where:
34:
35:1. `max_smooth_radius` is treated as Gaussian σ in **voxel units** (not Å)
36:2. The distance is computed as `dist_sq = di² + dj² + dk²` — isotropic in voxel-offset space
37:3. No `spacing` or cell geometry information is available to these functions
38:
39:This means:
40:- The Python docstring (`rust_fun.py:240`) says `max_smooth_radius` is in **Å**, but the Rust code uses it as **voxels** — a unit mismatch
41:- For non-orthogonal cells, voxels are skewed parallelepipeds (not cubes), so a 1-voxel offset along axis i ≠ 1-voxel offset along axis j in physical distance
42:- The grid axes are not perpendicular, so `di² + dj² + dk²` does not equal the physical squared distance
43:
44:### Evidence
45:
46:**Rust kernel** (`lib.rs:3452-3456`, identical pattern in all three functions):
47:```rust
48:let dist_sq = (di * di + dj * dj + dk * dk) as f32;  // voxel² — NOT Å²
49:let sigma_sq = smooth_radius * smooth_radius;          // voxel² — NOT Å²
50:let weight = (-0.5 * dist_sq / sigma_sq).exp();
51:```
52:
53:**Kernel size** (`lib.rs:3440`):
54:```rust
55:let kernel_size = (smooth_radius * 3.0).ceil() as usize;  // in voxels
56:```
57:
58:**No spacing parameter** in any of the three function signatures (`lib.rs:3383`, `3614`, `3857`).
59:
60:**Python docstring** (`rust_fun.py:240`):
61:```python
62::param max_smooth_radius: Maximum radius (Å) used when smoothing maxima
63:```
64:
65:**Python caller** (`rust_fun.py:352-362`) — does NOT pass spacing:
66:```python
67:merged_maxima, unmerged_maxima, raw_maxima = (
68:    mgr_rust.detect_maxima_selective_smoothing_parallel(
69:        dgrid,
70:        preserve_threshold=preserve_threshold,
71:        smooth_threshold=smooth_threshold,
72:        max_smooth_radius=max_smooth_radius,
73:        h=h,
74:        threshold_abs=threshold_abs,
75:        min_distance=min_distance,
76:        cell=cell,  # passed with .T — CORRECT for merge step
77:    )
78:)
79:```
80:
81:### Affected Functions
82:
83:| Function | File:Lines |
84:|---|---|
85:| `detect_maxima_selective_smoothing` | `lib.rs:3383-3470` |
86:| `detect_maxima_selective_smoothing_parallel` | `lib.rs:3614-3856` |
87:| `detect_maxima_selective_smoothing_parallel2` | `lib.rs:3857-4105` |
88:
89:### Proposed Fix
90:
91:**Approach: Per-axis anisotropic Gaussian (Option B)**
92:
93:1. Add a `voxel_lengths: Option<(f32, f32, f32)>` parameter to all three Rust functions (lx, ly, lz in Å)
94:2. When `voxel_lengths` is provided:
95:   - Interpret `max_smooth_radius` as Å (matching the docstring)
96:   - Compute per-axis kernel bounds: `kernel_i = ceil(3.0 · sigma / lx)`, etc.
97:   - Compute physical distance: `dist_sq = (di·lx)² + (dj·ly)² + (dk·lz)²`
98:   - Use `sigma` directly as Å
99:3. When `voxel_lengths` is `None` (backward compat): fall back to current voxel-unit behavior
100:
101:**Kernel change** (from):
102:```rust
103:let kernel_size = (smooth_radius * 3.0).ceil() as usize;
104:for di in -(kernel_size as isize)..=(kernel_size as isize) {
105:    for dj in -(kernel_size as isize)..=(kernel_size as isize) {
106:        for dk in -(kernel_size as isize)..=(kernel_size as isize) {
107:            let dist_sq = (di * di + dj * dj + dk * dk) as f32;
108:            let sigma_sq = smooth_radius * smooth_radius;
109:            let weight = (-0.5 * dist_sq / sigma_sq).exp();
110:```
111:
112:**Kernel change** (to):
113:```rust
114:// Per-axis kernel bounds (in voxels, derived from physical σ)
115:let (ki, kj, kk) = if let Some((lx, ly, lz)) = voxel_lengths {
116:    let ki = (smooth_radius * 3.0 / lx).ceil() as usize;
117:    let kj = (smooth_radius * 3.0 / ly).ceil() as usize;
118:    let kk = (smooth_radius * 3.0 / lz).ceil() as usize;
119:    (ki, kj, kk)
120:} else {
121:    let k = (smooth_radius * 3.0).ceil() as usize;
122:    (k, k, k)
123:};
124:
125:for di in -(ki as isize)..=(ki as isize) {
126:    for dj in -(kj as isize)..=(kj as isize) {
127:        for dk in -(kk as isize)..=(kk as isize) {
128:            let dist_sq = if let Some((lx, ly, lz)) = voxel_lengths {
129:                let dx = di as f32 * lx;
130:                let dy = dj as f32 * ly;
131:                let dz = dk as f32 * lz;
132:                dx * dx + dy * dy + dz * dz
133:            } else {
134:                (di * di + dj * dj + dk * dk) as f32
135:            };
136:            let sigma_sq = smooth_radius * smooth_radius;
137:            let weight = (-0.5 * dist_sq / sigma_sq).exp();
138:```
139:
140:**Python-side changes** (`rust_fun.py`):
141:- `do_periodic_segmentation_rust` (line 352-362): compute per-axis voxel lengths and pass them
142:- Per-axis voxel lengths from `M = cellpar_to_cell(...)` (without .T, since we want lattice vector magnitudes):
143:  `lx = ||M[0] / (nx-1)||`, `ly = ||M[1] / (ny-1)||`, `lz = ||M[2] / (nz-1)||`
144:- Note: the cell passed to the Rust smoothing function for the merge step STILL uses `.T` — that's correct and separate from voxel_lengths
145:
146:### Limitations
147:
148:- This approach neglects off-diagonal metric terms (triclinic skew). For most MOFs (orthorhombic or near-orthorhombic) this is negligible.
149:- A full metric tensor approach (`dist_sq = dv^T G dv` where `G = M_rac^T M_rac`) would be more correct for highly skewed cells but adds complexity.
150:
151:---
152:
153:## Bug 2: Cell Matrix — 3 Locations MISSING `.T` (REVERSED from initial analysis)
154:
155:### What's wrong
156:
157:ASE uses row-vector convention: `cart = frac @ M`. The Rust code computes `C @ frac`. For these to match: `C = M.T`.
158:
159:**3 locations in `rust_fun.py` pass `cellpar_to_cell(...)` WITHOUT `.T`** — this causes the Rust code to compute `M @ frac` instead of the correct `M.T @ frac`. For non-orthogonal cells, this produces wrong Cartesian coordinates.
160:
161:### Proof (numerical)
162:
163:Triclinic cell: a=10, b=12, c=14, α=80°, β=90°, γ=100°
164:
165:```
166:M = cellpar_to_cell(cellpar)   # rows = lattice vectors
167:
168:frac @ M = [2.166, 5.961, 6.890]   ← ASE expected result
169:M.T @ frac = [2.166, 5.961, 6.890] ← SAME (this is what .T achieves in Rust)
170:M @ frac = [3.000, 4.102, 7.878]   ← WRONG (this is what no .T gives in Rust)
171:```
172:
173:For orthogonal cells, M is diagonal so `M == M.T` and the bug is invisible.
174:
175:### Locations That Need `.T` ADDED
176:
177:| # | File | Line | Function | Current code | Fix |
178:|---|---|---|---|---|---|
179:| 1 | `mofography/rust_fun.py` | 34 | `dgrid_from_atoms_rust` | `cellpar_to_cell(...)` | Add `.T` |
180:| 2 | `mofography/rust_fun.py` | 65 | `dgrid_from_atoms_aabb_rust` | `cellpar_to_cell(...)` | Add `.T` |
181:| 3 | `mofography/rust_fun.py` | 94 | `dgrid_from_atoms_aabb_interp_rust` | `cellpar_to_cell(...)` | Add `.T` |
182:
183:### Locations That Are Already Correct (no change needed)
184:
185:| File | Line | Function | Code |
186:|---|---|---|---|
187:| `mofography/rust_fun.py` | 145 | `egrid_from_atoms_rust` | `cellpar_to_cell(...).T` ✓ |
188:| `mofography/rust_fun.py` | 350 | `do_periodic_segmentation_rust` → merge | `cellpar_to_cell(...).T` ✓ |
189:| `mofography/rust_fun.py` | 955 | `do_periodic_segmentation_rust` → window_connections | `cellpar_to_cell(...).T` ✓ |
190:| `mofography/rust_fun.py` | 1202 | `compute_pore_sizes_uniform_rust` | `cellpar_to_cell(...).T` ✓ |
191:| `mofography/rust_fun.py` | 1303 | `compute_pore_sizes_pocketwise_rust` | `cellpar_to_cell(...).T` ✓ |
192:| `mofography/rust_fun.py` | 1466 | `compute_psds_pocketwise_accessible` | `cellpar_to_cell(...).T` ✓ |
193:
194:### Impact
195:
196:- **The dgrid is the foundational computation** — wrong dgrid → wrong everything downstream
197:- For NU-1000 (hexagonal, γ=120°): individual-axis distance errors ~11-13%
198:- For monoclinic cells (α=80, β=110, γ=100): combined distance error ~10.5%
199:- Even when distances happen to be similar, the displacement **direction** is wrong, affecting minimum-image selection
200:- For orthogonal cells: no change (M == M.T)
201:
202:### Why this wasn't caught
203:
204:- `dgrid_from_atoms` (brute-force exact) and `dgrid_from_atoms_aabb` (R-tree accelerated) both have the bug, so their outputs agree with each other — but both are wrong for non-orthogonal cells
205:- Most test cases use orthorhombic/hypercubic cells where `M == M.T`
206:- The comment at `rust_fun.py:146` ("3×3 column-vector matrix") is actually CORRECT — it accurately describes what `.T` produces
207:
208:---
209:
210:## Implementation Sequence
211:
212:### Phase 1: Cell Matrix `.T` Fix (dgrid — highest priority)
213:
214:1. **Add `.T`** to lines 34, 65, 94 in `rust_fun.py` (the 3 dgrid functions)
215:2. Build: `cd mgr_rust && cargo build --release`
216:3. Verify: for orthogonal cells, dgrid values should be unchanged
217:4. Verify: for non-orthogonal cells, dgrid values will change — spot-check against Python reference
218:
219:### Phase 2: Smoothing Anisotropy Fix
220:
221:1. Add `voxel_lengths: Option<(f32, f32, f32)>` parameter to all three `detect_maxima_selective_smoothing*` functions in `lib.rs`
222:2. Modify the Gaussian kernel to use per-axis voxel lengths when provided
223:3. Update Python wrapper in `rust_fun.py` to compute and pass voxel lengths
224:4. Note: the cell matrix passed for the merge step keeps its `.T` (correct)
225:5. Build and test
226:
227:### Phase 3: Verification
228:
229:1. Run end-to-end on orthogonal cell (IRMOF-1) — results should match
230:2. Run end-to-end on non-orthogonal cell (NU-1000, hexagonal γ=120°) — results will change, verify they're physically reasonable
231:3. Cross-check dgrid values against a known Python reference for a non-orthogonal cell
232:
233:---
234:
235:## Open Questions for Reviewer
236:
237:1. **Smoothing approach**: Is per-axis anisotropic (Option B) sufficient, or should we use a full metric tensor for triclinic cells?
238:2. **Default behavior**: Should `max_smooth_radius` default remain 2.0 (now interpreted as Å instead of voxels)? At 0.5 Å spacing this means ~4 voxels σ instead of the current ~2 voxels.
239:3. **Backward compatibility**: The `.T` fix on dgrid will change results for ALL non-orthogonal cells. The smoothing fix (when `voxel_lengths` provided) will also change results. Is this acceptable?
240:4. **Test coverage**: Should we add a dedicated non-orthogonal test case to the test suite?

---

## Reviewer Feedback (2026-05-31)

**Overall verdict:** Plan is sound and ready to implement. Bug 2 reversal and convention note match the codebase. Phase ordering (dgrid first) is correct. Address items below before or during implementation.

### Confirmed correct

- **C = M.T** matches Rust `cart = C @ frac` vs ASE `frac @ M` (verified in `lib.rs` dgrid and merge paths).
- **Three buggy call sites** are real: `rust_fun.py` lines 34, 65, 94 pass `cellpar_to_cell(...)` without `.T`; egrid, segmentation, and PSD paths already use `.T`.
- **Bug 1** is real: Rust docstrings say "grid units"; Python docstring says Å; kernel uses `di² + dj² + dk²` in voxel space.
- **Senior review reversal** is well documented — prevents repeating the inverted Bug 2 analysis.

### Recommended additions

#### 1. Single helper for cell matrix (high value)

Add once in `rust_fun.py`:

```python
def _cell_matrix_for_rust(atoms):
  """ASE row-vector cell M; Rust expects C = M.T so C @ frac == frac @ M."""
  return cellpar_to_cell(atoms.cell.cellpar()).T
```

Use in all Rust wrappers (dgrid, egrid, segmentation, PSD). Prevents future `.T` drift. Plan lists 6 correct sites + 3 fixes — optional refactor but strongly recommended.

#### 2. Phase 3 reference path — be explicit

Plan says "Python reference" — specify:

- `mofography/compute.py` uses `complete_cell(get_cell()).T` (e.g. lines 146, 774, 6863), equivalent to `cellpar_to_cell(...).T` for standard ASE cells.
- Golden test: one triclinic or hexagonal CIF, one grid point, compare Rust dgrid vs `dgrid_from_atoms_cpu` or manual `||cell @ wrap(Δfrac)|| - r`.

#### 3. Sync Rust docstrings with Bug 1 fix

When `voxel_lengths` is provided, update misleading comments in `lib.rs` (e.g. lines 3378, 3389, 3620, 3863: "grid units" / "Maximum smoothing radius in grid units"). Otherwise doc/API mismatch returns immediately.

#### 4. Decide `max_smooth_radius` default now (answers Open Question 2)

| Mode | Effective behavior |
|------|-------------------|
| Today | ~2 voxel σ (spacing-dependent) |
| After fix (2.0 Å) | At 0.5 Å spacing ≈ 4 voxel support — **much wider** smoothing |

**Recommendation:** Keep default `2.0` as Å (matches Python docstring intent). Document as intentional behavior change in CHANGELOG. Optionally retune (e.g. 1.0 Å) after NU-1000 regression.

#### 5. Always pass `voxel_lengths` from production Python

`Option::None` fallback is fine for API backward compat, but `do_periodic_segmentation_rust` should **always** compute and pass voxel lengths — avoid silent voxel-unit smoothing in production.

#### 6. Metric tensor vs Option B (answers Open Question 1)

| Cell type | Recommendation |
|-----------|----------------|
| Orthorhombic / near-ortho MOFs | Per-axis Option B sufficient |
| High skew (γ far from 90°, monoclinic/triclinic) | Cross terms in `G = MᵀM` matter for smoothing kernel, not only merge |

For NU-1000 (hexagonal), run one comparison: Option B vs full metric `dist² = Δfᵀ G Δf` with `Δf = (di/(nx-1), dj/(ny-1), dk/(nz-1))` and `G = M Mᵀ` (M rows = lattice vectors). If delta is small, ship Option B and document limitation.

#### 7. Tests — mandatory (answers Open Question 4: yes)

| Test | Assertion |
|------|-----------|
| Orthogonal (IRMOF-1) | dgrid + segmentation unchanged within rtol |
| Non-orthogonal (NU-1000 or triclinic) | dgrid changes vs pre-fix; matches ASE / `complete_cell().T` reference |
| Smoothing | With `voxel_lengths`, σ in Å gives expected kernel width vs spacing |

Add to `mgr_rust/test_compute_pore_sizes.py` or new `test_cell_convention.py`.

#### 8. Fix all three smoothing variants

`detect_maxima_selective_smoothing_parallel2` is exported but has no Python caller — still patch for parity (low cost).

#### 9. Misleading Rust comments (same PR)

`lib.rs` ~2265 says "lattice vectors as columns" but code uses row-dot `C @ frac`. Fix in same PR to reduce future `.T` confusion.

#### 10. Breaking change communication (answers Open Question 3: yes, acceptable)

Treat as **bugfix**, not regression. Require CHANGELOG entry; consider semver minor if package promises numerical stability across cell types.

### Reviewer answers to open questions

| # | Question | Answer |
|---|----------|--------|
| 1 | Option B vs metric tensor? | Option B for most MOFs; add metric tensor if triclinic golden test fails tolerance |
| 2 | Default `max_smooth_radius`? | Keep 2.0 Å; document behavior change |
| 3 | Backward compatibility? | Yes — intentional fix; CHANGELOG + non-ortho tests |
| 4 | Dedicated non-orthogonal tests? | Yes — required |

### Minor nits

- Plan date `2026-05-31` may be a typo (2025?).
- `voxel_lengths` from `||M[i]|| / (n_i - 1)` is correct for linspace `-0.5…0.5` grid convention used in Rust.
- Verify `complete_cell` vs `cellpar_to_cell` equivalence once on a triclinic cell in tests.

### Implementation summary

| Item | Rating |
|------|--------|
| Diagnosis | Strong |
| Bug 2 fix | Approve |
| Bug 1 fix | Approve with doc sync + default σ decision |
| Sequencing | Approve |
| Gaps | Helper fn, golden tests, Rust doc sync, CHANGELOG |

**Approve implementation** after locking default-σ policy and test plan. Phase 1 alone is high-value and can merge independently of Phase 2.
241:
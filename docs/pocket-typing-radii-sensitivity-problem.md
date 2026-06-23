# Pocket typing under VDW radii choice: problem statement

**Repo:** `mofography` (MOF void segmentation + pocket typing)  
**Audience:** External technical review (LLM or human). Please propose *new* approaches; do not assume the directions listed under “Attempted directions” are correct or complete.  
**Date:** June 2026

---

## What we are trying to do

Assign each void pocket (cage, window, channel segment) in a periodic MOF structure to a **pocket type** — an equivalence class of void environments that chemists would treat as the same adsorption site / same cage family.

**Desired outputs:**

- Integer type label per pocket (PBC representative ID)
- Stable across **reasonable** numerical choices where possible
- Validated on a small benchmark set (see below)

**Non-goals for this document:**

- Perfect agreement with every external porosity code’s radii convention
- Unsupervised typing with no physical meaning

---

## Pipeline (current production path)

```
CIF + VDW radii  →  dgrid (3D distance field)  →  periodic watershed segmentation
                 →  pocket list + fingerprints  →  env_cluster typing
```

### `dgrid`

At each grid voxel, store the radius of the largest probe sphere centered there that does not overlap framework atoms (atoms inflated by the chosen VDW set). Computed in Rust (`dgrid_from_atoms_rust`). **The entire void geometry is a function of the radii table.**

Grid spacing default for typing benchmarks: **0.5 Å**.

### Segmentation

Watershed on `-dgrid` yields region labels; periodic copies are grouped into PBC pockets. Connections / saddles between pockets come from winding percolation (`collect_periodic_edges_canonical`).

### Fingerprints (per pocket)

Multiple channels exist; production `env_cluster` primarily uses:

| Channel | What it measures | Used in env_cluster? |
|---------|------------------|----------------------|
| `max_r` | Max `dgrid` on pocket fragment (Å) | **Yes — seeds ordered bands** |
| Pocket PSD | Histogram of probe-occupiable sizes (Rust sampler) | Indirectly via `d1_A` in older `shape_refined` path; in merge cosines for env_cluster |
| Framework SH | Spherical harmonics of framework atoms near void | Merge cosine similarity |
| Connectivity SH | SH of directions to periodic neighbors | Partially / audit |
| CSD | Histogram of \|r_i − r_c\| for void voxels, PBC circmean centroid | Computed (`add_csd_to_rag_with_pbc`) but **not** primary typing gate today |

### `env_cluster` (default typing method)

1. Sort pockets by `max_r`.
2. **Seed bands:** split when consecutive `max_r` gap **> 0.15 Å** (`ENV_CLUSTER_INITIAL_BAND_GAP_A`).
3. **Merge adjacent bands** when fingerprint cosine ≥ ~0.90 (and gap ≤ 0.5 Å), with a looser cosine for large gaps.

No split pass after merge in the basic env_cluster — merge-only correction of over-segmentation.

---

## Benchmark materials (expected pocket type counts)

Validated informally under **Bondi VDW + 0.5 Å spacing** with `env_cluster`:

| Material | CIF (in repo) | Expected types | Notes |
|----------|---------------|----------------|-------|
| NU-1000 | `examples/TestCIFs/NU-1000_S1.cif` | **3** | Micropores `{2,3,8,9}` must be **one** type |
| ZIF-20 | `~/repos/RASPA2/structures/mofs/cif/ZIF-20.cif` | **6** (Bondi) | Mid-window families ~2.56 Å vs ~2.74 Å `max_r` must stay **separate** |
| UiO-66 | `examples/TestCIFs/UiO-66-RASPA.cif` | **2** | Tetrahedral vs octahedral cages; similar accessible size, **different shape** |
| IRMOF-1 | `examples/TestCIFs/P1-IRMOF-1.cif` | **2** | Two cage-size families |
| Cu-BTC | `examples/TestCIFs/P1-Cu-BTC.cif` | **3** | Window / cage / large cage |

**Sensitivity target:** NU-1000 is the canonical failure when switching VDW table from Bondi to CRC or Alvarez while holding 0.5 Å spacing.

---

## The crux: NU-1000 micropores (measured @ 0.5 Å)

### What does **not** change across radii sets

- Watershed **pocket count:** 11
- **Pocket representative IDs:** same 1…11
- Micropore **identity** (still pockets 2, 3, 8, 9)
- No wholesale relabeling or domain fragmentation

### What **does** change

**Zirconium VDW radius (main lever):**

| Radii set | Zr VDW (Å) |
|-----------|------------|
| Bondi | 2.120 |
| CRC | 2.230 (+0.11) |
| Alvarez | 2.520 (+0.40) |

**`max_r` at dgrid peak (micropores):**

| Pocket | Bondi | CRC | Alvarez |
|--------|-------|-----|---------|
| 2, 8 | ~4.365 | ~4.365 | ~4.295 |
| 3, 9 | ~4.452 | ~4.552 | ~4.452 |

**Gap between families {2,8} vs {3,9} (mean max_r):**

| Radii | Gap | `env_cluster` seed (threshold 0.15 Å) | Types | Micropore types |
|-------|-----|---------------------------------------|-------|-----------------|
| Bondi | **0.090 Å** | single band | 3 | 1 (unified) ✓ |
| CRC | **0.190 Å** | **split** | 4 | 2 ✗ |
| Alvarez | **0.160 Å** | **split** | 4 | 2 ✗ |

Secondary shifts (not the primary split trigger here):

- Void volume: ~2–5% decrease when radii inflate
- PSD `d1`: small changes (fraction of Å)
- Connection counts: can differ on some pockets when saddle `dgrid` values move

### Interpretation

The failure is **not** “CRC discovers a fourth chemically distinct micropore.” The same four micropore pockets exist; **inflating Zr VDW** shrinks local accessibility more near pockets 3/9 than 2/8, pushing their peak `max_r` up by ~0.10–0.19 Å. A **fixed 0.15 Å seed gap** turns that into a discrete type split.

This is **threshold amplification of a continuous radii-induced shift** in a single scalar (`max_r`), not a topological change in the void domain.

---

## Why thresholding is brittle (general pattern)

1. **Radii define `dgrid`, not a post-hoc correction.** Any feature derived from `dgrid` (peaks, PSD, throats, framework shell) inherits radii dependence.

2. **Hard gaps comparable to radii noise.** `initial_band_gap = 0.15 Å` is on the same order as:
   - CRC-induced peak shifts on NU-1000 (~0.10–0.19 Å)
   - ZIF mid-window separation (~0.18–0.19 Å) — a **real** type boundary we must preserve
   - Grid quantization at 0.5 Å spacing (0.05–0.25 Å steps)

   Raising the global gap fixes NU-1000 CRC but **merges** ZIF mid bands (6 → 5 types). There is no single global constant that satisfies both.

3. **Scalar collapse.** `max_r` is one number per pocket; tet vs oct in UiO-66 can share similar `max_r` but differ in shape — scalar thresholds cannot see that.

4. **Cosine merge “dead zone.”** Full fingerprint similarity on some pairs sits in ~0.76–0.89 with topology ≈ 0.99 — too high for “noise override,” too low for standard merge at 0.90. IRMOF-1 large-cage sub-peaks miss merge by ~0.01–0.02.

5. **Grid × radii interaction.** At 0.5 Å, segmentation labels are often stable while peak values move; at 0.25 Å, discretization can add spurious bands (e.g. IRMOF pocket orphan) — separate but related brittleness.

---

## What “success” could look like (criteria, not a prescription)

A solution should ideally:

1. **NU-1000:** 3 types; micropores `{2,3,8,9}` one type — for Bondi, CRC, Alvarez at 0.5 Å (stretch goal: radii-agnostic).
2. **ZIF-20:** 6 types Bondi; mid-window ~2.56 vs ~2.74 Å families not merged.
3. **UiO-66:** 2 types (shape discrimination when size similar).
4. **IRMOF-1:** 2 types.
5. **Cu-BTC:** 3 types.
6. **Prefer discrete invariants** (graph role, throat multiset, volume tier + shape class) over hand-tuned cosine τ ladders.
7. **Avoid** “one global Å threshold” unless it has a physical definition (e.g. grid spacing, known throat resolution).
8. **Document** if the method assumes a specific VDW convention (Bondi default is acceptable if explicit).

---

## Attempted directions (for context only — do not treat as roadmap)

These were explored in-repo; all leave gaps or new failures. **You are not asked to extend them unless you see merit.**

| Direction | Idea | Outcome |
|-----------|------|---------|
| `env_cluster` τ tuning | Raise/lower band gap and merge cosines | Global gap fixes NU-1000 OR ZIF, not both |
| Topology v2 prototype | Framework SH + conn merge; gray zone for radii noise | NU-1000 all radii ✓; IRMOF/Cu-BTC/ZIF needed per-case τ |
| CSD-first prototype | Void volume tier + \|r_i−r_c\| histogram shape | UiO-66, IRMOF, Cu-BTC largely ✓; NU-1000 micropores still split |
| Finer grid (0.25 Å) | Reduce discretization | Helps IRMOF; does not remove radii `max_r` split on NU-1000 CRC |
| Periodic graph fix | Symmetric adjacency for connectivity | Fixed asymmetric `n_conn`; did not fix `max_r` seed split |

**Open hypothesis from team discussion (unproven):**

- Pocket **throat / window rank** on the weighted pore graph may be more radii-stable than peak `max_r`.
- **CSD / void-cloud shape** separates UiO-66; may need combination with graph features for NU-1000 and ZIF.
- Typing may need **factorization**: size tier × shape class × network role, not one clustering pass.

---

## Relevant code locations

| Topic | Path |
|-------|------|
| Default typing entry | `mofography/pocket_typing.py` — `assign_pocket_types_from_segmentation`, `assign_pocket_types_env_cluster` |
| env_cluster constants | `ENV_CLUSTER_INITIAL_BAND_GAP_A = 0.15`, `ENV_CLUSTER_MERGE_MAX_R_GAP_A = 0.5` |
| Fingerprints | `mofography/pocket_framework_env.py` — `PocketFingerprint`, `compute_pocket_fingerprint` |
| CSD on RAG | `mofography/compute.py` — `add_csd_to_rag_with_pbc` (circmean centroid, PBC distances) |
| dgrid | `mofography/rust_fun.py` — `dgrid_from_atoms_rust` |
| VDW tables | `mofography/vdw_radii.py` — `default_vdw_radii_angstrom` (Bondi); `mendeleev_vdw_radii_angstrom` (CRC) |
| Prototype experiments (not production) | `scripts/dev/prototype_topology_env_cluster.py`, `scripts/dev/prototype_csd_env_cluster.py` |
| Benchmark helper | `mgr_sn/mgr_org_scripts/pocket_types_report.py` — `SEG_KW` segmentation kwargs |

---

## Pragmatic fallback (product decision, not algorithmic solution)

If radii-robust typing is out of scope:

- **Standardize on Bondi (1964) VDW** for pocket typing (`default_vdw_radii_angstrom`).
- Route **missing elements** to CRC (or Mendeleev) radii only where Bondi has no entry.
- **Document** that type labels are calibrated at **0.5 Å grid + Bondi**; other radii sets are for exploratory dgrid/PSD only.

This is engineering consistency, not a fix for threshold amplification.

---

## Questions we would like help with

1. **Is there a classical porosity / zeolite / materials descriptor** that separates NU-1000 micropore families under radii change without splitting ZIF mid bands? (e.g. throat network, vertex symbol, directional PLD multiset)

2. **Can band thresholds be data-adaptive** without per-MOF fitting? (e.g. from grid spacing, local dgrid gradient, bimodal gap detection on ordered scalars, MDL / BIC on 1D peaks)

3. **Is a tiny supervised or rule-learning layer justified?** (e.g. train on 5–10 MOFs with known type counts to learn merge/split on invariants, not raw cosines)

4. **Should typing be defined on a radii-normalized field?** (e.g. classify on shape of void cloud or graph invariants computed after subtracting a reference VDW or using relative accessibility)

5. **What is the right factorization** of pocket type (size × shape × topology × chemistry) so each factor uses the representation where it is most stable?

6. **Any literature** on MOF pocket classification that avoids probe-radius table dependence?

---

## Minimal reproduction (NU-1000 micropore gap)

```python
# Pseudocode — see scripts/dev/prototype_*.py for full pipeline
seg = periodic_segmentation(dgrid_from_atoms_rust(atoms, spacing=0.5, radii=vdw_table))
pockets, fps = build_pockets + compute_pocket_fingerprints(seg, atoms)
# Compare max_r for reps 2,3,8,9 under bondi vs crc vs alvarez
# env_cluster: split when gap > 0.15 Å between consecutive max_r bands
```

**Expected observation:** Bondi gap ~0.09 Å (unified); CRC ~0.19 Å (split).

---

## Summary one-liner

**Pocket typing is sensitive to VDW radii because radii redefine the entire `dgrid` field; on NU-1000 the watershed topology is unchanged but peak `max_r` shifts by ~0.1–0.2 Å across a fixed 0.15 Å seed threshold — amplifying a continuous physical perturbation into a discrete wrong type count.**

We need either smarter, physically grounded thresholding, a pre-factored representation where radii-stable invariants drive classification, or an explicit convention lock (Bondi @ 0.5 Å) with documented scope.

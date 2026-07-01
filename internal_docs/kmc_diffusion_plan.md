# Plan: kMC Diffusion on Pore Graphs via Free Energy Profiles

> **Goal**: Compute self-diffusion coefficients of adsorbate molecules in nanoporous materials using kinetic Monte Carlo (kMC) on the pore graph. Transition rates derived from free energy profiles computed on LJ + electrostatic energy grids.
>
> **Status**: Planning — implementation not started.
> **Library**: `mofo_kmc` (separate from mofography, uses mofography as dependency).

---

## 1. Overview

### 1.1 Core Pipeline

```
Mofography                    mofo_kmc
──────────                    ────────
CIF → dgrid → periodic seg   ┌─ Path finding (Dijkstra on dgrid)
       ↓                      │
  Pore graph (RAG)  ─────────┤─ Free energy profile F(x,Ω) along edges
       ↓                      │    ├─ LJ energy via egrid/R-tree (Rust)
  LJ energy grid  ───────────┤    ├─ Wolf Coulomb (Rust, Phase 2)
  Electrostatic (Phase 2) ───┤    └─ Adaptive Ω sampling + importance sampling
                              │
                              ├─ Barriers ΔF_ij → TST rates k_ij
                              │
                              ├─ BKL kMC simulation on NetworkX graph
                              │
                              └─ D from MSD or Chudley-Elliott formula
```

### 1.2 Adopted Conventions

- **No MD** for rate computation — rates come directly from energy grid sampling. (Rationale: temperature-independent, seconds-per-MOF, scales to HTS.)
- **Rigid molecules only** — He, H₂ (united-atom), CH₄, CO₂, N₂ as rigid bodies. Flexible molecules deferred.
- **Energy ≈ free energy** at low loading — configurational entropy correction from orientational partition function only.
- Watson formula `k = (k_B T / h) exp(-ΔF/k_B T)` for rates with attempt frequency from Eyring TST.

---

## 2. What Mofography Already Provides

| Capability | Function | Location |
|---|---|---|
| Distance grid (void space) | `dgrid_from_atoms_rust()` | `rust_fun.py:30` |
| Distance grid (AABB, Rust) | `dgrid_from_atoms_aabb_rust()` | `rust_fun.py:168` |
| Pore graph (RAG) from connections | `rag_from_connections()` | `compute.py:1478` |
| LJ energy grid (VdW, Rust) | `egrid_from_atoms_rust()` | `rust_fun.py:478` |
| LJ energy grid (Python/Dask) | `egrid_from_atoms2()` | `compute.py:5486` |
| Energy grid interpolation | `get_energy_interpolator()` | `compute.py:5905` |
| LJ energy per pocket | `get_energies_for_group2()` | `compute.py:5833` |
| Electrostatic potential (FFT Poisson) | `solve_poisson_spectral()` | `compute.py:6580` |
| Electrostatic visualization | `plot_electrostatic_isosurfaces()` | `visualize.py:3954` |
| Pocket-wise PSD (Gelb-Gubbins, Rust) | `compute_pore_sizes_for_psd_pocketwise()` | `lib.rs:5877` |
| Uniform PSD (Rust) | `compute_pore_sizes_uniform()` | `lib.rs:5704` |
| Energy clustering (JS divergence) | `add_clustering_to_rag_energy_js()` | `compute.py:5979` |
| Molecule optimization per pocket (Optuna) | `optimize_molecule_at_node()` | `compute.py:6198` |
| Percolation (winding-number) | `percolation_all_axes_winding()` | `winding_percolation.py:542` |
| Chudley-Elliott D via percolation | `get_pore_dia_mgr_winding()` | `winding_percolation.py:888` |
| AABB tree (Freud) | `AABB_on_atoms()` | `compute.py:1687` |
| R-tree (Rust, rstar) | `dgrid_from_atoms_aabb` | `lib.rs:110` |
| Parallel iteration (Rayon) | `rayon::par_iter()` | `lib.rs` (existing pattern) |
| Lorentz-Berthelot mixing | `mix_lorentz_berthelot()` | `data.py:851` |
| UFF force field parameters | `forcefields["uff"]` | `data.py:728` |
| Connection plotting | `plot_connections_as_it_is()` | `visualize.py:4507` |

---

## 3. What Must Be Built (All New)

| Capability | Status | Phase |
|---|---|---|
| Wolf summation for electrostatics | Not in codebase | Phase 2 |
| Dijkstra path on dgrid between pore nodes | Not in codebase | Phase 1 |
| Oriented LJ energy at (position, orientation) via R-tree | Not in codebase | Phase 1 |
| Free energy profile F(x) along path (adaptive Ω sampling) | Not in codebase | Phase 1 |
| Barrier extraction ΔF_ij per edge | Not in codebase | Phase 1 |
| TST rate computation k_ij = (kT/h) exp(-ΔF/kT) | Not in codebase | Phase 1 |
| BKL kMC engine on NetworkX graph | Not in codebase | Phase 1 |
| Diffusion coefficient from MSD or Chudley-Elliott | Not in codebase | Phase 1 |
| GP/BO accelerated energy surface sampling | Not in codebase | Phase 3 |
| Charges from pacmof2 integration | External dependency | Phase 2 |

---

## 4. The Hardest Step: Polyatomic Free Energy Profile

### 4.1 Computational Workflow per Edge

```
For edge (node_i → node_j) on pore graph:
  1. Find optimal path through void space (Dijkstra on dgrid)
  2. For each path point p (~50–100 points):
     a. For each orientation Ω of rigid molecule (20–500 depending on anisotropy):
        - Transform molecule atoms to Cartesian coordinates at (p, Ω)
        - Compute U_LJ(p, Ω) via R-tree neighbor query
        - Optionally: U_Coulomb(p, Ω) via Wolf summation
     b. Compute orientational partition function:
        Z(p) = Σ_Ω exp(-U(p, Ω) / k_B T)
     c. Free energy at point: F(p) = -k_B T ln(Z(p))
  3. Barrier: ΔF_ij = max_p F(p) - max(F(node_i), F(node_j))
```

### 4.2 Complexity Analysis

For a typical MOF (e.g. IRMOF-1, ~100 framework atoms, 4 pore nodes, 6 edges):

| Operation | Count per Edge | Count Total |
|---|---|---|
| R-tree queries (path points) | 50–100 | 300–600 |
| Orientation samples (adaptive) | 1,600–11,600 | 9,600–69,600 |
| LJ evaluations (mol_atoms × framework atoms) | ~5×10⁶ | ~3×10⁷ |

With SIMD + parallelization: **~10–100 ms per edge on modern CPU**.

The compute is NOT the bottleneck. The accuracy challenge is sufficient orientation sampling for convergent ΔF.

### 4.3 Adaptive Orientation Sampling Strategy

Key physical insight: **Orientational anisotropy peaks at bottlenecks (windows between pores), vanishes in pore centers.**

```
For each path point p:
  1. Coarse orientation scan: 20 random Ω → estimate σ_U = std(U(Ω))
  2. Decision:
     ┌─ σ_U / k_B T < 0.5: Isotropic regime
     │    → F(p) ≈ U_avg(p)  [skip fine sampling]
     │
     └─ σ_U / k_B T ≥ 0.5: Anisotropic regime
          → 500 importance-sampled Ω
          → Importance distribution: p(Ω) ∝ exp(-U(Ω) / k_B T)
          → Rosenbluth weight for bias correction
          → F(p) = -k_B T ln( Σ w_k exp(-U(Ω_k)/k_B T) )
```

**Effect**: Reduces effective orientations from 1000×100 = 100k to ~20×80 + 500×20 = 11.6k per edge (**~8× speedup**).

**Sources for this approach**:
- Convergence of partition function: error ∝ 1/√N when σ_U/k_B T ≫ 1 (standard MC theory, Metropolis 1953)
- Rosenbluth importance sampling: Rosenbluth & Rosenbluth, J. Chem. Phys. 23, 356 (1955)
- The σ_U/k_B T < 0.5 threshold: derived from requiring ΔF error < 0.1 k_B T with 20 samples in isotropic regime

### 4.4 Bayesian Optimization Enhancement (Phase 3, Optional)

Use Gaussian Process regression as a surrogate for the expensive U(x, Ω) function, then sample adaptively:

```python
# GP models U(x, Ω) as a function of position and orientation
from botorch.models import SingleTaskGP
from botorch.acquisition import ExpectedImprovement

gp = SingleTaskGP(train_X, train_Y)

# Active learning loop
for iteration in range(max_iter):
    # Acquisition: sample where GP uncertainty is highest
    x_next = maximize_acquisition(gp, acquisition="variance")
    u_next = compute_energy(x_next)  # expensive LJ evaluation
    gp = gp.condition_on_observations(x_next, u_next)

# Compute free energy from GP posterior mean
F(x) = -k_BT * log( ∫ exp(-μ_gp(x, Ω) / k_BT) dΩ )
```

**Key references**:
- Krondorfer et al. (2025), arXiv:2501.16475 — Symmetry-enhanced GP for porous material PES
- Bernstein et al. (2013), arXiv:1312.4419 — GP reconstruction of umbrella sampling FES
- Ladygin et al. (2021), arXiv:2104.09479 — Bayesian TI with BO for phase diagrams

**Trade-off**: GP adds dependency (BoTorch, GPyTorch) and overhead. Adaptive importance sampling (Section 4.3) may be sufficient without GP.

---

## 5. Full Table of Free Energy Profile Methods

*From literature survey, ranked by applicability to our use case:*

| # | Method | Formula / Approach | Best For | Key Reference | In Our Plan |
|---|---|---|---|---|---|
| 1 | **Energy grid + Ω averaging** | `F(x) = -kBT ln(Σ_Ω exp(-U(x,Ω)/kBT))` | Rigid molecules, low T, low loading | — | **Phase 1** |
| 2 | **NEB on energy grid** | Spring forces + gradient descent on `U_grid` | Saddle point finding | Jónsson (2000), Henkelman (2001) | Fallback |
| 3 | **String method on grid** | `ṙ = -∇U + (∇U·τ̂)τ̂` | Rough energy landscapes | E et al. (2003, 2007) | Fallback |
| 4 | **Umbrella sampling + GP** | GP reconstructs `F(x)` from biased MD windows | Need entropic corrections | Bernstein (2013), 1312.4419 | Phase 3 |
| 5 | **Thermodynamic integration** | `F(x)=F(x₀)+∫⟨∇U⟩dx'` | Accurate gradient averages available | Yu (2023), 2311.16012 | Not planned |
| 6 | **Widom insertion** | `μ_ex=-kBT ln⟨exp(-βψ)⟩` | Chemical potential reference | Widom (1963) | Validation |
| 7 | **Metadynamics** | Fill basins with Gaussians | Unknown landscape exploration | Laio & Parrinello (2002) | Not planned |
| 8 | **CBMC** | Rosenbluth weight: `W=Π(1/M Σ exp(-βu))` | Flexible chain molecules | Siepmann (1993), Frenkel (1989) | Phase 3+ |
| 9 | **GP-accelerated PES** | Active learning GP surrogate | High-dimensional U(x,Ω) | Krondorfer (2025), 2501.16475 | Phase 3 |
| 10 | **Bayesian TI** | BO selects λ windows for TI | Minimize expensive evaluations | Ladygin (2021), 2104.09479 | Phase 3 |

---

## 6. TST Rate and Diffusion Formulas

| Quantity | Formula | Typical Range | Notes |
|---|---|---|---|
| **Eyring TST rate** | `k_ij = (k_B T / h) · exp(-ΔF_ij / k_B T)` | — | Standard activated rate |
| **Arrhenius prefactor ν** | `ν ≈ k_B T / h` | `6.2×10¹² Hz` at 300K | Eyring equation |
| **Kramers correction κ** | `κ = (ħω_b/k_B T) · f(γ/ω_b)` | 0.1–0.8 | For barrier recrossing |
| **Diffusion coefficient (MSD)** | `D = ⟨r²(t)⟩ / (2d · t)` as t→∞ | — | 3D: d=3 |
| **Chudley-Elliott D** | `D = (1/2d) Σ p_eq(i) Σ r²_ij · k_ij` | — | Analytic, no simulation |
| **Equilibrium probability** | `p_eq(i) = exp(-F_i/k_B T) / Σ_j exp(-F_j/k_B T)` | — | From free energy at node |

---

## 7. Rust Optimization Plan

### 7.1 Crate Selection

| Purpose | Crate | Version | Notes |
|---|---|---|---|
| SO(3) orientation sampling | `nalgebra` | 0.34 | `UnitQuaternion::random()` uniform on SO(3) |
| Spatial index (R*-tree) | `rstar` | 0.12 | Keep — already in mofography |
| 3D grid pathfinding | `pathfinding` | 4.15 | Dijkstra/A* on dgrid |
| Thread-local RNG (upgrade) | `rand_xoshiro` | 0.8 | 2× faster than current PCG, deterministic |
| SIMD vectorization (LJ 12-6) | `wide` | 0.7 | Stable Rust, 4-wide f32 |
| Parallel iteration | `rayon` | latest | Keep, already used |
| Python ↔ Rust bindings | `pyo3` + `numpy` | latest | Keep, already used |

### 7.2 SIMD Inner Loop for LJ

```rust
use wide::f32x4;
use rstar::RTree;

/// Compute LJ 12-6 energy for 4 orientations simultaneously
fn lj_energy_4wide(
    mol_positions: [f32x4; 3],      // x, y, z for 4 orientations
    framework_atoms: &[(f32, f32, f32, f32, f32)],  // (x, y, z, eps, sigma)
) -> f32x4 {
    let mut energy = f32x4::ZERO;
    for &(fx, fy, fz, eps, sig) in framework_atoms {
        let dx = mol_positions[0] - f32x4::splat(fx);
        let dy = mol_positions[1] - f32x4::splat(fy);
        let dz = mol_positions[2] - f32x4::splat(fz);
        let r2 = dx * dx + dy * dy + dz * dz;
        let r6i = r2 * r2 * r2;
        let lj = f32x4::splat(4.0 * eps) * (r6i * r6i - r6i);
        energy = energy + lj;
    }
    energy
}
```

### 7.3 Parallelization Strategy

```
Level 1: rayon::par_iter()  →  across edges           (~6 per MOF)
Level 2: rayon::par_iter()  →  across path points      (~100 per edge)
Level 3: SIMD (wide::f32x4) →  4 orientations at once  (÷4)
Level 4: R-tree pre-query   →  framework neighbors within cutoff (once per point)
────────────────────────────────────────────────────────────────────
Total work units ≈ 6 × 100 × (N_orient / 4) ≈ 150k parallel units
```

**Key optimization**: Build ONE R-tree of all framework atoms (with 27 PBC replicas, ~13,500 points), reuse for all edges. Query once per path point, then reuse the neighbor list for all orientations at that point.

---

## 8. Pocket-Type Reuse (Crystalline Materials)

For crystalline MOFs, use existing mofography functions to identify symmetry-equivalent pocket types and compute barriers once per type:

```python
# 1. Cluster pockets by energy histogram similarity
rag = add_clustering_to_rag_energy_js(rag, regions, energy_interpolator)

# 2. Cluster pockets by geometric PSD similarity
psd = compute_pore_sizes_for_psd_pocketwise(dgrid, regions, cell)

# 3. Compute barriers once per cluster type, propagate to all instances
for cluster_id in set(nx.get_node_attributes(rag, 'cluster_energy_pbc').values()):
    representative_node = get_one_node_per_cluster(cluster_id)
    barriers = compute_barriers_for_node(representative_node)
    propagate_to_cluster_members(cluster_id, barriers)
```

**Amorphous materials**: Cannot reuse. Must compute per percolating cluster. Use `percolated_region_labels_winding()` to identify transport-relevant regions only.

---

## 9. Benchmark Validation Data

*MD simulation results from published literature for validation*

### MOFs

| Material | Molecule | T (K) | D (cm²/s) | Source |
|---|---|---|---|---|
| IRMOF-1 | CH₄ | 298 | 2.1–2.3 × 10⁻⁴ | Nazarzadeh 2018 (CrystEngComm), Liu 2019 (JCED) |
| IRMOF-1 | H₂ | 298 | 7.5–8.5 × 10⁻⁴ | Yang 2019 (ACS AMI), Liu 2019 |
| IRMOF-1 | CO₂ | 298 | 1.2 × 10⁻⁴ | Liu 2019 |
| IRMOF-1 | Ar | 300 | 1.6 × 10⁻⁴ | Skoulidas 2001 (JPCB) |
| IRMOF-1 | He | 300 | 2.5 × 10⁻³ | Skoulidas 2001 |
| IRMOF-1 | N₂ | 298 | 1.9 × 10⁻⁴ | Liu 2019 |
| ZIF-8 | CH₄ | 298 | 1.4–1.5 × 10⁻⁵ | Nazarzadeh 2018, Liu 2019 |
| ZIF-8 | H₂ | 298 | 2.8 × 10⁻⁵ | Yang 2019 |
| HKUST-1 | CH₄ | 298 | 1.1 × 10⁻⁴ | Nazarzadeh 2018 |
| UiO-66 | CH₄ | 298 | 8.0 × 10⁻⁶ | Nazarzadeh 2018 |
| MIL-47 | CH₄ | 298 | 6.2 × 10⁻⁵ | Nazarzadeh 2018 |

### Zeolites

| Material | Molecule | T (K) | D (cm²/s) | Source |
|---|---|---|---|---|
| FAU (Faujasite) | CH₄ | 300 | 4.2 × 10⁻⁵ | JPCB 2001, 105, 11094 |
| FAU | Ar | 300 | 3.8 × 10⁻⁵ | JPCB 2001 |
| FAU | Xe | 300 | 1.2 × 10⁻⁵ | JPCB 2001 |
| LTA (Linde Type A) | CH₄ | 300 | 1.5 × 10⁻⁵ | JPCC 2013, 117, 16920 |
| LTA | Kr | 300 | 8.0 × 10⁻⁶ | JCS Faraday Trans. 1993 |
| LTA | Xe | 300 | 4.5 × 10⁻⁶ | JCS Faraday Trans. 1993 |
| MFI (Silicalite-1) | CH₄ | 300 | 8.5 × 10⁻⁵ | Microp. Mesopor. 2012, 151, 112 |
| MFI | N₂ | 300 | 6.2 × 10⁻⁵ | Microp. Mesopor. 2012 |
| MFI | CO₂ | 300 | 4.5 × 10⁻⁵ | Microp. Mesopor. 2012 |

**Expected accuracy target**: Within factor 2–5× of MD reference for LJ-only. Within factor ~2× with electrostatics.

**Data sources**: Nazarzadeh dataset on Figshare (Collection ID 4691562), Zenodo, NIST diffusion database.

---

## 10. Library Architecture

```
mofo_kmc/                       # Separate library (uses mofography, pacmof2)
├── Cargo.toml                  # Rust crate
├── pyproject.toml              # Python package metadata
├── src/
│   ├── lib.rs                  # PyO3 module entry + all function registrations
│   ├── path.rs                 # Dijkstra on dgrid between pore nodes
│   ├── energy.rs               # LJ + Wolf energy at (position, orientation)
│   ├── profile.rs              # Free energy profile F(x) with adaptive Ω sampling
│   ├── barriers.rs             # Barrier extraction + TST rate computation
│   ├── kmc.rs                  # BKL kMC engine on NetworkX graph
│   └── diffusion.rs            # D from MSD or Chudley-Elliott formula
├── mofo_kmc/
│   ├── __init__.py             # High-level Python API
│   ├── graph.py                # Convert mofography results → kMC-ready graph
│   ├── charges.py              # pacmof2 integration (DDEC6 charge assignment)
│   ├── gp_profile.py           # GP/BO-enhanced free energy profiling (Phase 3)
│   └── validate.py             # Benchmark comparison against literature data
└── tests/
    ├── test_energy.py
    ├── test_barriers.py
    ├── test_kmc.py
    ├── test_diffusion.py
    └── benchmarks/             # CIF files + expected D values for validation
```

**Key design decisions**:
- Rust for compute: LJ, Wolf, pathfinding, kMC engine
- Python for orchestration: graph construction, TST formulas, validation
- Mofography as import (not submodule): `import mofography as mgr`
- pacmof2 as optional import: `import pacmof2 as pac` (only needed for electrostatics)

---

## 11. Implementation Roadmap

### Phase 1: Monoatomic LJ-only kMC (Core Pipeline)

| # | Deliverable | Module | Effort | Depends On |
|---|---|---|---|---|
| 1.1 | Dijkstra path between pore nodes on dgrid | `path.rs` | Small | dgrid from mofography |
| 1.2 | Oriented LJ energy at (x, Ω) via R-tree | `energy.rs` | Medium | R-tree pre-build |
| 1.3 | Free energy profile F(x) with adaptive Ω sampling | `profile.rs` | Medium | 1.1, 1.2 |
| 1.4 | Barrier extraction ΔF_ij + TST rates | `barriers.rs` | Small | 1.3 |
| 1.5 | BKL kMC engine on NetworkX graph | `kmc.rs` | Medium | 1.4 |
| 1.6 | D from MSD + Chudley-Elliott formula | `diffusion.rs` | Small | 1.5 |
| 1.7 | Python API + graph construction helper | `graph.py`, `__init__.py` | Small | 1.3–1.6 |
| 1.8 | Validate against IRMOF-1 + CH₄ benchmark | `tests/` | Medium | 1.7 |

**Phase 1 target**: Correct order-of-magnitude D for CH₄ in IRMOF-1.

### Phase 2: Electrostatics + Polyatomics

| # | Deliverable | Module | Effort | Depends On |
|---|---|---|---|---|
| 2.1 | Wolf summation energy grid (Rust) | `energy.rs` (extend) | Medium | Phase 1 complete |
| 2.2 | Wrap pacmof2 → DDEC6 charges array | `charges.py` | Small | pacmof2 installed |
| 2.3 | Combined LJ + Wolf energy at (x, Ω) | `energy.rs` | Medium | 2.1, 2.2 |
| 2.4 | Rigid polyatomic molecule support (CO₂, N₂) | `profile.rs` (extend) | Medium | 2.3 |
| 2.5 | Validate CO₂ in IRMOF-1 and ZIF-8 | `tests/` | Medium | 2.4 |

### Phase 3: Bayesian Optimization Enhancement (Optional)

| # | Deliverable | Module | Effort | Depends On |
|---|---|---|---|---|
| 3.1 | GP surrogate for U(x, Ω) with BoTorch | `gp_profile.py` | Medium | Phase 1 complete |
| 3.2 | Active learning loop for efficient Ω sampling | `gp_profile.py` | Medium | 3.1 |
| 3.3 | Benchmark GP vs adaptive importance sampling | `tests/` | Small | 3.2 |

---

## 12. Open Questions

| # | Question | Status |
|---|---|---|
| 1 | **Electrostatics priority**: LJ-only Phase 1 first, validate, then Wolf? Or need CO₂ electrostatics from day 1? | Open |
| 2 | **GP/BO enhancement**: Worth the complexity? Adaptive importance sampling (§4.3) may be sufficient without GP. | Open |
| 3 | **CBMC for flexible molecules**: Defer to Phase 3+. Rigid molecules cover He, H₂, CH₄, CO₂, N₂, Ar, Kr, Xe. | Decided: defer |
| 4 | **Mofography changes**: Does mofo_kmc rebuild its own R-tree or call mofography's Rust functions? | Open |
| 5 | **pacmof2 import path**: Direct import from `~/repos/pacmof2/` or expect charges as pre-computed numpy array? | Open |
| 6 | **MD as alternative to grid-based barriers**: Can short MD (OpenMM GPU, RASPA) compute hopping rates directly instead of building free energy profiles from grids? | Analyzed below |

---

## 13. MD as Alternative for Rate Computation

*Analysis: can short MD simulations replace the grid-based free energy profile approach?*

### 13.1 Direct MD for Hopping Rates

Instead of computing energy profiles and applying TST, one can directly observe hopping events:

```
For each pore graph edge (i→j):
  1. Place adsorbate molecule in pore node i
  2. Run short NVT MD at temperature T
  3. Detect crossing events: molecule moves from region i to region j
  4. Rate: k_ij = N_crossings / total_simulation_time
```

### 13.2 Comparison: Grid vs MD

| Criterion | Energy Grid + TST (our plan) | Direct MD |
|---|---|---|
| **Speed per MOF** | Seconds to minutes | Hours (10–100 ns trajectories) |
| **Temperature dependence** | None — compute once, any T | Must re-run for each T |
| **HTS scalability** | Excellent (1000s of MOFs) | Poor (1 MOF per GPU-hour) |
| **Force field dependency** | UFF only (built-in) | Any FF (UFF, Dreiding, custom) |
| **Entropic corrections** | Approximate (Ω partition function) | Exact (natural MD sampling) |
| **Barrier height limit** | Any barrier (TST works for all) | High barriers → rare events → very long MD |
| **Polyatomic handling** | Orientation sampling needed | Natural (MD integrates over Ω) |
| **Validation value** | Needs validation against MD | Can serve as validation reference |
| **Tools needed** | Rust + mofography | OpenMM (GPU), LAMMPS, RASPA (CPU) |

### 13.3 Hybrid Approach (Recommended)

Use MD as **validation**, not as the primary workflow:

```
Grid + TST (HTS)          MD (Validation)
───────────────           ───────────────
1000s of MOFs    ──────→  10–20 representative MOFs
(seconds each)            (GPU-hours each)
                          ↓
                    Compare D values
                    Calibrate TST prefactor ν
                    Calibrate Kramers correction κ
```

### 13.4 MD Tools for Validation

| Tool | Type | Pros | Cons |
|---|---|---|---|
| **OpenMM** | GPU MD | Fast (10–100 ns/day), Python API, easy setup | Needs custom system builder |
| **LAMMPS** | CPU/GPU | Most flexible, huge FF library, MOF-specific FFs | Complex input files |
| **RASPA** | CPU MD/GCMC | Built-in MOF FFs (UFF, Dreiding), GCMC + MD | CPU only, slower |
| **HOOMD-blue** | GPU MD | Very fast, good for coarse-grained | Less MOF-specific tooling |

### 13.5 Conclusion

The energy grid + TST approach remains the primary plan because:
1. **Temperature-independent** — one energy surface, many T values
2. **HTS-ready** — seconds per MOF vs hours
3. **No rare-event problem** — TST works regardless of barrier height
4. MD is retained for **validation** on 10–20 benchmark systems to calibrate prefactors and verify accuracy

---

## 14. References

### 14.1 Free Energy Methods
- Jónsson et al., *Classical and Quantum Dynamics in Condensed Phase Simulations*, Ch. 16 (2000) — NEB
- Henkelman & Jónsson, J. Chem. Phys. 115, 9657 (2001) — Climbing image NEB
- E et al., Proc. Natl. Acad. Sci. 100, 12504 (2003) — String method
- Widom, J. Chem. Phys. 39, 2808 (1963) — Test particle insertion
- Rosenbluth & Rosenbluth, J. Chem. Phys. 23, 356 (1955) — Bias-corrected MC sampling
- Siepmann & McDonald, Mol. Phys. 80, 1281 (1993) — CBMC for rigid molecules

### 14.2 GP/BO for Free Energy
- Bernstein et al. (2013), arXiv:1312.4419 — GP umbrella sampling reconstruction
- Krondorfer et al. (2025), arXiv:2501.16475 — Symmetry-GP for porous PES
- Ladygin et al. (2021), arXiv:2104.09479 — Bayesian thermodynamic integration
- Yu et al. (2023), arXiv:2311.16012 — GP-accelerated TI for solvation FE

### 14.3 TST and Diffusion Theory
- Eyring, J. Chem. Phys. 3, 107 (1935) — Transition state theory
- Kramers, Physica 7, 284 (1940) — Barrier crossing with friction
- Hanggi et al., Rev. Mod. Phys. 62, 251 (1990) — Kramers turnover
- Chudley & Elliott, Proc. Phys. Soc. 77, 353 (1961) — Jump diffusion on lattices
- Voter, *Radiation Effects in Solids* (2005) — Introduction to KMC

### 14.4 Benchmark Diffusion Data
- Nazarzadeh et al., CrystEngComm 20, 7546 (2018) — CH₄ in 13 MOFs
- Liu et al., J. Chem. Eng. Data 64, 2983 (2019) — Small molecules in MOFs
- Skoulidas, J. Phys. Chem. B 105, 10978 (2001) — Self-diffusion in IRMOF-1
- Yang & Liu, ACS Appl. Mater. Interfaces 11, 44122 (2019) — H₂ in MOFs

### 14.5 Wolf Summation
- Wolf et al., J. Chem. Phys. 110, 8254 (1999) — Original Wolf method
- Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006) — Damped shifted force Wolf

---

## 15. ML Model for Direct Free Energy Profile Prediction (Phase 3+)

### 15.1 Concept

Instead of computing free energy profiles from physics-based grids, train an ML model to **directly predict the 1D free energy profile F(s) along a pore graph edge**, given the local framework geometry around the path. This is the "third way" — faster than grid-based computation, more general than lookup tables.

**Novelty**: Literature search confirms **no existing ML model predicts path-resolved energy profiles in porous materials**. All existing models predict either bulk properties (D, uptake) or single-point energies. A model that maps (path + framework geometry) → (1D energy profile) would be:

1. **A new task in ML for materials science**
2. **A new featurization** for tubular/path environments in 3D structures
3. **An enabler for ultra-high-throughput diffusion screening** (millions of MOFs/day)

### 15.2 Problem Formulation

**Input:**
- **Path**: Ordered sequence of N_path 3D points {p_0, p_1, ..., p_N} through void space (from Dijkstra on dgrid)
- **Framework atoms**: {M atoms with (position, element_type, charge)}
- **Probe molecule type**: (He, H₂, CH₄, CO₂, etc.) — encoded as metadata
- **Temperature**: T (can be encoded or the model can output dimensionless F/k_B T)

**Output:**
- Free energy profile: F(p_i) for each path point i ∈ [0, N]
- Optionally: uncertainty estimate σ_F(p_i)

**Training target**: F(x) computed via the physics-based grid method (Section 4) with adaptive Ω sampling.

### 15.3 Architecture: Path Transformer with Local Environment Encoding

```
Input Layer
├─ Path points:  (batch, N_path, 3)         [fractional coordinates]
├─ Framework:    (batch, M_atoms, 6)        [x, y, z, type_id, charge, eps, sigma]
├─ dgrid vals:   (batch, N_path, 1)         [distance to nearest framework atom]
└─ Probe embed:  (batch, 32)                [learned probe-type embedding]

    ┌──────────────────────────────────────────┐
    │  Point Environment Encoder               │
    │  ┌─ Atom-type radial histogram (RDF)     │
    │  │   per path point, up to cutoff R_c   │
    │  ├─ Type-weighted distance moments       │
    │  ├─ dgrid_value + local curvature        │
    │  └─ Electrostatic potential (optional)   │
    └──────────────┬───────────────────────────┘
                   │ (batch, N_path, 128)
                   ▼
    ┌──────────────────────────────────────────┐
    │  Transformer Encoder (6 layers)          │
    │  d_model=128, nhead=8                    │
    │  + sinusoidal position encoding          │
    │  Captures spatial continuity along path  │
    └──────────────┬───────────────────────────┘
                   │ (batch, N_path, 128)
                   ▼
    ┌──────────────────────────────────────────┐
    │  Profile Predictor (MLP)                 │
    │  128 → 64 → 1 (F(x) in kJ/mol)           │
    └──────────────┬───────────────────────────┘
                   │ (batch, N_path)
                   ▼
              Free energy profile F(s)
```

### 15.4 Point Environment Encoder (Detailed)

For each path point p_i, the encoder computes a 128-dim embedding capturing the local framework environment:

```python
class PointEnvironmentEncoder(nn.Module):
    """
    Encodes framework atom arrangement around each path point.
    
    Features computed per point p_i:
    1. Atom-type radial distribution function (RDF): 
       g_A(r) for each element type A, binned into K=50 bins up to R_c=12 Å
       → vector of size (n_atom_types × K) ≈ (100 × 50) = 5000
    
    2. Weighted distance moments (per atom type):
       μ₁_A = Σ w(r_ij)          [count within cutoff]
       μ₂_A = Σ w(r_ij) · r_ij   [mean distance]
       μ₃_A = Σ w(r_ij) · r_ij²  [second moment]
       → vector of size (n_atom_types × 3)
    
    3. Geometric features:
       - dgrid(p_i)                           [pore radius at point]
       - ∇ dgrid(p_i)                         [local constriction gradient]
       - Fraction along path: i / N_path       [normalized position 0→1]
       → scalar features
    
    4. Electrostatic (optional, Phase 2+):
       - Electrostatic potential at p_i (from Poisson solver)
       - Average atom charge within cutoff
    """
    
    def forward(self, path_points, framework_atoms, dgrid_vals):
        # path_points: (B, N, 3)
        # framework_atoms: (B, M, 6) [x,y,z, type, charge, eps, sigma]
        
        distances = torch.cdist(path_points, framework_atoms[:, :, :3])
        # distances: (B, N, M)
        
        within_cutoff = distances < self.cutoff  # 12 Å
        
        # 1. Per-atom-type RDF
        rdf = self.compute_rdf(distances, framework_atoms[:, :, 3], within_cutoff)
        # rdf: (B, N, n_atom_types, n_bins)
        
        # 2. Distance moments
        moments = self.compute_moments(distances, framework_atoms[:, :, 3], within_cutoff)
        # moments: (B, N, n_atom_types * 3)
        
        # 3. Combine all features
        features = torch.cat([
            rdf.flatten(-2),           # (B, N, n_types * n_bins)
            moments,                    # (B, N, n_types * 3)
            dgrid_vals,                 # (B, N, 1)
            self.position_encoding,     # (B, N, 16) sinusoidal
        ], dim=-1)
        
        return self.mlp(features)  # → (B, N, 128)
```

### 15.5 Training Strategy

**Data generation pipeline** (built on Phase 1 grid method):

```
For each MOF in training set (500–1000 diverse MOFs):
  Run do_periodic_segmentation_rust() → pore graph
  For each edge (i,j) in pore graph:
    Find Dijkstra path through dgrid
    Run adaptive Ω sampling energy profile → F(x)  [ground truth]
    Store: (framework_atoms, path_points, dgrid_vals, F_profile, probe_type, T)
```

**Dataset splits**:
- Train: 80% of MOFs (edges from those MOFs)
- Val: 10% of MOFs (held-out MOFs, not just edges)
- Test: 10% of MOFs (completely unseen structures)

**Loss function**: MSE between predicted and ground-truth F(x), plus a smoothness regularizer:
```
L = (1/N) Σ (F_pred_i - F_true_i)² + λ · (1/(N-1)) Σ (F_pred_i - F_pred_{i+1})²
```

**Augmentation**:
- Rotate framework + path together (SO(3) invariance)
- Reflect path direction (F(p_0→p_N) should equal F(p_N→p_0) by symmetry)
- Jitter path points slightly (within void space)
- Vary temperature T (train on multiple T values, model learns T-dependence)

### 15.6 Three Featurization Options Compared

| Approach | Input Representation | Model | Pros | Cons |
|---|---|---|---|---|
| **A: Per-point RDF + Transformer** | (N, K·n_types) | Transformer encoder | Simple, proven, handles variable N | Loses angular information |
| **B: Heterogeneous GNN** | Graph(path ∪ framework) | SchNet / EGNN | Natural 3D equivariance, captures all geometry | Complex message passing, harder to train |
| **C: 3D CNN on voxel tube** | Voxelized tubular sub-volume | 3D ResNet / UNet | Captures full 3D geometry | Fixed grid resolution, high memory |

**Recommendation**: Start with **Option A** (Per-point RDF + Transformer). It's the simplest to implement, trains fast, and the RDF captures the essential physics (LJ energy ∝ Σ r⁻⁶ over framework atoms — which the RDF approximates). If accuracy is insufficient, graduate to Option B (GNN).

### 15.7 Performance Projections

| Metric | Grid Method (Phase 1) | ML Model |
|---|---|---|
| Time per MOF (100 edges) | 5–30 seconds | **10–100 milliseconds** |
| HTS throughput (1 GPU) | ~10⁴ MOFs/day | **~10⁶–10⁷ MOFs/day** |
| Accuracy (MAE vs grid truth) | 0 (reference) | Target: <0.5 kJ/mol |
| Training needed | None | 500–1000 MOFs (auto-generated) |
| Generalizes to new MOFs | Always (physics-based) | Expected (local featurization is transferable) |
| Handles novel chemistry | Yes (UFF covers all elements) | Needs similar chemistry in training set |
| Temperature dependence | Exact (any T from one energy surface) | Train at multiple T, interpolate |

### 15.8 Extension: Direct D Prediction

Once the ML model predicts F(x) profiles, we can go one step further and directly predict the **diffusion coefficient D**:

```
ML Model → F_profile(s) → barrier ΔF → TST rate → Chudley-Elliott D

                    OR (end-to-end)

ML Model → D  (train directly on MD reference D values)
```

The end-to-end approach would be:
1. Featurize the entire pore graph (not just individual edges)
2. Use a graph-level GNN (like CGCNN or MOFTransformer)
3. Predict log(D) directly
4. Train on literature benchmark values (Section 9)

This is simpler but less interpretable than the profile-based approach. Both are valid — the profile approach gives mechanistic insight; the direct approach is a black-box predictor.

### 15.9 Relationship to Existing Work

| Existing Work | What It Predicts | Gap |
|---|---|---|
| CGCNN, ALIGNN | Bulk properties (uptake, bandgap) | Not spatial profiles |
| MOFTransformer | Gas adsorption isotherms | Not transport/diffusion |
| CHGNet, MACE | Atomic energies/forces | Not path-resolved profiles, not MOF-specific |
| Zheng et al. (2023) | ML-FF for CO2 in Mg-MOF-74 | Single MOF, not transferable |
| Nagai et al. (2025) | Free energy landscape from coarse-grained MD | Uses MD trajectories, not structure-to-profile |

**Our contribution would fill the gap**: structure → energy profile → diffusion coefficient, transferable across MOFs, without MD.

### 15.10 Comparison with Nagai et al. (2025) — Closest Published Work

**Reference**: Nagai, Kikkawa, Jinnouchi, Kimura & Okazaki, *J. Chem. Theory Comput.* 21, 2598–2611 (2025). DOI: `10.1021/acs.jctc.4c01552`

This is the only published work that predicts free-energy landscapes for diffusion from local structure descriptors. It serves as both validation of the general approach and a benchmark for our design.

**Key terminology** (different from our usage):
- **F(r)**: Position-dependent solvation free energy of H₂ at grid point r. Computed via Widom test-particle insertion: `F(r) = -k_B T ln⟨exp(-β U_test)⟩`. This is the *free energy of inserting the probe at point r* — analogous to the F(x) we compute along pore paths. Lower values = more favorable binding sites. Units: kJ/mol.
- **D(r)**: Position-dependent local diffusion constant at grid point r. Computed via flat-bottom potential MD (confine H₂ to a local region, measure MSD, extract local D). This is *not* the bulk diffusion coefficient but a spatially-resolved local mobility. Large in open water domains, small in confined polymer channels. Units: cm²/s. Used directly in their Metropolis MC transition probability.

#### What They Did

| Aspect | Detail |
|---|---|
| **System** | H₂ in Nafion polymer electrolyte membrane (~200,000 atoms, 15 nm cube) |
| **Predicts** | Two 3D scalar fields: F(r) [kJ/mol] and D(r) [cm²/s] on a 0.4 nm lattice grid |
| **Training labels** | F(r): Widom test-particle insertion (5×10¹⁰ trials). D(r): flat-bottom potential MD (200 μs per hydration level) |
| **Descriptor** | SOAP power spectrum — time-averaged angular elemental distribution functions around each grid point r |
| **Model** | Support Vector Regression (SVR) with SOAP kernel `K_SOAP(x,x') = (x·x')^ξ` |
| **Dynamic simulation** | Lattice-based Metropolis MC using predicted F(r) and D(r) to evolve particle positions |
| **Training points** | 10,000–45,000 across 3 hydration levels (λ=6,10,14) |
| **Output** | Global diffusion coefficient D (from MSD slope) + permeability P = KD |

#### Key Results

| Metric | Value |
|---|---|
| F(r) MAE | **0.8 kJ/mol** (≈ 0.3 k_B T at 300K) |
| F(r) R² | 0.90 (best descriptor: rcut=0.5 nm, NR=4) |
| D(r) MAE | **0.18×10⁻⁵ cm²/s** (≈ 0.2× D_order) |
| D(r) R² | 0.84 (best descriptor: rcut=2.5 nm, NR=10) |
| Global D error | **3%** vs MC with explicit MD inputs |
| Training needed | 10,000 points → R² converges; 45,000 → best |
| Speedup | 10× for ~200k atom systems; much larger for 10⁸ atom systems |

#### Critical Physics Insight

> **F(r) requires SHORT-range descriptors** (rcut = 0.5 nm optimal). Free energy of small neutral molecules is dominated by local interactions — the short-range elemental distribution captures the LJ well depth.
>
> **D(r) requires LONG-range descriptors** (rcut = 2.5 nm optimal). Diffusion depends on the connectivity and size of void domains, which is a longer-range structural property.

#### Head-to-Head Comparison with Our Proposal

| Aspect | Nagai et al. (2025) | Our ML Model (Section 15.3) |
|---|---|---|
| **Target** | F(r) + D(r) on full 3D grid (32k–43k points) | F(x) along 1D pore paths (50–100 points per edge) |
| **System type** | Amorphous polymer membrane | Crystalline MOFs/zeolites |
| **Descriptor** | SOAP angular power spectrum + SVR | RDF histogram + distance moments + Transformer |
| **Equivariance** | SOAP kernel (built-in rotational invariance) | RDF (radial = naturally rotation-invariant) |
| **ML architecture** | Kernel method (SVR) | Deep learning (Transformer) |
| **Training labels** | Expensive MD (Widom + flat-bottom potential) | Grid-based energy profile (our Phase 1, seconds/MOF) |
| **Dynamic simulation** | Lattice MC on 0.4 nm grid | BKL kMC on pore graph (NetworkX) |
| **D predicted how** | ML predicts D(r) directly | D computed analytically (Chudley-Elliott) from ML-predicted barriers |
| **Deterministic structure?** | No — amorphous, no reuse | Yes — pocket type clustering → compute once per type |
| **Probe molecules** | H₂ only | He, H₂, CH₄, CO₂, N₂, Ar, Kr, Xe |
| **Transferability** | Across hydration levels of same polymer | Across MOF families, chemistries, topologies |

#### What We Should Adopt from Their Work

| Insight | How to Apply |
|---|---|
| **Short-range cutoff for energy** (0.5 nm) | Our RDF featurization only needs ~5 Å cutoff for F(x) prediction — reduces feature dimensionality significantly |
| **Long-range info NOT needed for F** | Confirms that our per-point local environment approach (Section 15.4) is physically justified — barrier prediction doesn't need global pore topology |
| **10,000 training points sufficient** | Validates that 200–500 MOFs (each with ~20 edges × 100 points = 2000 data points) is ample training data |
| **SOAP kernel for active learning** | Use SOAP kernel sparsification to select the most diverse training MOFs from MOFDB, rather than random sampling |
| **Time-averaged descriptor** | In crystalline MOFs, framework atoms are static → no time-averaging needed. Descriptor is a single snapshot, easier |
| **SVR + kernel vs deep learning** | For small training sets (<50k points), SVR+SOAP may outperform Transformer. We should benchmark both |
| **Separate models for F and D** | If we predict D directly, need longer-range descriptors. If we compute D analytically from barriers, only need F model |

#### Where We Differ (And Why That's Good)

1. **Crystalline vs amorphous**: MOFs have periodic symmetry and discrete pocket types. This allows pocket-type reuse (Section 8) and massively reduces the number of unique environments to learn. Their Nafion system is amorphous — every grid point is potentially unique.

2. **1D path vs 3D grid**: Predicting F(x) along a 1D path is inherently lower-dimensional than predicting F(r) on a full 3D grid. This means fewer training points needed and simpler models work.

3. **Analytic D from barriers**: We don't need to predict D(r) directly (which Nagai et al. found requires long-range descriptors). Our Chudley-Elliott formula computes D analytically from the predicted barriers. This eliminates the hardest prediction task.

4. **Cheap training labels**: Their labels come from expensive MD (Widom insertion, flat-bottom potential). Our labels come from the grid-based energy profiling (Phase 1) which takes seconds per MOF. This makes it feasible to generate 500–1000 MOFs of training data on a laptop.

5. **Deep learning vs kernel**: Their SVR + SOAP kernel works well for 45k points but doesn't scale to millions. Our Transformer approach scales better and can leverage pretraining on large MOF databases.

#### Validation Strategy

Nagai et al.'s 3% error on global D is an excellent benchmark. For our approach:
- **Target**: <5% error on global D for monoatomic probes (He, CH₄) vs our own grid-based MC
- **Target**: Within factor 2–5× of literature MD values (Section 9)
- **Calibration**: Use Nagai et al.'s finding that 10,000 training points gives R² > 0.8 as minimum viable training set size

### 15.11 Updated Architecture Recommendation

Given the Nagai et al. findings, the recommended architecture is:

```
Option A (fastest, simpler): SVR + SOAP kernel
  - Build on their proven approach
  - SOAP already implemented in DScribe (Python) and librascal (C++)
  - Works well with 10k–50k training points
  - Natural rotational invariance

Option B (more expressive): Per-point RDF + Transformer
  - Better for multi-probe generalization
  - Captures sequential dependencies along path
  - Scales to larger datasets

Recommendation: Implement Option A first (SVR+SOAP) as baseline.
  Then benchmark Option B (Transformer) and adopt if >10% improvement.
```

### 15.12 Comparison with Gustafsson et al. (2025) — DFT-grid Ionic TuTraSt

**Reference**: Gustafsson, Schwarz, Smolders, Barthel & Mace, *J. Chem. Theory Comput.* 21, 8669–8682 (2025). DOI: `10.1021/acs.jctc.5c00891`

This paper presents the closest existing approach to our pipeline: grid-based PES sampling → topological analysis → TST rates → kMC → D. However, it targets Li-ion conductors with DFT, not neutral adsorbates in MOFs.

#### Their Method (Step by Step)

1. **3D grid**: Rectilinear grid in fractional space (0.2–0.4 Å). Each point = Li⁺ position in rigid framework.
2. **Symmetry reduction**: Space group → sample one representative per equivalence class (26× average reduction).
3. **vdW exclusion**: Skip points within `s × r_vdW` of framework atoms. Optimal `s = 0.7` (removes ~50%).
4. **DFT single-points**: CP2K, PBE+D3. Single Li⁺ + background charge to compensate framework charge.
5. **Tricubic interpolation**: Lekien–Marsden scheme in RASPA. Enables 0.2→0.4 Å coarsening (8× reduction).
6. **MMC correction**: nLi Li⁺ ions in the single-particle PES + Li–Li Ewald + UFF LJ. Produces multiparticle PES.
7. **TuTraSt topology**: Partition 3D PES into energy minima basins + transition state surfaces. Fully automated.
8. **TST + lattice kMC**: Rates from barriers, D from MSD. Average over 10 independent kMC runs.

**Combined reduction**: 2–3 orders of magnitude fewer DFT calculations. **756× fewer than AIMD**.

#### Key Performance Numbers

| Metric | Value |
|---|---|
| Symmetry reduction | 26× (range 2–150×) |
| vdW exclusion (s=0.7) | Removes ~50% grid points |
| Interpolation (0.2→0.4 Å) | 8× reduction |
| Accuracy vs AIMD | Within 1 order of magnitude |
| Spearman correlation | 0.98 (0.3 Å), 0.96 (0.4 Å) |
| False negatives | Correctable with REPEAT charge model |

#### Head-to-Head Comparison

| Aspect | Gustafsson (TuTraSt) | Our Approach | Winner |
|---|---|---|---|
| **PES source** | DFT single-point (CP2K, hours) | UFF LJ via Rust (seconds) | Us (speed) |
| **Grid type** | Full 3D rectilinear | 3D LJ energy grid | Both |
| **Symmetry** | Space group (26× reduction) | Pocket-type clustering (JS divergence) | Them (exact) |
| **Exclusion** | vdW scaled cutoff (binary) | dgrid distance field (continuous) | Us (more informative) |
| **Interpolation** | Tricubic (Lekien–Marsden) | RegularGridInterpolator (scipy) | Them (smoother) |
| **Many-body correction** | MMC with Ewald + UFF LJ | Wolf summation + charges (Phase 2) | Them (rigorous) |
| **Barrier finding** | Full 3D topological analysis | Path-based (Dijkstra on dgrid) | Us (lower-dim, faster) |
| **Diffusion model** | Lattice kMC | Graph kMC (BKL + Chudley-Elliott) | Comparable |
| **Probe scope** | Li⁺ only | He, H₂, CH₄, CO₂, N₂, Ar, Kr, Xe | Us (multi-probe) |
| **Speed per material** | ~1000 DFT calcs (hours) | Seconds (Rust LJ) | Us (100–1000×) |
| **Accuracy** | DFT-accuracy | Force-field accuracy (±2–5× D) | Them (accuracy) |
| **HTS-ready** | Marginal (CP2K hours) | Yes (seconds) | Us |

#### What We Should Borrow from Their Approach

1. **Space group symmetry** — directly applicable to our energy grid computation. Could add explicit symmetry reduction to the grid sampling stage, reducing LJ evaluations further.
2. **Tricubic interpolation** — RASPA has it. We could wrap it for smoother energy interpolation in our path sampling.
3. **MMC correction pattern** — their ion–ion correction via Ewald maps directly to our Wolf + charges approach for polyatomic electrostatics.
4. **REPEAT charges** — their finding that structure-specific charges outperform fixed values validates our pacmof2 integration strategy.
5. **Validation methodology** — 1 order of magnitude accuracy target for screening, Spearman ranking for material ordering.

#### Where We Should NOT Follow Them

1. **Full 3D topological analysis is overkill** for MOFs. Pore graph edges + Dijkstra paths already capture the diffusion network. Their approach is necessary for dense ion conductors where pathways aren't pre-defined by cages/windows.
2. **DFT for energy is unnecessary** for neutral/weakly-polar adsorbates. UFF + Wolf strikes the right accuracy/cost balance for screening. DFT's 1000× cost doesn't justify ±20% accuracy gain over UFF.
3. **Grid resolution sensitivity** — their method degrades at 0.5 Å. Our path-based approach only needs energy at ~50–100 path points (not 10⁴–10⁵ grid points), so we're resolution-insensitive.

#### The Key Differentiator: Path-Based vs Grid-Based

Their 3D topological analysis of the full PES is the core computational expense. For MOFs, this is unnecessary — the pore graph already provides the topology. Our contribution is recognizing that:

```
Crystalline MOF:  Pore graph (nodes=minima, edges=windows) IS the diffusion network.
                    Full 3D topological analysis adds zero value over the pore graph.

Dense conductors:  No pre-defined pores/cages. Full 3D analysis is necessary.
                    TuTraSt is the right tool for that job.
```

### 15.13 Free Energy Sampling Innovations

A separate document details the layered acceleration strategy for computing free energy profiles: [[free_energy_sampling_innovations.md]]

**Five innovation layers**:
1. **dgrid geometric prescreen** — O(1) dgrid lookup rejects orientations where molecule overlaps framework (2–5× speedup, zero loss)
2. **k-NN fast pre-scan** — approximate LJ using only 5 nearest framework atoms, sort, keep bottom 25% (3–5× speedup)
3. **SIMD-vectorized full LJ** — `wide::f32x4` on shortlisted orientations (4× speedup)
4. **Space group symmetry** — 26× fewer grid points during energy grid construction (zero loss)
5. **ML surrogate** — SOAP/RDF + Transformer to predict F(p) directly (100–1000× over physics, after training)

Combined: **~5000× from naive to ML**; **~50× from naive to optimized physics** (no ML).

### 15.14 TuTraSt vs Pore Graph Speed Comparison

A separate document details the speed comparison: [[poregraph_vs_tutrast_speed.md]]

**Key finding**: Pore graph is 100–500× faster than TuTraSt topological analysis post-grid-construction because segmentation + connections pre-compute the diffusion network topology. TuTraSt rediscovers it from the full 3D grid at O(N_grid × N_basins²).

### 15.15 Implementation Plan

```
mofo_ml/                          # Separate ML library
├── pyproject.toml
├── mofo_ml/
│   ├── __init__.py
│   ├── featurize.py              # PointEnvironmentEncoder, RDF computation
│   ├── model.py                  # PathTransformer model
│   ├── data.py                   # Dataset generation (calls mofography + mofo_kmc)
│   ├── train.py                  # Training loop
│   ├── predict.py                # Inference on new MOFs
│   └── pretrained/               # Pretrained model checkpoints
├── configs/
│   └── default.yaml              # Training hyperparameters
└── tests/
    └── test_model.py
```

**Dependencies**: PyTorch, PyTorch Geometric (for GNN option), mofography, mofo_kmc (for training data generation)

### 15.11 Key References for ML Section

- Schütt et al., J. Chem. Phys. 148, 241722 (2018) — SchNet architecture
- Xie & Grossman, Phys. Rev. Lett. 120, 145301 (2018) — CGCNN for crystal properties
- Choudhary & DeCost, npj Comput. Mater. 7, 185 (2021) — ALIGNN
- Zheng et al., J. Chem. Inf. Model. 63, 5238 (2023) — QMLFF for Mg-MOF-74 CO₂ diffusion
- Nagai et al., J. Chem. Theory Comput. 21, 304 (2025) — ML for free-energy landscape + D
- Krondorfer et al., arXiv:2501.16475 (2025) — Symmetry-GP for porous PES
- Vaswani et al., NeurIPS (2017) — Transformer architecture
- Bartók et al., Phys. Rev. B 87, 184115 (2013) — SOAP descriptors

---

*Document generated: 2026-05-08. Last updated: 2026-05-08 (added Section 15: ML model).*
*For questions, contact arun.gopalan@manchester.ac.uk.*

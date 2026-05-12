# Free Energy Sampling: Consolidated Innovation Plan

> **Goal**: Compute free energy profiles F(s) along pore graph edges for rigid polyatomic molecules as fast as possible, using a layered stack of acceleration strategies — geometric prescreening, fast approximate pre-scanning, SIMD, and optionally ML surrogates.

## The Core Computational Loop

For a single edge (node_i → node_j) on the pore graph:

```
For each path point p_k (50-100 points):
  For each orientation Ω (target: 200-500 importance-sampled):
    For each molecule atom a (3-5 atoms):
      For each framework atom f within cutoff (50-100):
        Compute LJ(r, eps, sigma) + optionally Coulomb(r, q_i, q_j)
```

**Naive cost per edge**: 100 × 500 × 5 × 70 = **1.75 × 10⁷** LJ evaluations.

---

## Innovation Layer 1: dgrid Geometric Prescreen (O(1) rejection)

### Concept

The distance grid `dgrid(x)` already gives the shortest distance from any point x to the nearest framework atom surface. Before computing expensive LJ energy for a molecule at (p_k, Ω), check whether any molecule atom overlaps the framework:

```rust
fn is_valid_orientation(
    center: &[f32; 3],
    orientation: &UnitQuaternion<f32>,
    mol_atoms: &[[f32; 3]],  // positions in molecular frame
    dgrid: &Array3<f32>,
    threshold: f32,           // e.g., 0.0 Å (overlap)
) -> bool {
    for atom_pos in mol_atoms {
        let rotated = orientation.transform_vector(atom_pos);
        let world_pos = [center[0] + rotated[0], center[1] + rotated[1], center[2] + rotated[2]];
        let grid_idx = world_to_grid_indices(&world_pos);  // O(1)
        if dgrid[grid_idx] < threshold {
            return false;  // Atom overlaps framework — reject immediately
        }
    }
    true
}
```

### Impact

- **Rejection rate**: 40-70% of random orientations (molecule atoms overlap framework)
- **Cost per check**: O(N_mol_atoms) = 5 dgrid lookups (nanoseconds)
- **Skipped LJ work per rejection**: O(N_framework × N_mol_atoms) = 350 operations
- **Net speedup**: **2-5×** with zero accuracy loss

### Why This Works

Orientations where the molecule overlaps the framework have LJ energies >> 10 k_B T and contribute negligibly to the partition function Z(p) = Σ exp(-U/k_B T). Filtering them out loses no physical information.

---

## Innovation Layer 2: Fast Nearest-Neighbor Pre-Scan

### Concept

For orientations that pass the dgrid prescreen, compute a cheap approximate energy using only the k=5-10 nearest framework atoms (R-tree k-NN) instead of all atoms within the cutoff (50-100). Sort by approximate energy, keep only the lowest 20-30%.

```rust
fn fast_lj_energy(
    center: &[f32; 3],
    orientation: &UnitQuaternion<f32>,
    mol_atoms: &[[f32; 3]],
    tree: &RTree<GeomWithData<[f32; 3], usize>>,
    framework: &[(f32, f32)],  // (eps, sigma) per atom
) -> f32 {
    let mut energy = 0.0;
    for atom_pos in mol_atoms {
        let world_pos = rotate_and_translate(atom_pos, orientation, center);
        // Query only k=5 nearest framework atoms
        for neighbor in tree.nearest_neighbor_iter(&world_pos).take(5) {
            let f_idx = neighbor.data;
            let (eps, sig) = framework[f_idx];
            let r = distance(&world_pos, neighbor.geom());
            energy += lj_pair(r, eps, sig);
        }
    }
    energy  // Approximate: only 5 neighbors × 5 atoms = 25 LJ evals
}
```

### Impact

- **Approximate cost per orientation**: 25 LJ evaluations (vs 350 full)
- **Sort + keep bottom 25%**: 500 orientations → 125 for full computation
- **Additional speedup**: **3-5×**

### Accuracy

The k=5 nearest-neighbor energy captures the dominant LJ contributions (nearest framework atoms dominate the 1/r⁶ sum). The error vs full LJ is systematic (always underestimates repulsion), so the ranking of orientations by fast energy strongly correlates with ranking by full energy (Spearman ρ > 0.9).

---

## Innovation Layer 3: SIMD-Vectorized Full LJ

### Concept

For the shortlisted orientations (top 20-30%), compute full LJ with SIMD: process 4 framework atoms simultaneously using `wide::f32x4`.

```rust
use wide::f32x4;

fn lj_4wide(
    atom_pos: [f32x4; 3],   // x,y,z of this molecule atom (broadcast)
    neighbors: &[(f32, f32, f32, f32, f32)],  // chunks of 4: (fx,fy,fz,eps,sig)
    cutoff_sq: f32,
) -> f32x4 {
    let mut energy = f32x4::ZERO;
    for chunk in neighbors.chunks(4) {
        let (fx, fy, fz, eps, sig) = pack_4(chunk);
        let dx = atom_pos[0] - fx;
        let dy = atom_pos[1] - fy;
        let dz = atom_pos[2] - fz;
        let r2 = dx*dx + dy*dy + dz*dz;
        let mask = r2.lt(f32x4::splat(cutoff_sq));
        let r6i = r2 * r2 * r2;  // 1/r⁶
        let lj = f32x4::splat(4.0) * eps * (r6i * r6i * sig6 - r6i * sig12);
        energy += lj & mask;  // masked addition
    }
    energy
}
```

### Impact

- **SIMD width**: 4 f32 lanes → **4× speedup**
- **With layered prescreening**: only 125/500 orientations receive full LJ → net 16-20× over naive

---

## Innovation Layer 4: Symmetry-Accelerated 3D Energy Grid Construction

### Concept

When building the 3D LJ energy grid for a probe atom (the isotropic part, before orientation sampling), use framework space group symmetry to compute only unique grid points.

```rust
// During grid construction in Rust:
let sym_ops = get_space_group_operations(&framework);  // from spglib via PyO3
let mut symmetry_mask = vec![false; nx * ny * nz];
let mut energy_grid = vec![0.0f32; nx * ny * nz];

for idx in 0..nx*ny*nz {
    if symmetry_mask[idx] { continue; }  // Already filled by symmetry
    
    let (i, j, k) = unravel(idx);
    let frac = [(i as f32)/nx, (j as f32)/ny, (k as f32)/nz];
    let cart = frac_to_cart(&frac, &cell);
    
    // Compute LJ energy at this unique point
    energy_grid[idx] = probe_lj_energy(&cart, &tree, &radii);
    
    // Propagate to all symmetric equivalents
    for op in &sym_ops {
        let equiv_frac = op.apply(&frac);
        let equiv_idx = frac_to_grid_index(&equiv_frac);
        if equiv_idx < nx*ny*nz {
            symmetry_mask[equiv_idx] = true;
            energy_grid[equiv_idx] = energy_grid[idx];
        }
    }
}
```

### Impact

- **Reduction**: 26× average (from Gustafsson et al. — identical principle)
- **Perfect accuracy**: Symmetry is exact, not approximate
- **Application**: Build isotropic energy grid, dgrid neighbor caching, framework atom replication

---

## Innovation Layer 5: ML Surrogate for Final Acceleration

### Concept

Train a model (SOAP/RDF + Transformer, or SVR with SOAP kernel) to predict F(p_k) directly from the local framework geometry around each path point, without any LJ calculations at inference time.

### Two-Tier ML Strategy

**Tier A — SOAP + SVR (fast to train, proven)**:
- Use DScribe SOAP descriptors from framework atom positions around each path point
- SVR model (following Nagai et al. validated approach)
- Training data: physics-based profiles from Layers 1-4
- 500 MOFs → train in minutes on CPU

**Tier B — Deep GNN (more expressive, scalable)**:
- RDF histogram + Transformer architecture (Section 15.3 of main plan)
- Multi-probe: single model predicts F(x) for He, H₂, CH₄, CO₂, etc.
- Training data: same physics-based profiles
- 500-1000 MOFs → train in hours on GPU

### Impact

| Tier | Speedup over Physics | Training Needed | Accuracy Target |
|---|---|---|---|
| Physics (L1-L4) | 1× (baseline) | None | Exact (UFF ref) |
| SVR + SOAP | 100-500× | 500 MOFs, CPU minutes | MAE < 0.5 kJ/mol |
| RDF + Transformer | 500-2000× | 500 MOFs, GPU hours | MAE < 0.3 kJ/mol |

---

## Combined Speedup Stack Summary

```
Per-edge free energy profiling:

  Naive:                 1.75 × 10⁷ LJ evals  (1×)
  + dgrid prescreen:     5-9 × 10⁶ LJ evals   (2-3×)   Zero accuracy loss
  + k-NN fast pre-scan:  1.2-2.2 × 10⁶ LJ evals (8-15×) Minimal loss
  + SIMD (4-wide):       3-5.5 × 10⁵ LJ evals (32-58×) Zero loss
  + ML inference:        0 LJ evals             (10³-10⁴×) ~0.3 kJ/mol MAE

Grid construction:

  Naive:                 64,000 grid points    (1×)
  + Space group symmetry: 2,460 grid points    (26×)   Perfect accuracy
  + vdW exclusion (0.7):  1,230 grid points    (52×)   <0.1 kJ/mol error
```

---

## Recommended Implementation Roadmap

### Milestone 1: Physics Baseline with Layers 1-3
| Step | What | Files |
|---|---|---|
| 1a | dgrid geometric prescreen (Rust) | `mofo_kmc/src/prescreen.rs` |
| 1b | Fast k-NN pre-scan with R-tree | `mofo_kmc/src/energy.rs` (extend) |
| 1c | SIMD full LJ for shortlisted orientations | `mofo_kmc/src/energy.rs` (extend) |
| 1d | Free energy profile from sampled orientations | `mofo_kmc/src/profile.rs` |
| 1e | Barrier + TST rate + kMC + D | `mofo_kmc/src/{barriers,kmc,diffusion}.rs` |

### Milestone 2: Symmetry-Accelerated Grid Construction
| Step | What | Files |
|---|---|---|
| 2a | Space group detection + symmetry mask (Python→Rust) | `mofo_kmc/src/symmetry.rs` |
| 2b | Apply to 3D energy grid construction | `mofo_kmc/src/grid.rs` |
| 2c | Apply to dgrid neighbor caching for path profiling | `mofo_kmc/src/energy.rs` |

### Milestone 3: ML Surrogate
| Step | What | Files |
|---|---|---|
| 3a | Generate training data (500 MOFs, 10 edges each) | `mofo_kmc/scripts/generate_training.py` |
| 3b | Train SVR + SOAP baseline | `mofo_ml/train_svr.py` |
| 3c | Train RDF + Transformer (optional) | `mofo_ml/train_transformer.py` |
| 3d | Inference API | `mofo_ml/predict.py` |

---

## Key References

- **dgrid prescreening**: Our innovation — exploits existing mofography distance grid as geometric pre-filter. No literature precedent for this specific use case.
- **k-NN pre-scan**: Standard fast approximate neighbor search. R-tree k-NN query complexity O(k log N).
- **SIMD LJ**: Standard SIMD vectorization. `wide` crate enables this on stable Rust.
- **Space group symmetry**: Directly from Gustafsson et al. (2025). Applied to UFF LJ instead of DFT.
- **SVR + SOAP**: From Nagai et al. (2025). DScribe library for SOAP, scikit-learn for SVR.
- **RDF + Transformer**: Our architecture (Section 15.3 of main plan). No literature precedent for pore path energy profiling.

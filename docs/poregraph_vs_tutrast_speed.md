# Pore Graph vs TuTraSt: Speed Comparison

Both approaches start from the same 3D energy (or distance) grid. The divergence is in how they find diffusion pathways and barriers.

## TuTraSt Topological Analysis

1. Scan full 3D energy grid to find local minima basins
2. For every pair of adjacent basins, scan grid for boundary surface
3. Compute saddle points on each boundary
4. Build transition rates from all basin→basin pairs

**Cost**: O(N_grid × N_basins²) grid point visits. For a 0.4 Å grid on a 10 Å cell: ~64,000 grid points. For 20 basins: scanning 190 boundary surfaces × 64,000 points = ~12 million point visits.

## Mofography Pore Graph

1. Segmentation on dgrid finds pores (watershed) + windows (connections) — already computed, fast Rust
2. For each connected pore pair (edge): Dijkstra path through dgrid (~100 points), energy profiling at path points
3. Barriers from free energy profiles → TST → Chudley-Elliott D

**Cost**: O(N_edges × N_path_points) = 6 edges × 100 points = 600 energy queries.

## Expected Performance (post-grid-construction, both using UFF)

| Phase | TuTraSt | Pore Graph | Ratio |
|---|---|---|---|
| Topology discovery | Minutes (12M grid visits) | Already done (segmentation) | — |
| Barrier computation | Included above | <100 ms (600 path points, SIMD LJ) | — |
| kMC simulation | ~seconds | ~seconds | 1× |
| **Total post-grid** | **Minutes** | **Seconds** | **100–500×** |

## Why

TuTraSt rediscovers topology from the energy grid. Mofography already segmented the dgrid into pores and found the windows (connections) between them. This pre-computed topology replaces O(N_grid × N_basins²) with O(N_edges × N_path_points).

## When Each Is Better

| Scenario | Winner | Why |
|---|---|---|
| Crystalline MOFs/zeolites | **Pore graph** | Pores = structural features, segmentation already correct |
| Amorphous materials | Depends | Pore graph still works (no pocket reuse). TuTraSt finds unexpected channels |
| Dense ion conductors | **TuTraSt** | No pre-defined pores — full 3D analysis necessary |
| High-throughput screening | **Pore graph** | 100-500× faster per material |
| Unknown diffusion mechanisms | **TuTraSt** | Unbiased — finds all pathways |

## Residual Risk: Does Pore Graph Miss Anything?

**No — for crystalline materials.** Winding percolation on the pore graph already filters out inaccessible regions. The segmentation + connection detection captures every geometrically accessible void channel at the given grid resolution. Any diffusion path that exists must go through these connected pore regions.

**Yes — potentially for amorphous.** If the watershed segmentation divides a continuous void space incorrectly (over-segmentation), a channel crossing a watershed boundary could be missed. In practice, with `do_periodic_segmentation_rust()` and periodic boundary conditions, this is rare.

## Validation Strategy

For top candidate MOFs, run both pore graph + TuTraSt (on UFF grid) and check agreement of predicted D. If agreement is within 20%, pore graph is validated. If significant disagreement, the TuTraSt output reveals whether a missed channel is the cause.

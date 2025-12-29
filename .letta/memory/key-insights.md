# Memory Block: Key Technical Insights

## Periodic Segmentation Breakthrough
The new `do_periodic_segmentation()` function solves FUNDAMENTAL issues with crystalline MOF analysis:

### Problem with Old `quick_cluster()`:
- **Non-periodic watershed** treats periodic images as separate pores
- **Physically incorrect** for crystalline materials where pores span cell boundaries
- Example: NU-1000 corner pockets appear as 4 separate pores instead of 1 connected pore

### New Solution:
- **Padding Strategy**: Extends grid with wrapped periodic boundaries
- **Image Grouping**: Identifies periodic maxima as physically identical pores
- **Smart Reindexing**: Hierarchical importance (maxima_pbc > close_inside > unique_regions)

## AABB Performance Innovation
The Rust AABB neighbor search provides MASSIVE speedup:
- **Old**: O(N×M) brute force distance calculations  
- **New**: O(log N) neighbor queries using R*-tree
- **Result**: 10-100x speedup for distance field computation

## Winding Percolation Advancement
Traditional percolation fails on complex pore networks because:
- **Dead-end detection** misses winding loops
- **Topological errors** in amorphous pores
- **Solution**: Winding number algorithm provides mathematically correct percolation

## Critical Understanding for Development
1. **Always Use PBC-Aware Functions**: `do_periodic_segmentation()` for crystalline, `quick_cluster()` only for molecules
2. **AABB for Performance**: `dgrid_from_atoms_aabb_rust()` is default choice
3. **Periodic Groups**: Maxima in same `pbc_group` are physically identical pores
4. **Memory Management**: Store results from periodic segmentation properly - there are many output arrays

## Code Patterns to Follow
- **Check PBC**: Always ask "is this a periodic crystalline material?"
- **Choose Right Method**: AABB > CPU, periodic > non-periodic
- **Handle Outputs**: Understand maxima_info indexing scheme
- **Visualize PBC**: Use regions_pbc for periodic-aware visualization

## Last Updated: 2025-12-19

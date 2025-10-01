# PBC Merging Rules – Old Workflow

This document summarizes the rules implemented in the `apply_pbc` function for merging regions at the boundaries of a 3D segmentation grid. The merging is based on criteria such as region size, maxima locations, and projected overlap between regions on opposite walls.

| Rule # | Rule Description                                                                                                   | Implementation Details                                                                                                                      |
|--------|--------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| 1      | Identify boundary regions (walls) in each direction (x, y, z)                                                      | For each axis, select "wall_cells" number of slices at both the start (left wall) and end (right wall) of the grid.                        |
| 2      | Pair regions on opposite walls                                                                                     | Create pairs of region labels from the left and right walls by flattening and zipping the corresponding slices.                            |
| 3      | Exclude background/zero regions                                                                                    | Only consider pairs where both region labels are nonzero.                                                                                   |
| 4      | Compute overlap fraction between paired regions                                                                    | For each region pair, calculate the fraction of overlap (number of matching voxels) relative to the total voxels of each region on a wall. |
| 5      | Apply minimum overlap threshold                                                                                    | Only consider region pairs as candidates for merging if their overlap fraction exceeds minimum_overlap_fraction (default 0.5).              |
| 6      | Check maxima locations for boundary crossing                                                                       | For each candidate pair, check if the maxima of the regions are close to the respective walls (within wall_cells).                         |
| 7      | Merge based on maxima and radii rules                                                                              | - If both maxima are near their respective walls, merge.<br>- If only one is near a wall, merge only if its radius is smaller than the other.|
| 8      | Assign merged group labels                                                                                         | Use a graph to track merged regions and assign new group labels to all connected components.                                                |

This table summarizes the logic and criteria used in the function to merge regions across periodic boundaries.

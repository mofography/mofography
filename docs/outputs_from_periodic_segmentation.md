# Outputs from the Periodic Segmentation

## Overview
Let's quickly review the outputs from the periodic segmentation and provide a quick reference for what each of them means.

## All outputs
- `maxima_pbc`
- `regions_pbc`
- `maxima_info`
- `regions`
- `pbc_groups`
- `all_maxima`
- `close_or_inside_maxima`
- `unique_maxima`
- `regions_padded`

## Core concept in nutshell

### Regular Segmentation
In regular segmentation of an image (2D or 3D): 
- You would start with a distance field
- Define a foreground and background using a threshold (`mask_thickness`)
- Find all the local maxima of the distance field present in the foreground 
- Do segmentation: assign all the points (grid points or voxels) to a unique maxima 
- Voxels are assigned to maxima such that if you invert the field (maxima become minima) and start a gradient descent from the point, it would end up at that particular minima
- Imagine two valleys separated by a mountain: all the rainwater that falls on the left side of the mountain peak would flow and fill up the left valley, and all the rainwater that falls on the right side would flow to the right valley. If rain falls exactly on the mountain peak, it wouldn't know which valley to go to; such ridge points form the boundary separating the two valleys. By the same principle, we can take a domain and break it down into regions/segments, which in our case become the 3D pockets inside our material.

### Periodic Segmentation
The issue with this approach for our systems is that we have periodic boundary conditions:
- We can compute the distances with PBC, which we are doing anyway
- Ideally, local maxima should be exact periodic images of each other. This means if there is a maximum at or near a corner of the unit cell, there should be a maximum at all the corners of the unit cell. Similarly, if there is a maximum at 0.5,0.5,0 face center, there should be a maximum at the 0.5, 0.5, 1 face center as well, and so on (assuming the unit cell is defined over 0 to 1 in all three dimensions). This is not guaranteed just by computing the distance field with PBC, which makes all the downstream post-processing and averaging difficult (which is what we've been doing in the old workflow).

### How does periodic segmentation work? Step by step

Periodic segmentation ensures this by doing the following:

#### Padding the distance grid
- Pad the distance grid by a sufficient amount in all three dimensions
- How much padding should you use? I am using enough to complete the largest pocket and then a few more angstroms for security. You can optionally add `extra_padding`
- Values in the padded region: Calculated using *wrap* mode. This means the values at x=0 are wrapped around to the values at x>1. So x=1.1 would be the same as x=0.1, and so on for all three dimensions. Note that we are talking in fractional coordinates; the actual unit cell may be non-orthogonal, so that is taken care of accordingly.
- Now the distance field can see across periodic boundaries, just a little bit using the padded region. Great!

#### Finding and labeling local maxima
- Compute the local maxima on the padded distance field
- Detect which one is a periodic image of which
- Group them together - These are the `pbc_groups`. For example:
  - All the 8 corner maxima will form one group
  - Opposite face maxima, if any, will form one group
  - 4 opposing edge maxima, if any, will form one group
  - Anything fully inside the unit cell will form a group of its own, etc.
- For each group, pick a maxima that is inside the unit cell (original shape before padding) as the representative. Collection of these representatives is `maxima_pbc`, which is a numpy array of shape (n_maxima_pbc, 3) where n_maxima is the number of representative maxima found in the unit cell. Of course, `n_maxima_pbc = len(pbc_groups)`. You can think of a `pbc_group` as the collection of regions that would form a complete pocket if the pocket is chopped up across periodic boundaries. This is super important for averaging properties correctly.
- Now, we need a way to keep track of all the local maxima and their relation to the `maxima_pbc`. That's where the `maxima_info` array comes in.

#### `maxima_info` array
Let's say we give separate labels 1 to M for all the maxima, including `maxima_pbc` and their images, including the ones in the padded region. You can relate these maxima to `maxima_pbc` using simple offsets.
- Offsets tell you in which directions to shift the representative maxima (`maxima_pbc`) to get the labeled maxima. This is tabulated in `maxima_info`.
- A 2D array of shape `(N_all_maxima, 7)`:

| Column         | Description                                      |
|----------------|--------------------------------------------------|
| `maxima_label` | Label of the maxima                              |
| `pbc_x/y/z`    | Coordinates of the corresponding PBC maxima      |
| `offset_x/y/z` | Grid shifts to locate the periodic image         |

- To compute the actual indices of a maxima `k`:
  ```python
  maxima_indices = maxima_info[k-1, 1:4] + maxima_info[k-1, 4:7] * np.array(dgrid.shape)
  ```

#### Tracking maxima location with respect to unit cell
Note that although a `maxima_pbc` is inside the unit cell, its periodic image could be in the padded region, sometimes close to the unit cell wall but slightly outside the unit cell (since the grid spacing is not infinitesimal). We need to keep track of this information as well.

Also, we don't often need all the maxima from the padded grid; pockets and connections fully inside or just close outside the unit cell are enough for most applications. So it would be useful to trim them down quickly. Hence, this information is provided among the outputs:

- `close_or_inside_maxima`: 1D numpy array that lists the indices of all the maxima that are fully inside the unit cell or close to the wall on the outside.
- `unique_maxima`: 1D numpy array that lists the indices of all the maxima that have at least one voxel assigned to it inside the unit cell. Now this could just be noise from the segmentation, or it could be a real physical influence where there is a pocket segment in the unit cell, but its maxima is slightly outside the unit boundary (hard to visualize, I know! But imagine segments along large smooth channels, for example).

#### Clever indexing scheme
There is a clever indexing scheme put in place to make the selection of these maxima easier. Maxima are reindexed in the order of decreasing importance:

1. `maxima_pbc` get indices `1, 2, ..., N_maxima_pbc`
2. `close_or_inside_maxima` get the next indices if not already assigned
3. `unique_regions` maxima influencing the unit cell get the next indices
4. Maxima in the padded regions connected to `close_or_inside_maxima` get indices
5. Remaining maxima in `all_maxima` (a.k.a. unlabeled ones in the padded regions) get indices

This indexing convention simplifies downstream analyses greatly. For example, to filter maxima corresponding to `unique_regions`, you can do:
```python
maxima = all_maxima[:regions.max()] # or 
maxima = all_maxima[:(len(unique_regions)+1)]
```

#### Segmentation(s)
Note that some of the information above requires you to know the segmentation. So that is indeed done first, and then the maxima are reindexed to match the above scheme.

Segmentation is done twice on the padded distance field:
- Once with `maxima_pbc` markers: the same marker (say, `i`) is used for a `maxima_pbc` and all its images. So you'd get that wrap-around effect for regions across periodic boundaries.
- Once with `all_maxima`: separate markers are used for `maxima_pbc` and its periodic images, so you would get different labels at the 8 corner regions, even if these maxima belong to the same `pbc_group`. This is necessary for computing connections to construct the pore graph.

In both cases, you would get a segmentation array of shape (nx+2*padding, ny+2*padding, nz+2*padding), where the padding is the amount of padding you added to the distance field. From these, you would just grab the central region corresponding to the original unit cell to get `regions_pbc` and `regions`, respectively.

#### Connections
Not going to say much about this here, since it's pretty much the same as in the old workflow, except that we compute them on the `regions_padded` instead of `regions`, so we can get the connections in the padded regions as well. If you don't need these, you have a few options:

- Trim down the `connections`: `[c for c in connections if np.all(np.isin(c[0], unique_regions))]`, meaning keep only those connections where both regions are in the unit cell. You can customize using `np.any` or using `close_or_inside_maxima` instead of `unique_regions` as you need.
- You can also ignore the connections from the periodic segmentation workflow and instead just compute them fresh using the unpadded `regions` array and the distance field, like in the old workflow.
- Both of these should give the same results.

## Periodic Segmentation Workflow Diagram

The following diagram illustrates the complete workflow of the periodic segmentation algorithm, showing how the different steps connect to produce the final outputs:

```mermaid
graph TD
    subgraph Input
        A[Distance Grid] --> C[Padding]
        B[ASE Atoms] --> C
    end

    subgraph Maxima Processing
        C --> D[Padded Distance Grid]
        D --> E[Find Local Maxima]
        E --> F[Identify Periodic Images]
        F --> G[Group into pbc_groups]
        G --> H[Select maxima_pbc Representatives]
        H --> I[Create maxima_info Array]
    end

    subgraph Classification & Indexing
        I --> J[Classify Maxima Location]
        J --> K[Reindex by Importance]
        K --> L1[close_or_inside_maxima]
        K --> L2[unique_maxima]
        K --> L3[all_maxima]
    end

    subgraph Segmentation
        L1 --> M1[Segmentation with maxima_pbc]
        L2 --> M1
        L3 --> M2[Segmentation with all_maxima]
        M1 --> N1[regions_pbc]
        M2 --> N2[regions]
    end

    subgraph Final Processing
        N1 --> O[Compute Connections]
        N2 --> O
        O --> P[Final Outputs]
    end

    style Input fill:#f9f,stroke:#333,stroke-width:2px
    style Maxima Processing fill:#bbf,stroke:#333,stroke-width:2px
    style Classification & Indexing fill:#bfb,stroke:#333,stroke-width:2px
    style Segmentation fill:#fbb,stroke:#333,stroke-width:2px
    style Final Processing fill:#fbf,stroke:#333,stroke-width:2px
```

### Key Workflow Steps:

1. **Input Preparation**: Start with distance grid and ASE atoms object
2. **Padding**: Add periodic boundary regions to the distance grid
3. **Maxima Detection**: Find local maxima using h_maxima algorithm
4. **Image Identification**: Identify periodic images of each maxima
5. **Grouping**: Group periodic images into `pbc_groups`
6. **Representative Selection**: Choose one representative per group (`maxima_pbc`)
7. **Information Array**: Create `maxima_info` with offsets and labels
8. **Classification**: Classify maxima by their location relative to unit cell
9. **Reindexing**: Smart reindexing by importance priority
10. **Dual Segmentation**: Perform two segmentations with different marker strategies
11. **Extraction**: Extract central regions from padded grids
12. **Connections**: Compute connections between regions
13. **Final Output**: Generate all output arrays and data structures

## Summary

This document explains the periodic segmentation approach for analyzing porous materials with periodic boundary conditions. The key innovations include:

**Core Problem**: Regular segmentation doesn't handle periodic boundary conditions properly, leading to inconsistent maxima placement across unit cell boundaries and making downstream analysis difficult.

**Solution**: Periodic segmentation addresses this by:
1. **Padding**: Adding wrapped boundary regions to allow the distance field to see across periodic boundaries
2. **Grouping**: Identifying and grouping periodic images of maxima together (`pbc_groups`)
3. **Representative Selection**: Choosing one representative maxima per group (`maxima_pbc`) for consistent analysis
4. **Smart Indexing**: Implementing a hierarchical indexing scheme that prioritizes important maxima

**Key Outputs**:
- `maxima_pbc`: Representative maxima for each periodic group
- `regions_pbc` and `regions`: Segmented regions with and without periodic wrapping
- `maxima_info`: Detailed mapping between all maxima and their representatives
- `pbc_groups`: Groups of periodic image maxima
- Various filtering arrays (`close_or_inside_maxima`, `unique_maxima`) for different analysis needs

**Benefits**: This approach enables proper averaging of properties across periodic boundaries and provides a solid foundation for pore network analysis in crystalline porous materials.

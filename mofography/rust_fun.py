import numpy as np 
import ase
import logging
logger = logging.getLogger(__name__)


def dgrid_from_atoms_rust(ase_atoms, spacing, radii=None):
    import mgr_rust
    from ase.geometry import cellpar_to_cell
    from ase.data import vdw_radii
    if radii is None:
        radii = vdw_radii[ase_atoms.get_atomic_numbers()]
    fpoints = ase_atoms.get_scaled_positions() - 0.5
    cell = cellpar_to_cell(ase_atoms.cell.cellpar())
    nx, ny, nz = (ase_atoms.cell.cellpar()[:3] / spacing).astype(int) + 1
    dgrid = mgr_rust.dgrid_from_atoms(fpoints, cell, radii, int(nx), int(ny), int(nz))
    return dgrid
def egrid_from_atoms_rust(ase_atoms, spacing, forcefield_mixed=None, cutoff=12.5,probe_symbol='He'):
    import mgr_rust
    if forcefield_mixed is None:
        forcefield_mixed = mgr.mix_lorentz_berthelot(ase_atoms, probe_symbol=probe_symbol)
    mat_atoms_repeat, ff_repeat = mgr.repeat_to_meet_cutoff(ase_atoms, forcefield_mixed=forcefield_mixed, cutoff=cutoff)
    eps, sig = ff_repeat[:, 0], ff_repeat[:, 1]
    frac = mat_atoms_repeat.get_scaled_positions() - 0.5  # centre to [-0.5,0.5]
    cell = cellpar_to_cell(mat_atoms_repeat.cell.cellpar()).T
    # cell = mat_atoms_repeat.get_cell().complete().T        # 3×3 column-vector matrix
    nx, ny, nz = (mat_atoms_repeat.cell.cellpar()[:3] / spacing).astype(int) + 1
    egrid = mgr_rust.egrid_from_atoms(frac, cell, eps, sig, nx, ny, nz, cutoff = cutoff)
    return egrid

def do_periodic_segmentation_rust( dgrid, mat_atoms, mask_thickness=0.5, h=0.2, min_distance=1.5, threshold_abs=1.5, return_padded_regions=False, spacing=None, minimum_window_separation=None, wall_thickness=0.05, extra_padding=5, compactness=0, min_region_size=0.001, disable_tqdm=False, compute_connections=True, return_dgrid = False, ):
    import numpy as np
    # from skimage.segmentation import watershed
    from ase.geometry import cellpar_to_cell
    import copy
    # from ase.geometry import complete_cell
    # import sparse
    # from skimage.segmentation import find_boundaries
    # from concurrent.futures import ThreadPoolExecutor
    from time import time
    # from tqdm import tqdm
    import mgr_rust
    start_time = time()
    from mofography import compute_wall_cells, move_maximum_to_peak, process_wall_maxima, is_inside_or_close_to_boundary, get_fractional_coordinates, interpolate_labels, find_touching_region_pairs, find_end_channel_pairs, merge_small_regions, merge_small_regions_pbc, find_pore_limiting_radius, calculate_padding_thickness, detect_periodic_maxima_direct, get_pbc_groups
    # print(f"timer started")
    logger.info(f"timer started")
    # h = 0.11
    # note: uncomment these imports when refactoring
    # from .helper import find_touching_region_pairs,compute_wall_cells, move_maximum_to_peak, process_wall_maxima, is_inside_or_close_to_boundary, get_fractional_coordinates, interpolate_labels
    # smooth = false
    if spacing is None:
        logger.warning('warning: no spacing provided. estimating. provide the dgrid spacing, for more accurate calculation!\n')
        spacing = np.min(mat_atoms.get_cell_lengths_and_angles()[:3]/ (np.array(dgrid.shape)-1))
        logger.info(f"spacing estimated from dgrid and ase_atoms is {spacing}")

    if minimum_window_separation is None:
        minimum_window_separation = spacing * 3

    logger.info(f"h is {h}")
    plr = find_pore_limiting_radius(dgrid)

    logger.info(f"time(s): {np.round(time() - start_time, decimals=2)} pore limiting radius is {np.round(plr, decimals=2)}\n")
    if (plr-0.25) < mask_thickness:
        logger.warning('warning: mask thickness is too high, resetting to pore limiting radius - 0.25\n')

    mask_thickness = np.min([plr - 0.25, mask_thickness])
    # min_distance = 1.1
#     threshold_abs = 1
# # merged_maxima, all_maxima = detect_periodic_maxima_dgrid(dgrid, ase_atoms=mat_atoms, spacing=spacing, h=h)

    # dt, merged_maxima, all_maxima =detect_periodic_maxima_direct(
    #     dgrid,
    #     mat_atoms,
    #     spacing,
    #     h=h,
    #     threshold_abs=threshold_abs,
    #     min_distance=min_distance,
    #     # smooth=smooth,
    #     maxima_mask=mask_thickness,
    # )

    cell = cellpar_to_cell(mat_atoms.cell.cellpar()).T
    # merged_maxima, unmerged_maxima, raw_maxima = mgr_rust.detect_maxima_h_morph(dgrid, spacing=spacing, h=h, threshold_abs=threshold_abs, min_distance=min_distance, sigma_smooth=1.0, cell=cell)
    merged_maxima, unmerged_maxima, raw_maxima = mgr_rust.detect_maxima_selective_smoothing_parallel(dgrid, preserve_threshold=3.0, smooth_threshold=6.0, max_smooth_radius=2.0, h=h, threshold_abs=threshold_abs, min_distance=min_distance, cell=cell)
    dt = dgrid.copy()
    temp = np.ceil(merged_maxima)
    # compute wall cells
    logger.info(f"time(s): {np.round(time() - start_time, decimals=2)} computing maxima with pbc\n")
    #computed merged maxima nad regions_pbc!
    # padding_thickness=5
    padding_thickness = calculate_padding_thickness(dgrid, mat_atoms, spacing) + np.ceil(extra_padding/spacing).astype(int)
    logger.info(f"padding thickness is {padding_thickness}")
    wall_cells = compute_wall_cells(mat_atoms, spacing, wall_thickness)
    wall_maxima_to_remove = process_wall_maxima(
        merged_maxima,
        dgrid,
        padding_thickness=padding_thickness,
        wall_cells=wall_cells,
        ase_atoms=mat_atoms,
    )
    wall_maxima_to_remove = []
    #logger.info("maxima to remove are", wall_maxima_to_remove)
    # optimized_positions = []
    # optimized_values = []
    # original_indices = []
    # for i, pos in enumerate(merged_maxima):
    #     if i not in wall_maxima_to_remove:
    #         # padded_pos = pos + padding_thickness
    #         opt_pos, opt_val = move_maximum_to_peak(pos, dt, padding_thickness)
    #         optimized_positions.append(opt_pos)
    #         optimized_values.append(opt_val)
    #         original_indices.append(i)
    # optimized_positions = np.vstack(optimized_positions)
    # optimized_positions.shape
    # use the rust functions instead
    optimized_positions, optimized_values = mgr_rust.move_maxima_to_peaks(dt.astype(np.float32), merged_maxima.astype(np.uintp))
    optimized_positions = np.asarray(optimized_positions)

    all_removals = set(wall_maxima_to_remove)
    merged_maxima = np.vstack(
        [m for i, m in enumerate(merged_maxima) if i not in all_removals]
    )


    logger.info(f"time(s): {np.round(time() - start_time, decimals=2)} computing regions with pbc \n")
    # compactness = 0
    # compute markers for padded, periodic segmentation
    markers = np.zeros_like(dgrid, dtype=int)
    if len(merged_maxima) > 0:
        maxima_coords = tuple(merged_maxima.T)
        markers[maxima_coords] = np.arange(1, len(merged_maxima) + 1)
    else:
        logger.info("problem: no valid maxima detected after merging.")
    # void mask for segmentation to use the mask_thickness
    void_mask = dgrid > mask_thickness
    padded_grid = np.pad(dt, padding_thickness, mode="wrap")
    padded_markers = np.pad(markers, padding_thickness, mode="wrap")
    padded_void = np.pad(void_mask, padding_thickness, mode="wrap")

    # # run the watershed segmentation using the merged markers.
    # padded_regions = watershed(
    #     -padded_grid,
    #     padded_markers,
    #     mask=padded_void,
    #     compactness=compactness,
    #     connectivity=2,  # 26-connectivity for full 3d
    # )

    # run the watershed segmentation using the merged markers.
    padded_regions = mgr_rust.seeded_watershed_bucket(
        padded_grid.astype(np.float32),
        padded_markers.astype(np.uint32),
        mask=padded_void.astype(np.uint8),
        # # compactness=compactness,
        # connectivity=2,  # 26-connectivity for full 3d
    )


    # clean up the regions_pbc from small regions, remove less than 0.1% of the total size, this is trickier than one would think. We want to find the volume % of just the box, then make our selection to remove, and then pad back up to get the cleaned padded grid.
    padded_regions_box_trimmed = padded_regions[padding_thickness:-padding_thickness, padding_thickness:-padding_thickness, padding_thickness:-padding_thickness]
    padded_regions_box_trimmed_cleaned = merge_small_regions(padded_regions_box_trimmed, min_size=padded_regions_box_trimmed.size*min_region_size)
    padded_regions_cleaned = np.pad(padded_regions_box_trimmed_cleaned, padding_thickness, mode='wrap')
    # padded_regions_cleaned = merge_small_regions(padded_regions, min_size=padded_regions.size*min_region_size)
    unique_padded_regions= np.unique(padded_regions)
    unique_padded_regions_cleaned = np.unique(padded_regions_cleaned)
    removed_labels_pbc =unique_padded_regions[~np.isin(unique_padded_regions,unique_padded_regions_cleaned)]
    # only keep removed labels that correspond to original maxima
    max_label = len(merged_maxima)
    removed_maxima_labels = removed_labels_pbc[(removed_labels_pbc >= 1) & (removed_labels_pbc <= max_label)]
    logger.info(f"removing these small maxima labels{removed_maxima_labels}")
    # remove those rows from merged_maxima
    merged_maxima = np.delete(merged_maxima, removed_maxima_labels - 1, axis=0)

    # replace regions_pbc with regions_pbc_cleaned
    padded_regions=padded_regions_cleaned

    # extract central (unpadded) region to obtain regions_pbc.
    regions_pbc = padded_regions[
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
    ]
    # print("computed merged maxima nad regions_pbc!!")
    logger.info('removing the small regions, with counts < 0.001 * regions_pbc.size')

    #
    logger.info(f"time(s): {np.round(time() - start_time, decimals=2)} computing maxima and regions without pbc\n")

# step 2: get unmerged maxima info using valid maxima
    original_shape = np.array(dgrid.shape)
    padded_shape = np.array(padded_grid.shape)

    new_maxima_info_list = []
    counter = 1
    for m in merged_maxima:
            m_pad = m + padding_thickness
            offset_ranges = [
                range(-((p - 1) // o), (p - 1) // o + 1)
                for p, o in zip(padded_shape, original_shape)
            ]
            for offs in np.array(np.meshgrid(*offset_ranges)).T.reshape(-1, 3):
                candidate = m_pad + offs * original_shape
                if np.all(candidate >= 0) and np.all(candidate < padded_shape):
                    candidate_unpadded = candidate - padding_thickness
                    offset = candidate_unpadded // original_shape
                    image_coord = np.mod(candidate_unpadded, original_shape)
                    new_maxima_info_list.append(
                        [
                            counter,
                            image_coord[0],
                            image_coord[1],
                            image_coord[2],
                            offset[0],
                            offset[1],
                            offset[2],
                        ]
                    )
                    counter += 1


    maxima_info = np.array(new_maxima_info_list,dtype=int)

    # rust version goes here--
    # maxima_info = mgr_rust.build_maxima_info(merged_maxima.astype(np.uintp), np.array(dgrid.shape, dtype=np.uintp), int(padding_thickness))
    # maxima_info = np.asarray(maxima_info, dtype=int)

    # print(maxima_info)
    # merged_maxima
# create unmerged markers and regions using valid maxima
    unmerged_markers = np.zeros_like(padded_grid, dtype=int)
    for info in maxima_info:
        # warn: should we add the offset times shape here?  fixed
        pos = info[1:4] + padding_thickness + info[4:7] * original_shape
        unmerged_markers[tuple(pos)] = info[0]

# # compute unmerged regions
#     unmerged_padded_regions = watershed(
#         -padded_grid,
#         unmerged_markers,
#         mask=padded_void,
#         compactness=compactness,
#         connectivity=2,
#     )

# compute unmerged regions
    unmerged_padded_regions = mgr_rust.seeded_watershed_bucket(
        padded_grid.astype(np.float32),
        seeds = unmerged_markers.astype(np.uint32),
        mask=padded_void.astype(np.uint8),
        # compactness=compactness,
        # connectivity=2,
    )


# grab the central block
    unmerged_regions = unmerged_padded_regions[
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
    ]

    logger.info(f"time(s): {np.round(time() - start_time, decimals=2)} computed regions without pbc\n")
    endchannel_pairs = find_end_channel_pairs(regions_pbc, unmerged_regions, size_ratio=0.3)
    for p in endchannel_pairs:
        # print(f"replacing {p[0]} with {p[1]}")
        unmerged_regions[unmerged_regions == p[0]] = p[1]
    

    # print('looking for noisy regions, counts < 0.001 * regions.size')
    # unmerged_regions_cleaned =  merge_small_regions(unmerged_regions, min_size=unmerged_regions.size * 0.001)
    # unique_regions = np.unique(unmerged_regions)
    # unique_regions_cleaned = np.unique(unmerged_regions_cleaned)
    # # find the labels that were removed, if any.
    # removed_labels = unique_regions[~np.isin(unique_regions, unique_regions_cleaned)]
    # # remove the corresponding maxima from maxima_info array
    # maxima_info = np.delete(maxima_info, removed_labels - 1, axis=0)
    # # finally replace unmerged regions with the cleaned version
    # unmerged_regions = unmerged_regions_cleaned
    #
    # set_trace()
    logger.info(f"time(s):', {np.round(time() - start_time, decimals=2)} ,'reindexng\n")

# note: we need to reindex here,
# index mapping for the merged_maxima in unmerged parts
# ensure merged maxima indices run from 1 to n
    wall_region_thicknesses = compute_wall_cells(mat_atoms, spacing, wall_thickness=wall_thickness)
# let's start by keeping only those maxima that are there in the regions array
    maxima_indices = maxima_info[:, 1:4] + \
        maxima_info[:, 4:7] * np.array(dgrid.shape)
# get location classifications for all maxima
    inside_close_outside = is_inside_or_close_to_boundary(
        maxima_indices, dgrid.shape, wall_region_thicknesses
    )
# create a boolean mask for maxima that are either inside or close_outside
    valid_maxima_mask = np.array(
        [loc in ["inside", "close_outside"] for loc in inside_close_outside]
    )
# find the core maxima ids (region labels), don't forget to add 1 to get label from index
    valid_maxima_ids = np.where(valid_maxima_mask)[0] + 1
    # print("computed maxima close to the cell")
    valid_maxima_indices = maxima_indices[valid_maxima_mask]
# reindexing logic
    merged_maxima_indices = {tuple(m): i + 1 for i, m in enumerate(merged_maxima)}
    logger.info('constructing remap dict\n')
    used_indices = set()
    remap_dict = {}

    # note: assign indices to periodic maxima
    for info in maxima_info:
        original_maxima = tuple(info[1:4])
        if np.all(info[4:7] == 0) and original_maxima in merged_maxima_indices:
            # remap_dict[info[0]] = merged_maxima_indices[original_maxima]
            # if info[0] in remap_dict and remap_dict[info[0]] != merged_maxima_indices[original_maxima]:
            target_index = merged_maxima_indices[original_maxima]
            # check if the target index is already used
            if target_index in used_indices:
                logger.info(f"conflict: maxima {info[0]} maps to already-used index {target_index}")
            else:
                remap_dict[info[0]] = target_index
                used_indices.add(target_index)

            # remap_dict[info[0]] = merged_maxima_indices[original_maxima]

    # print('stage 1 maxima_pbc', remap_dict)
    unique_unmerged_regions = np.delete(np.unique(unmerged_regions), 0)
    # print(
    #     "all closeby maxima are in unmerged regions",
    #     np.all(np.isin(valid_maxima_ids, unique_unmerged_regions)),
    # )
    # note: assign indices to valid maxima
    # new_index = len(merged_maxima) + 1
    next_index = max(used_indices) + 1 if used_indices else 1
    for vid in valid_maxima_ids:
        if vid not in remap_dict:
            remap_dict[vid] =next_index
            used_indices.add(next_index)
            next_index += 1
    # print('stage 2 valid_maxima', remap_dict)

    # note: assign indices to remaining unique regions in unmerged_regions
    next_index = max(used_indices) + 1 if used_indices else 1
    for r in unique_unmerged_regions:
        if r not in remap_dict:
            remap_dict[r] =next_index
            used_indices.add(next_index)
            next_index += 1
    # print('stage 3 unique_regions', remap_dict)

    # note: assign indices to remaining unique regions in unmerged_regions
    unique_unmerged_padded_regions = np.delete(
        np.unique(unmerged_padded_regions), 0)
    # new_index = len(unique_unmerged_regions) + 1
    next_index = max(used_indices) + 1 if used_indices else 1

    # we need to reindex the padded regions after computing connections
    padding_indices_start_here = next_index
    for r in unique_unmerged_padded_regions:
        if r not in list(remap_dict.keys()):
            remap_dict[r] =next_index
            used_indices.add(next_index)
            next_index += 1
    # print('stage 4 unique_unmerged_padded_regions', remap_dict)

    # print('remap_dict constructed')
    maxima_info_copy = copy.deepcopy(maxima_info)

# great! now we have a complet map. time to remap the maxima_info, unmerged_regions, and unmerged_padded_regions
    regions = copy.deepcopy(unmerged_regions)
    regions_padded= copy.deepcopy(unmerged_padded_regions)
    for k, v in remap_dict.items():
        # print(f"replace label {k} with {v}")
        maxima_info_copy[k - 1][0] = v
        regions[unmerged_regions == k] = v
        regions_padded[unmerged_padded_regions == k] = v
# sort maxima based on first column
    maxima_info_copy = maxima_info_copy[maxima_info_copy[:, 0].argsort()]
# maxima_info[:10, :], maxima_info_copy[:10, :]
    logger.info(f"time(s): {np.round(time() - start_time, decimals=2)}.  done with reindexing-1/2\n")
    # print('remap_dict is ', remap_dict)

    pbc_groups = [v for k, v in get_pbc_groups(regions, maxima_info_copy).items()]
    len(pbc_groups[0])
# find touching pairs in the regions padded
    # pairs = find_touching_region_pairs(regions_padded)
    # print(pairs, pairs.shape)
    wall_region_thicknesses = compute_wall_cells(mat_atoms, spacing, wall_thickness=wall_thickness)
# let's start by keeping only those maxima that are there in the regions array
    maxima_indices = maxima_info_copy[:, 1:4] + maxima_info_copy[:, 4:7] * np.array(
        dgrid.shape
    )
# get location classifications for all maxima
    inside_close_outside = is_inside_or_close_to_boundary(
        maxima_indices, dgrid.shape, wall_region_thicknesses
    )
# create a boolean mask for maxima that are either inside or close_outside
    valid_maxima_mask = np.array(
        [loc in ["inside", "close_outside"] for loc in inside_close_outside]
    )
# lets find the core maxima ids (region labels), dont' forget to add 1 to get label from index
    valid_maxima_ids = np.where(valid_maxima_mask)[0] + 1
    # print("computed maxima close to the cell")
    # pairs_filtered = [p for p in pairs if np.any(np.isin(p, valid_maxima_ids))]


   # --- [copy everything up to and including maxima_info_copy = maxima_info_copy[maxima_info_copy[:, 0].argsort()] ] ---
    # (Paste the code from the original function up to the first reindexing, just before pbc_groups = ...)
    # --- End of copy ---

    # Only compute connections and further reindexing if requested
    if not compute_connections:

        all_maxima = maxima_info_copy[:, 1:4] + maxima_info_copy[:, 4:7] * np.array(dgrid.shape)
        logger.info(f"time(s): { np.round(time() - start_time, decimals=2)} compute_connections is set to False. No connections computed. Calculations Ending! :) \n")
        outputs = dict(
            # dgrid=dgrid,
            all_maxima = all_maxima,
            maxima_pbc=merged_maxima,
            regions_pbc=regions_pbc,
            maxima_info=maxima_info_copy,
            regions=regions,
            pbc_groups = pbc_groups,
            close_or_inside_maxima=valid_maxima_indices,
            close_or_inside_regions=valid_maxima_ids,
            unique_regions=np.delete(np.unique(regions), 0),
            maxima_pbc_dict=merged_maxima_indices,
            maxima_relative_locations=np.array(inside_close_outside))
        outputs['input_params'] = {
            'mask_thickness': mask_thickness,
            'h': h,
            'min_distance': min_distance,
            'threshold_abs': threshold_abs,
            'spacing': spacing,
            'minimum_window_separation': minimum_window_separation,
            'compactness': compactness,
            'wall_thickness': wall_thickness,
            'extra_padding': extra_padding,
            'padding_thickness': padding_thickness,
            'min_region_size': min_region_size,
        }
        if return_padded_regions:
            outputs['regions_padded'] = regions_padded
    else:
        logger.info(f"time(s): { np.round(time() - start_time, decimals=2)} compute_connections is True--this might take a while! \n")

        # def compute_boundary(reg):
        #     return reg, sparse.COO(
        #         find_boundaries(regions_padded == reg, mode="outer").astype(int)
        #     )


# # find uniques and remove background use the reindexed regions here
#         regions_to_check = np.delete(np.unique(regions_padded), 0)
#         with ThreadPoolExecutor() as executor:
#             results = list(executor.map(compute_boundary, regions_to_check))
#         outer = dict(results)
# # don't screw this line up

        # A_unit = complete_cell(mat_atoms.cell).T
        A_unit = cellpar_to_cell(mat_atoms.cell.cellpar()).T
        minimum_window_separation = minimum_window_separation
        # padding_thickness = mgr.calculate_padding_thickness(dgrid, mat_atoms, spacing)
        padded_grid = np.pad(dgrid, padding_thickness, mode="wrap")
        # dinterp = interpolate_labels(dgrid)

        logger.info(f"TIME(s): {np.round(time() - start_time, decimals=2)} Adjusting for periodic connections\n")


        connections_list = mgr_rust.window_connections(regions_padded.astype(np.uint32), padded_grid.astype(np.float32), padding_thickness,A_unit, 1.5)
        connections_list = [np.array([c[0], c[1], c[2], c[3], c[4]], dtype="object") for c in connections_list] 
# let's mark a connection as periodic if any maxima is far from the box, no reindexing just yet.
        for c in connections_list:
            if np.any(~np.isin(c[0], valid_maxima_ids)):
                # print(f"This connection is {c[0]}")
                c[2] = "periodic"
        logger.info(f"TIME(s): {np.round(time() - start_time, decimals=2)} Done with computing connections\n")
        logger.info(f"TIME(s): {np.round(time() - start_time, decimals=2)} Reindexing part 2/2-- Connected nodes getting priority\n")

# NOTE: if padded regions are in the connections list, we need to prioritize those in the indexing over other regions in regions padded
# make a remap dictionary
        remap_dict2 = {}
        unique_regions_padded = np.delete(np.unique(regions_padded), 0)
# here we removed 0 for clarity
        unique_regions = np.delete(np.unique(regions), 0)
        # new_index = len(unique_regions) + 1  # can also do regions.max() + 1
        next_index = padding_indices_start_here
        unique_connected_indices = np.unique( np.vstack([conn[0] for conn in connections_list]))
        for c in unique_connected_indices:
            # if c not in unique_regions:  # this is node that needs priority
            if c >= padding_indices_start_here:  # this is node that needs priority
                remap_dict2[c] =next_index
                next_index+= 1
        for r in unique_regions_padded:
            if np.logical_and(r not in remap_dict2, r >= padding_indices_start_here):
                remap_dict2[r] = next_index
                next_index+= 1
        # print('remap dict 2 is', remap_dict2)
        # set_trace()
        regions_padded_copy = copy.deepcopy(regions_padded)
        for k, v in remap_dict2.items():
            # print(f"replace label {k} with {v}")
            regions_padded_copy[regions_padded == k] = v
        connections_list_copy = copy.deepcopy(connections_list)
        for i, c in enumerate(connections_list):
            for j, ci in enumerate(c[0]):
                if ci in remap_dict2:
                    # reindx the connections list
                    connections_list_copy[i][0][j] = remap_dict2[ci]
        maxima_info_copy_copy = copy.deepcopy(maxima_info_copy)
        for k, v in remap_dict2.items():
            # print(f"replace label {k} with {v}")
            maxima_info_copy_copy[k - 1][0] = v
        maxima_info_copy_copy = maxima_info_copy_copy[maxima_info_copy_copy[:, 0].argsort(
        )]

        # outputs=dict(maxima_pbc=merged_maxima, regions_pbc=regions_pbc, maxima_info=maxima_info_copy_copy, regions=regions,conections=connections_list_copy, pbc_groups=pbc_groups, close_or_inside_maxima=valid_maxima_indices, close_or_inside_regions= valid_maxima_ids, unique_regions = unique_regions)
        all_maxima = maxima_info_copy_copy[:, 1:4] + maxima_info_copy_copy[:, 4:7] * np.array(dgrid.shape)
        outputs=dict(dgrid=dgrid, maxima_pbc=merged_maxima, regions_pbc=regions_pbc, maxima_info=maxima_info_copy_copy, regions=regions,connections=connections_list_copy, pbc_groups=pbc_groups, close_or_inside_maxima=valid_maxima_indices, close_or_inside_regions= valid_maxima_ids, unique_regions = unique_regions, maxima_pbc_dict=merged_maxima_indices, all_maxima=all_maxima, maxima_relative_locations=np.array(inside_close_outside), )
        # we need to save the input params to the results, for future reference.
        outputs['input_params'] = {
            'mask_thickness': mask_thickness,
            'h': h,
            'min_distance': min_distance,
            'threshold_abs': threshold_abs,
            'spacing': spacing,
            'minimum_window_separation': minimum_window_separation,
            'compactness': compactness,
            'wall_thickness': wall_thickness,
            'extra_padding': extra_padding,
            'padding_thickness': padding_thickness,
            'min_region_size': min_region_size,
        }
        if return_padded_regions:
            outputs['regions_padded'] = regions_padded_copy
            
    if return_dgrid:
        outputs['dgrid'] = dgrid
    return outputs



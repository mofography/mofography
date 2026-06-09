#!/usr/bin/env python3
"""Test blended-priority watershed (watershed_blended crate) vs rust_det / rust_det_metric.

Usage:
    python test_blended_watershed.py                                    # NU-1000, blended
    python test_blended_watershed.py --backend rust_det                 # unweighted baseline
    python test_blended_watershed.py --backend rust_det_metric          # old metric (dgrid*inv_w)
    python test_blended_watershed.py --backend rust_det_blended         # new blended priority
    python test_blended_watershed.py --backend rust_det_blended --scale 100
    python test_blended_watershed.py --cif examples/TestCIFs/VIGNOQ_clean.cif --html
"""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path
from time import time

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from ase.geometry import cellpar_to_cell
from ase.io import read

import mofography as mgr
import mgr_rust
import watershed_blended as wb

from mofography import (
    compute_wall_cells,
    is_inside_or_close_to_boundary,
    find_end_channel_pairs,
    merge_small_regions,
    find_pore_limiting_radius,
    calculate_padding_thickness,
    get_pbc_groups,
    trim_connections_to_unique_regions,
)
from mofography.rust_fun import (
    _cell_matrix_for_metric_weights,
    _voxel_lengths_for_smoothing,
    _compute_metric_inv_weights,
    dgrid_from_atoms_rust,
)

import logging

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent

# ── helpers ──────────────────────────────────────────────────────────


def _label_stats(arr, name):
    labels = np.unique(arr)
    nz = labels[labels != 0]
    print(f"  {name}: shape={arr.shape}, max={arr.max()}, unique(non0)={len(nz)}")
    return nz


def _check_indexing(regions, maxima_info, label=""):
    mi = np.asarray(maxima_info, dtype=int)
    in_cell = np.all(mi[:, 4:7] == 0, axis=1)
    mismatches = 0
    for row in mi[in_cell]:
        pos = tuple(row[1:4])
        if all(0 <= p < s for p, s in zip(pos, regions.shape)):
            if regions[pos] != row[0]:
                mismatches += 1
                print(
                    f"  MISMATCH {label}: pos={pos}, mi_label={row[0]}, regions={regions[pos]}"
                )
    print(f"  Indexing mismatches ({label}): {mismatches}")
    return mismatches


# ── watershed dispatch ───────────────────────────────────────────────


def _run_watershed(padded_grid, padded_markers, padded_void, backend, **kw):
    """Dispatch watershed to the correct backend."""
    img = padded_grid.astype(np.float32)
    seeds = padded_markers.astype(np.uint32)
    mask = padded_void.astype(np.uint8) if padded_void is not None else None

    if backend == "rust_det":
        return mgr_rust.seeded_watershed_det(img, seeds, mask)

    elif backend == "rust_det_metric":
        inv_w = kw.get("inv_metric_weights")
        if inv_w is None:
            raise ValueError("rust_det_metric requires inv_metric_weights")
        return mgr_rust.seeded_watershed_det_metric(img, seeds, inv_w, mask)

    elif backend == "rust_det_blended":
        edge_lengths = kw.get("edge_lengths")
        scale = kw.get("scale", 200.0)
        if edge_lengths is None:
            raise ValueError("rust_det_blended requires edge_lengths")
        return wb.seeded_watershed_det_blended(img, seeds, edge_lengths, scale, mask)

    else:
        raise ValueError(f"Unknown backend: {backend}")


# ── main pipeline (verbatim copy of do_periodic_segmentation_rust) ──


def do_periodic_segmentation_blended(
    dgrid,
    mat_atoms,
    *,
    backend="rust_det_blended",
    scale=200.0,
    spacing=0.5,
    mask_thickness=0.05,
    h=0.2,
    min_distance=1.5,
    threshold_abs=1.5,
    preserve_threshold=3.0,
    smooth_threshold=6.0,
    max_smooth_radius=2.0,
    extra_padding=1.0,
    compactness=0,
    min_region_size=0.0,
    watershed_connectivity=2,
    skip_endchannel_pairs=False,
    compute_connections=True,
    return_padded_regions=True,
):
    """Full periodic segmentation with selectable watershed backend."""

    start_time = time()
    if spacing is None:
        spacing = np.min(
            mat_atoms.get_cell_lengths_and_angles()[:3] / (np.array(dgrid.shape) - 1)
        )
    minimum_window_separation = spacing * 3

    plr = find_pore_limiting_radius(dgrid)
    print(f"PLR = {plr:.4f} A")

    if (plr - 0.25) < mask_thickness:
        print(f"WARNING: mask_thickness too high, resetting")
    mask_thickness = np.min([plr - 0.25, mask_thickness])

    cell = cellpar_to_cell(mat_atoms.cell.cellpar()).T
    voxel_lengths = _voxel_lengths_for_smoothing(mat_atoms, spacing, grid_shape=dgrid.shape)

    merged_maxima, unmerged_maxima, raw_maxima = (
        mgr_rust.detect_maxima_selective_smoothing_parallel(
            dgrid,
            preserve_threshold=preserve_threshold,
            smooth_threshold=smooth_threshold,
            max_smooth_radius=max_smooth_radius,
            h=h,
            threshold_abs=threshold_abs,
            min_distance=min_distance,
            cell=cell,
            voxel_lengths=voxel_lengths,
        )
    )
    print(f"\n=== Step 1: Maxima detection ({backend}) ===")
    print(f"  raw={len(raw_maxima)}, unmerged={len(unmerged_maxima)}, merged={len(merged_maxima)}")

    # ── Precompute metric weights / edge lengths ──
    cell_for_metric = _cell_matrix_for_metric_weights(mat_atoms)
    inv_metric_weights = _compute_metric_inv_weights(
        cell_matrix=cell_for_metric, grid_shape=dgrid.shape, spacing=spacing
    )
    edge_lengths = wb.compute_edge_lengths(cell_for_metric, dgrid.shape)

    if backend == "rust_det_blended":
        el = np.array(edge_lengths)
        print(f"  edge_lengths: range=[{el.min():.4f}, {el.max():.4f}], mean={el.mean():.4f}")
        print(f"  scale={scale}")
    elif backend == "rust_det_metric":
        iw = np.array(inv_metric_weights)
        print(f"  inv_metric_weights: range=[{iw.min():.4f}, {iw.max():.4f}], mean={iw.mean():.4f}")

    dt = dgrid.copy()
    padding_thickness = calculate_padding_thickness(dgrid, mat_atoms, spacing) + np.ceil(
        extra_padding / spacing
    ).astype(int)
    print(f"\n=== Step 2: padding_thickness={padding_thickness} ===")

    wall_cells = compute_wall_cells(mat_atoms, spacing, 0.05)
    wall_maxima_to_remove = []

    optimized_positions, optimized_values = mgr_rust.move_maxima_to_peaks(
        dt.astype(np.float32), merged_maxima.astype(np.uintp)
    )
    optimized_positions = np.asarray(optimized_positions)

    all_removals = set(wall_maxima_to_remove)
    merged_maxima = np.vstack(
        [m for i, m in enumerate(merged_maxima) if i not in all_removals]
    )

    # ── Step 3: PBC flood ──
    markers = np.zeros_like(dgrid, dtype=int)
    if len(merged_maxima) > 0:
        markers[tuple(merged_maxima.T)] = np.arange(1, len(merged_maxima) + 1)
    void_mask = dgrid > mask_thickness
    padded_grid = np.pad(dt, padding_thickness, mode="wrap")
    padded_markers = np.pad(markers, padding_thickness, mode="wrap")
    padded_void = np.pad(void_mask, padding_thickness, mode="wrap")

    ws_kw = {}
    if backend == "rust_det_metric":
        ws_kw["inv_metric_weights"] = inv_metric_weights
    elif backend == "rust_det_blended":
        ws_kw["edge_lengths"] = edge_lengths
        ws_kw["scale"] = scale

    padded_regions = _run_watershed(
        padded_grid, padded_markers, padded_void, backend, **ws_kw
    )

    # ── Step 4: merge_small_regions ──
    padded_regions_box_trimmed = padded_regions[
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
    ]
    padded_regions_box_trimmed_cleaned = merge_small_regions(
        padded_regions_box_trimmed,
        min_size=padded_regions_box_trimmed.size * min_region_size,
    )
    padded_regions_cleaned = np.pad(
        padded_regions_box_trimmed_cleaned, padding_thickness, mode="wrap"
    )

    unique_padded_regions = np.unique(padded_regions)
    unique_padded_regions_cleaned = np.unique(padded_regions_cleaned)
    removed_labels_pbc = unique_padded_regions[
        ~np.isin(unique_padded_regions, unique_padded_regions_cleaned)
    ]
    print(f"\n=== Step 4: merge_small_regions ===")
    print(f"  Removed PBC labels: {removed_labels_pbc.tolist() or 'none'}")

    max_label = len(merged_maxima)
    removed_maxima_labels = removed_labels_pbc[
        (removed_labels_pbc >= 1) & (removed_labels_pbc <= max_label)
    ]
    merged_maxima = np.delete(merged_maxima, removed_maxima_labels - 1, axis=0)

    if len(removed_maxima_labels) > 0:
        print(f"  Removed maxima labels: {removed_maxima_labels.tolist()}, redoing PBC flood...")
        markers = np.zeros_like(dgrid, dtype=int)
        if len(merged_maxima) > 0:
            markers[tuple(merged_maxima.T)] = np.arange(1, len(merged_maxima) + 1)
        void_mask = dgrid > mask_thickness
        padded_grid = np.pad(dt, padding_thickness, mode="wrap")
        padded_markers = np.pad(markers, padding_thickness, mode="wrap")
        padded_void = np.pad(void_mask, padding_thickness, mode="wrap")
        padded_regions = _run_watershed(
            padded_grid, padded_markers, padded_void, backend, **ws_kw
        )

    regions_pbc = padded_regions[
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
    ]

    print(f"\n=== Step 5: After PBC flood ===")
    _label_stats(regions_pbc, "regions_pbc")
    print(f"  merged_maxima: {len(merged_maxima)}")
    for i, m in enumerate(merged_maxima):
        cnt = int(np.sum(regions_pbc == i + 1))
        print(
            f"    label {i+1} at ({m[0]:3d},{m[1]:3d},{m[2]:3d}): dgrid={dgrid[m[0],m[1],m[2]]:.3f} A, {cnt} voxels"
        )

    # ── Step 6: Build maxima_info ──
    original_shape = np.array(dgrid.shape)
    padded_shape = np.array(padded_grid.shape)
    merged_maxima = np.asarray(merged_maxima, dtype=int)
    if merged_maxima.size:
        order = np.lexsort(
            (merged_maxima[:, 2], merged_maxima[:, 1], merged_maxima[:, 0])
        )
        merged_maxima = merged_maxima[order]

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
    maxima_info = np.array(new_maxima_info_list, dtype=int)

    in_cell_mi = np.all(maxima_info[:, 4:7] == 0, axis=1)
    print(f"\n=== Step 6: maxima_info ===")
    print(f"  rows={len(maxima_info)}, in-cell (offset=0)={np.sum(in_cell_mi)}")

    # ── Step 7: Unmerged flood ──
    unmerged_markers = np.zeros_like(padded_grid, dtype=int)
    for info in maxima_info:
        pos = info[1:4] + padding_thickness + info[4:7] * original_shape
        unmerged_markers[tuple(pos)] = info[0]

    unmerged_padded_regions = _run_watershed(
        padded_grid, unmerged_markers, padded_void, backend, **ws_kw
    )
    unmerged_regions = unmerged_padded_regions[
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
    ]

    print(f"\n=== Step 7: Unmerged flood ===")
    _label_stats(unmerged_regions, "unmerged_regions")

    # ── Step 8: endchannel_pairs ──
    if skip_endchannel_pairs:
        endchannel_pairs = []
    else:
        endchannel_pairs = find_end_channel_pairs(
            regions_pbc, unmerged_regions, size_ratio=0.3
        )
    print(f"\n=== Step 8: endchannel_pairs ===")
    print(f"  endchannel_pairs: {len(endchannel_pairs)} pairs")
    for p in endchannel_pairs[:5]:
        print(f"    {p}")

    for p in endchannel_pairs:
        j_max = np.where(maxima_info[:, 0] == p[1])[0][0]
        i_max = np.where(maxima_info[:, 0] == p[0])[0][0]
        wall_region_thicknesses = compute_wall_cells(mat_atoms, spacing, wall_thickness=0.05)
        maxima_indices = maxima_info[:, 1:4] + maxima_info[:, 4:7] * np.array(dgrid.shape)
        inside_close_outside = is_inside_or_close_to_boundary(
            maxima_indices, dgrid.shape, wall_region_thicknesses
        )
        valid_maxima_ids = np.array(
            [
                idx + 1
                for idx, loc in enumerate(inside_close_outside)
                if loc in ["inside", "close_outside"]
            ]
        )
        if p[0] not in valid_maxima_ids:
            print(f"  endchannel relabel: {p[0]} -> {p[1]} (far outside)")
            unmerged_regions[unmerged_regions == p[0]] = p[1]

    # ── Step 9: 4-stage remap ──
    print(f"\n=== Step 9: 4-stage remap ===")
    wall_region_thicknesses = compute_wall_cells(mat_atoms, spacing, wall_thickness=0.05)
    maxima_indices = maxima_info[:, 1:4] + maxima_info[:, 4:7] * np.array(dgrid.shape)
    inside_close_outside = is_inside_or_close_to_boundary(
        maxima_indices, dgrid.shape, wall_region_thicknesses
    )
    valid_maxima_ids = np.array(
        [
            idx + 1
            for idx, loc in enumerate(inside_close_outside)
            if loc in ["inside", "close_outside"]
        ]
    )
    valid_maxima_indices = maxima_indices[
        np.array([loc in ["inside", "close_outside"] for loc in inside_close_outside])
    ]

    merged_maxima_indices = {
        tuple(m): i + 1 for i, m in enumerate(merged_maxima)
    }
    used_indices = set()
    remap_dict = {}

    # Stage 1
    for info in maxima_info:
        original_maxima = tuple(info[1:4])
        if np.all(info[4:7] == 0) and original_maxima in merged_maxima_indices:
            target_index = merged_maxima_indices[original_maxima]
            if target_index not in used_indices:
                remap_dict[info[0]] = target_index
                used_indices.add(target_index)
    print(f"  Stage 1 (pbc maxima): {len(remap_dict)} mappings")

    # Stage 2
    unique_unmerged_regions = np.delete(np.unique(unmerged_regions), 0)
    next_index = max(used_indices) + 1 if used_indices else 1
    s2 = 0
    for vid in valid_maxima_ids:
        if np.sum(unmerged_regions == vid) > 0 and vid not in remap_dict:
            remap_dict[vid] = next_index
            used_indices.add(next_index)
            next_index += 1
            s2 += 1
    print(f"  Stage 2 (close/inside w/ voxels): +{s2}")

    # Stage 3
    next_index = max(used_indices) + 1 if used_indices else 1
    s3 = 0
    for r in unique_unmerged_regions:
        if r not in remap_dict:
            remap_dict[r] = next_index
            used_indices.add(next_index)
            next_index += 1
            s3 += 1
    print(f"  Stage 3 (unique_regions remainder): +{s3}")

    # Stage 4
    unique_unmerged_padded_regions = np.delete(np.unique(unmerged_padded_regions), 0)
    next_index = max(used_indices) + 1 if used_indices else 1
    padding_indices_start_here = next_index
    s4 = 0
    for r in unique_unmerged_padded_regions:
        if r not in list(remap_dict.keys()):
            remap_dict[r] = next_index
            used_indices.add(next_index)
            next_index += 1
            s4 += 1
    print(f"  Stage 4 (padding): +{s4}, padding_start={padding_indices_start_here}")

    # Apply remap
    maxima_info_copy = copy.deepcopy(maxima_info)
    regions = copy.deepcopy(unmerged_regions)
    regions_padded = copy.deepcopy(unmerged_padded_regions)
    for k, v in remap_dict.items():
        maxima_info_copy[k - 1][0] = v
        regions[unmerged_regions == k] = v
        regions_padded[unmerged_padded_regions == k] = v
    maxima_info_copy = maxima_info_copy[maxima_info_copy[:, 0].argsort()]

    print(f"\n=== Step 10: After remap ===")
    _label_stats(regions, "regions (unit cell)")
    in_cell = np.all(maxima_info_copy[:, 4:7] == 0, axis=1)
    print(f"  In-cell mi labels: {maxima_info_copy[in_cell, 0].tolist()}")
    _check_indexing(regions, maxima_info_copy, "remap1")

    pbc_groups = [v for k, v in get_pbc_groups(regions, maxima_info_copy).items()]
    print(f"  pbc_groups: {len(pbc_groups)}, sizes={[len(g) for g in pbc_groups]}")

    # ── Step 11: connections + remap2 ──
    wall_region_thicknesses = compute_wall_cells(mat_atoms, spacing, wall_thickness=0.05)
    maxima_indices = maxima_info_copy[:, 1:4] + maxima_info_copy[:, 4:7] * np.array(
        dgrid.shape
    )
    inside_close_outside = is_inside_or_close_to_boundary(
        maxima_indices, dgrid.shape, wall_region_thicknesses
    )
    valid_maxima_ids = np.array(
        [
            idx + 1
            for idx, loc in enumerate(inside_close_outside)
            if loc in ["inside", "close_outside"]
        ]
    )
    valid_maxima_indices = maxima_indices[
        np.array([loc in ["inside", "close_outside"] for loc in inside_close_outside])
    ]

    if not compute_connections:
        return dict(
            maxima_pbc=merged_maxima,
            regions_pbc=regions_pbc,
            maxima_info=maxima_info_copy,
            regions=regions,
            pbc_groups=pbc_groups,
            unique_regions=np.delete(np.unique(regions), 0),
            maxima_pbc_dict=merged_maxima_indices,
            maxima_relative_locations=np.array(inside_close_outside),
        )

    A_unit = cellpar_to_cell(mat_atoms.cell.cellpar()).T
    padded_grid = np.pad(dgrid, padding_thickness, mode="wrap")
    connections_list = mgr_rust.window_connections(
        regions_padded.astype(np.uint32),
        padded_grid.astype(np.float32),
        padding_thickness,
        A_unit,
        float(minimum_window_separation),
    )
    connections_list = [
        np.array([c[0], c[1], c[2], c[3], c[4]], dtype="object") for c in connections_list
    ]
    unit_cell_region_labels = np.delete(np.unique(regions), 0)
    for c in connections_list:
        if np.any(~np.isin(c[0], unit_cell_region_labels)):
            c[2] = "periodic"

    print(f"\n=== Step 11: Connections ===")
    print(f"  connections: {len(connections_list)}")

    # remap2
    remap_dict2 = {}
    unique_regions_padded = np.delete(np.unique(regions_padded), 0)
    next_index = padding_indices_start_here
    unique_connected_indices = np.unique(
        np.vstack([conn[0] for conn in connections_list])
    )
    for c in unique_connected_indices:
        if c >= padding_indices_start_here:
            remap_dict2[c] = next_index
            next_index += 1
    for r in unique_regions_padded:
        if r not in remap_dict2 and r >= padding_indices_start_here:
            remap_dict2[r] = next_index
            next_index += 1
    print(f"  remap_dict2: {len(remap_dict2)} mappings")

    regions_padded_copy = copy.deepcopy(regions_padded)
    for k, v in remap_dict2.items():
        regions_padded_copy[regions_padded == k] = v
    connections_list_copy = copy.deepcopy(connections_list)
    for i, c in enumerate(connections_list):
        for j, ci in enumerate(c[0]):
            if ci in remap_dict2:
                connections_list_copy[i][0][j] = remap_dict2[ci]
    maxima_info_copy_copy = copy.deepcopy(maxima_info_copy)
    for k, v in remap_dict2.items():
        maxima_info_copy_copy[k - 1][0] = v
    maxima_info_copy_copy = maxima_info_copy_copy[maxima_info_copy_copy[:, 0].argsort()]

    print(f"\n=== Step 12: Final output ({backend}) ===")
    _label_stats(regions, "regions (unit cell)")
    in_cell = np.all(maxima_info_copy_copy[:, 4:7] == 0, axis=1)
    print(f"  In-cell mi labels: {maxima_info_copy_copy[in_cell, 0].tolist()}")
    _check_indexing(regions, maxima_info_copy_copy, "final")

    all_maxima = maxima_info_copy_copy[:, 1:4] + maxima_info_copy_copy[:, 4:7] * np.array(
        dgrid.shape
    )
    return dict(
        dgrid=dgrid,
        mat_atoms=mat_atoms,
        maxima_pbc=merged_maxima,
        regions_pbc=regions_pbc,
        maxima_info=maxima_info_copy_copy,
        regions=regions,
        connections=connections_list_copy,
        pbc_groups=pbc_groups,
        close_or_inside_maxima=valid_maxima_indices,
        close_or_inside_regions=valid_maxima_ids,
        unique_regions=np.delete(np.unique(regions), 0),
        maxima_pbc_dict=merged_maxima_indices,
        all_maxima=all_maxima,
        maxima_relative_locations=np.array(inside_close_outside),
        regions_padded=regions_padded_copy,
    )


# ── HTML visualization ────────────────────────────────────────────────


def write_html(result: dict, mat_atoms, dgrid, out_path: Path) -> None:
    """Write a Plotly HTML visualization from the segmentation result dict."""
    cscale = px.colors.qualitative.Pastel
    isoval = 1.0

    regions_pbc = np.asarray(result["regions_pbc"], dtype=np.uint32)
    regions = np.asarray(result["regions"], dtype=np.uint32)
    maxima_info = np.asarray(result["maxima_info"])
    maxima_pbc = np.asarray(result["maxima_pbc"])
    all_maxima = np.asarray(result["all_maxima"])
    connections = result.get("connections", [])
    unique_regions = result["unique_regions"]
    grid_shape = np.array(dgrid.shape)

    in_cell_rows = []
    for i, row in enumerate(maxima_info):
        pos = (row[1:4] + row[4:7] * grid_shape).astype(int)
        if all(0 <= pos[j] < dgrid.shape[j] for j in range(3)):
            in_cell_rows.append(i)
    in_cell_rows = np.asarray(in_cell_rows, dtype=int)
    maxima_in_cell = all_maxima[in_cell_rows]
    region_labels = maxima_info[in_cell_rows, 0].astype(int)

    cage_radii = dgrid[tuple(np.mod(maxima_in_cell, dgrid.shape).astype(int).T)]
    maxima_hover = [
        f"{lab}  [{' '.join(map(str, maxima_in_cell[i]))}]  {np.round(cage_radii[i], 2)} A"
        for i, lab in enumerate(region_labels)
    ]

    connections_trimmed = trim_connections_to_unique_regions(connections, unique_regions)
    print(f"  connections: {len(connections)} total, {len(connections_trimmed)} trimmed")

    mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
    frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)
    cell_trace = mgr.plot_cell(mat_atoms, color="orange")

    isosurface_pbc = mgr.plot_surface_dgrid(
        dgrid, isovalue=isoval, ase_atoms=mat_atoms,
        region_labels=regions_pbc, opacity=1, cscale=cscale,
        stylize=True, name="Regions PBC",
    )
    isosurface_regions = mgr.plot_surface_dgrid(
        dgrid, isovalue=isoval, ase_atoms=mat_atoms,
        region_labels=regions, opacity=1, cscale=cscale,
        stylize=True, name="Regions",
    )
    maxima_trace_pbc = mgr.plot_maxima(
        dgrid, maxima_pbc, mat_atoms, show_spheres=False,
        marker_size=7.5, marker_colorscale=cscale, name="Maxima (pbc)",
    )
    maxima_trace = mgr.plot_maxima(
        dgrid, maxima_in_cell, mat_atoms, show_spheres=False,
        marker_size=7.5, marker_colorscale=cscale, name="Maxima (in-cell)",
        marker_color=region_labels, hovertext=maxima_hover,
    )

    connections_traces = mgr.plot_connections_as_it_is(
        connections_trimmed, all_maxima, regions.shape, mat_atoms,
        cscale=cscale, plot_disconnected=False, plot_periodic=True,
    )
    if len(connections_traces) >= 2:
        connections_traces[-2].marker.size = 1

    windows_trace = mgr.plot_windows(
        connections_trimmed, name="Windows", marker_color="yellow", marker_size=10,
    )

    fig = go.Figure(
        data=[
            *frame_trace, cell_trace, maxima_trace_pbc, isosurface_pbc,
            *connections_traces, windows_trace, isosurface_regions, maxima_trace,
        ],
        layout=mgr.dark_layout,
    )
    fig.update_layout(scene_aspectmode="data", legend_groupclick="toggleitem")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(out_path))
    print(f"Wrote visualization: {out_path}")


# ── entry point ─────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cif", type=Path, default=ROOT / "examples/TestCIFs/NU-1000_S1.cif"
    )
    parser.add_argument(
        "--backend",
        choices=["rust_det", "rust_det_metric", "rust_det_blended"],
        default="rust_det_blended",
    )
    parser.add_argument("--scale", type=float, default=200.0, help="Blended priority scale (default 200)")
    parser.add_argument(
        "--html", type=Path, nargs="?", const="AUTO", default=None,
        help="Write Plotly HTML (auto-name if flag alone)",
    )
    parser.add_argument("--min-region-size", type=float, default=0.0)
    args = parser.parse_args()

    atoms = read(str(args.cif.expanduser().resolve()))
    radii = mgr.default_vdw_radii_angstrom(atoms)
    dgrid = dgrid_from_atoms_rust(atoms, spacing=0.5, radii=radii).astype(np.float32)

    print(f"{'=' * 72}")
    print(f"CIF: {args.cif.name}")
    print(f"Cell: {atoms.cell.cellpar()}")
    print(f"dgrid shape: {dgrid.shape}")
    print(f"Backend: {args.backend}")
    if args.backend == "rust_det_blended":
        print(f"Scale: {args.scale}")
    print(f"{'=' * 72}")

    result = do_periodic_segmentation_blended(
        dgrid, atoms, backend=args.backend, scale=args.scale,
        spacing=0.5, mask_thickness=0.05, h=0.2,
        extra_padding=1.0, compute_connections=True,
        return_padded_regions=True, min_region_size=args.min_region_size,
    )

    # HTML
    if args.html is not None:
        if args.html == Path("AUTO"):
            stem = args.cif.stem
            html_path = ROOT / f"{stem}_{args.backend}_s{args.scale:.0f}.html"
        else:
            html_path = args.html
        print(f"\n{'=' * 72}")
        print(f"Generating HTML: {html_path}")
        write_html(result, atoms, dgrid, html_path)

    elapsed = time() - start_time if "start_time" in dir() else 0
    print(f"\n{'=' * 72}")
    print("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())

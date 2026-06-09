#!/usr/bin/env python
"""
Standalone test for spill-connection filtering.

For each connection from window_connections, checks if window center positions
are near framework atoms. Removes positions (and connections) that are
floating in empty pore space (spill-over artifacts).

Usage:
    python test_spill_filter.py
"""

import sys
import time
import logging
import numpy as np

logging.basicConfig(level=logging.WARNING)

from ase.io import read as ase_read
from ase.geometry import complete_cell

from mofography.rust_fun import (
    do_periodic_segmentation_rust,
    dgrid_from_atoms_rust,
)
from mofography import default_vdw_radii_angstrom
from mofography.winding_percolation import (
    get_pore_dia_mgr_winding2,
    make_directed_graph_from_segmentation_results,
    build_base_region_labels_from_results,
    find_pore_limiting_radius_winding,
    percolation_all_axes_winding,
    _dfsd_threshold,
    _expand_dfsd_labels_with_window_fingerprint,
    _window_frac_from_edge_data,
)


# ---------------------------------------------------------------------------
# Helper: nearest atom distance (PBC-aware)
# ---------------------------------------------------------------------------

def _precompute_atom_pbc(atoms):
    """Precompute cell matrix and fractional positions for fast nearest-atom lookups."""
    cell_matrix = complete_cell(atoms.get_cell()).T  # (3,3) frac→cart
    atom_frac = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
    return cell_matrix, atom_frac


def _nearest_atom_distances(frac_positions, atom_frac, cell_matrix):
    """Vectorized nearest-atom distance for multiple fractional positions.

    Parameters
    ----------
    frac_positions : (M, 3) float array
    atom_frac : (N, 3) float array
    cell_matrix : (3, 3) float array  (frac→cart)

    Returns
    -------
    (M,) float array of nearest distances
    """
    # frac_positions: (M,1,3) - atom_frac: (1,N,3) → (M,N,3)
    dr_frac = frac_positions[:, None, :] - atom_frac[None, :, :]
    dr_frac -= np.round(dr_frac)
    dr_cart = np.einsum("mni,ij->mnj", dr_frac, cell_matrix)
    dists = np.linalg.norm(dr_cart, axis=2)  # (M, N)
    return dists.min(axis=1)  # (M,)


def _cart_to_frac_padded(cart_pos, cell_matrix, dgrid_shape, padding_thickness, spacing):
    """Convert padded-grid Cartesian coords to wrapped unit-cell fractional coords.

    padded_grid coords: (i*spacing, j*spacing, k*spacing) where i goes from
    -padding_thickness*spacing to (nx+padding_thickness)*spacing.
    Unit cell frac: ((pos/spacing - padding_thickness + 0.5) / dgrid_shape) % 1
    """
    # Convert to padded-grid indices
    grid_idx = np.asarray(cart_pos, dtype=float) / spacing

    # Remove padding offset → unit-cell grid indices (0-based)
    uc_idx = grid_idx - padding_thickness

    # Convert to fractional coordinates
    shape = np.asarray(dgrid_shape, dtype=float)
    frac = (uc_idx + 0.5) / shape
    return np.mod(frac, 1.0)


# ---------------------------------------------------------------------------
# Main filter function
# ---------------------------------------------------------------------------

def filter_spill_connections(
    connections_list,
    results,
    atoms,
    spill_margin=2.0,
    spacing=0.5,
):
    """Remove spill-over windows from connections_list.

    For each connection, examines every window position. If the nearest
    framework atom is farther than ``window_radius + spill_margin``, the
    position is a spill artifact and is removed.

    If all positions are removed, the connection is dropped entirely.

    Parameters
    ----------
    connections_list : list of np.ndarray
        Output of ``mgr_rust.window_connections`` or the processed list
        from ``do_periodic_segmentation_rust``.
    results : dict
        Segmentation results (must contain 'dgrid', 'input_params').
    atoms : ase.Atoms
        Framework atoms.
    spill_margin : float
        Extra tolerance (Å) beyond window radius. Default 2.0.
    spacing : float
        Grid spacing. Default 0.5.

    Returns
    -------
    filtered : list of np.ndarray
        Connections with spill positions removed.
    stats : dict
        Filtering statistics.
    """
    import logging
    logger = logging.getLogger(__name__)

    dgrid = np.asarray(results["dgrid"])
    dgrid_shape = tuple(int(x) for x in dgrid.shape)
    params = results.get("input_params", {})
    pt = int(params.get("padding_thickness", 30))

    cell_matrix, atom_frac = _precompute_atom_pbc(atoms)

    filtered = []
    n_connections_dropped = 0
    n_positions_removed = 0
    n_positions_total = 0

    for conn in connections_list:
        regs = conn[0]
        connected = conn[1]
        type_str = conn[2]
        centers = conn[3]  # (M, 3) Cartesian padded-grid coords
        radii_arr = conn[4]  # (M,) window radii

        if centers.ndim != 2 or centers.shape[0] == 0:
            filtered.append(conn)
            continue

        n_positions_total += centers.shape[0]

        # Convert centers to fractional coordinates
        frac_pos = _cart_to_frac_padded(centers, cell_matrix, dgrid_shape, pt, spacing)

        # Compute nearest atom distances
        nearest_dists = _nearest_atom_distances(frac_pos, atom_frac, cell_matrix)

        # Check each position: keep if nearest atom is within wr + margin
        keep_mask = nearest_dists <= radii_arr + spill_margin

        n_removed = int((~keep_mask).sum())
        n_positions_removed += n_removed

        if not keep_mask.any():
            # All positions removed → drop entire connection
            n_connections_dropped += 1
            logger.info(
                "Spill filter: dropped connection %s (%d positions, "
                "nearest atoms %.1f-%.1f Å, wr %.2f-%.2f Å)",
                regs, centers.shape[0],
                nearest_dists.min(), nearest_dists.max(),
                float(radii_arr.min()), float(radii_arr.max()),
            )
            continue

        if n_removed > 0:
            logger.info(
                "Spill filter: removed %d/%d positions from connection %s "
                "(nearest %.1f-%.1f Å)",
                n_removed, centers.shape[0], regs,
                nearest_dists[~keep_mask].min(), nearest_dists[~keep_mask].max(),
            )

        # Rebuild connection with filtered positions
        new_centers = centers[keep_mask]
        new_radii = radii_arr[keep_mask]
        filtered.append(np.array([
            regs, connected, type_str, new_centers, new_radii
        ], dtype=object))

    stats = {
        "connections_before": len(connections_list),
        "connections_after": len(filtered),
        "connections_dropped": n_connections_dropped,
        "positions_total": n_positions_total,
        "positions_removed": n_positions_removed,
    }

    logger.info(
        "Spill filter summary: %d/%d connections kept, "
        "%d positions removed out of %d total",
        len(filtered), len(connections_list),
        n_positions_removed, n_positions_total,
    )

    return filtered, stats


# ---------------------------------------------------------------------------
# Copy of pore diameter computation (uses filtered connections)
# ---------------------------------------------------------------------------

def compute_pore_dia_with_filter(results, atoms, spacing=0.5, spill_margin=2.0):
    """Compute LCD/PLD/DFSD with spill-filtered connections."""
    from mofography import mgr

    # Filter connections
    conn_orig = results["connections"]
    conn_filtered, stats = filter_spill_connections(
        conn_orig, results, atoms, spill_margin=spill_margin, spacing=spacing,
    )

    # Replace connections in results for downstream use
    results_filtered = dict(results)
    results_filtered["connections"] = conn_filtered
    # Also need to re-wrap as the pipeline expects
    results_filtered["connections"] = [
        np.array([c[0], c[1], c[2], c[3], c[4]], dtype=object)
        for c in conn_filtered
    ]

    # Build graph from filtered results
    graph = make_directed_graph_from_segmentation_results(results_filtered)
    base_region_labels = build_base_region_labels_from_results(results_filtered)

    dgrid = np.asarray(results["dgrid"])
    lcd = float(dgrid.max() * 2)

    plr = find_pore_limiting_radius_winding(
        results_filtered,
        graph=graph,
        base_region_labels=base_region_labels,
    )
    if plr is None:
        return lcd, 0.0, 0.0, dgrid, stats

    pld = float(plr * 2)

    # DFSD
    dfsd_threshold = _dfsd_threshold(plr, results_filtered, {})
    axis_clusters, _, _ = percolation_all_axes_winding(
        results_filtered,
        dfsd_threshold,
        graph=graph,
        base_region_labels=base_region_labels,
    )
    percolating_labels = set()
    for groups in axis_clusters.values():
        for labels in groups:
            percolating_labels.update(labels)
    percolating_labels = _expand_dfsd_labels_with_window_fingerprint(
        graph=graph,
        base_region_labels=base_region_labels,
        percolating_labels=percolating_labels,
        results=results_filtered,
        dgrid=dgrid,
        plr=float(plr),
        kwargs={},
        atoms=atoms,
    )

    dfsd = 0.0
    if percolating_labels:
        regions = np.asarray(results_filtered["regions"])
        mask = np.isin(regions, list(percolating_labels))
        if mask.any():
            dfsd = float(dgrid[mask].max() * 2)

    return lcd, pld, dfsd, dgrid, stats


# ---------------------------------------------------------------------------
# Main test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    structures = {
        "TOVKOG":  "examples/TestCIFs/TOVKOG_clean.cif",
        "VIGNOQ":  "examples/TestCIFs/VIGNOQ_clean.cif",
        "REHHOD":  "examples/TestCIFs/REHHOD_clean.cif",
        "IRMOF-1": "outputs/raspa3_md_benzene_eval/IRMOF-1.cif",
        "CuBTC":   "examples/TestCIFs/P1-Cu-BTC.cif",
        "NU-1000": "examples/TestCIFs/NU-1000.cif",
        "WUPNIF":  "examples/TestCIFs/WUPNIF_clean.cif",
    }

    print(f"{'Structure':>10s} | {'LCD':>8s} | {'PLD(orig)':>9s} | {'PLD(w2)':>8s} | {'PLD(filt)':>9s} | {'DFSD(orig)':>10s} | {'DFSD(filt)':>10s} | {'Conn-':>5s} | {'Pos-':>5s}")
    print("-" * 100)

    for name, path in structures.items():
        try:
            atoms = ase_read(path)
            radii = default_vdw_radii_angstrom(atoms)
            dgrid = dgrid_from_atoms_rust(atoms, spacing=0.5, radii=radii)

            res = do_periodic_segmentation_rust(
                dgrid, atoms, spacing=0.5, mask_thickness=0.05, h=0.2,
                watershed_backend="rust_det", return_dgrid=True,
                return_padded_regions=True,
            )

            # Original winding2 (with resegmentation guardrail)
            lcd2, pld2, dfsd2, _ = get_pore_dia_mgr_winding2(
                results=res, mat_atoms=atoms, auto_rerun=True
            )

            # Filtered (no resegmentation)
            lcd_f, pld_f, dfsd_f, _, stats = compute_pore_dia_with_filter(
                res, atoms, spacing=0.5, spill_margin=2.0,
            )

            conn_dropped = stats["connections_dropped"]
            pos_removed = stats["positions_removed"]

            print(
                f"{name:>10s} | {lcd_f:8.2f} | {0.0:9.2f} | {pld2:8.2f} | "
                f"{pld_f:9.2f} | {dfsd2:10.2f} | {dfsd_f:10.2f} | "
                f"{conn_dropped:5d} | {pos_removed:5d}"
            )

        except Exception as e:
            print(f"{name:>10s} | ERROR: {e}")

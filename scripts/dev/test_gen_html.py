#!/usr/bin/env python3
"""Generate two Plotly HTML visualizations of CuBTC watershed segmentation
with skimage backend at connectivity=2 and connectivity=3.

Shows: cell, regions surface, regions_pbc surface, maxima_pbc,
       pore graph (internal nodes + connections).
"""

import sys, time, os, numpy as np
import plotly.graph_objects as go

sys.path.insert(0, "/Users/arungopalan/repos/mofography")
import mofography as mgr
from mofography.visualize import (
    plot_surface_dgrid, plot_cell, scatter3d, plot_pore_graph, dark_layout,
)
from mofography.helper import find_coord_from_indices, compute_wall_cells, is_inside_or_close_to_boundary
from mofography import do_periodic_segmentation_rust
from ase.io import read
from ase.data import covalent_radii

import networkx as nx

# ── Parameters ──────────────────────────────────────────────────
CIF_PATH = "/Users/arungopalan/repos/mofography/mofography/examples/TestCIFs/P1-Cu-BTC.cif"
SPACING = 0.5
H = 0.2
MASK_THICKNESS = 0.05
THRESHOLD_ABS = 1.5
COMPACTNESS = 0
MIN_DISTANCE = 1.5
ISOVALUE = 2.0   # for surface rendering

# ── Load ────────────────────────────────────────────────────────
print("Loading CuBTC ...")
mat_atoms = read(CIF_PATH)
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
radii = covalent_radii[mat_updated.get_atomic_numbers()]

# ── Compute dgrid once ──────────────────────────────────────────
print("Computing dgrid ...")
t0 = time.time()
grid = mgr.dgrid_from_atoms(mat_atoms, radii=radii, spacing=SPACING, block_size=20000)
dgrid = grid.compute().astype(np.float32)
print(f"  shape={dgrid.shape}  range=[{dgrid.min():.3f}, {dgrid.max():.3f}]  time={time.time()-t0:.1f}s")

# ── Run for both connectivities ─────────────────────────────────
for conn in [2, 3]:
    label = f"skimage_conn{conn}"
    print(f"\n{'='*70}")
    print(f"Running do_periodic_segmentation_rust with {label} ...")
    t0 = time.time()

    # Monkey-patch: temporarily override connectivity in skimage backend
    # We'll call with watershed_backend="skimage" and then manually set connectivity
    # by post-patching the result.

    # Actually, let's run the full pipeline. But do_periodic_segmentation_rust
    # uses the hardcoded connectivity=2. We can't override it easily.
    # So let's compute watershed ourselves and then build the same outputs.

    # ── Alternative: manual computation ──────────────────────────
    from ase.geometry import cellpar_to_cell
    import mgr_rust

    cell = cellpar_to_cell(mat_atoms.cell.cellpar()).T

    # Maxima detection (same for both)
    merged, unmerged, raw = mgr_rust.detect_maxima_selective_smoothing_parallel(
        dgrid,
        preserve_threshold=3.0, smooth_threshold=6.0, max_smooth_radius=2.0,
        h=H, threshold_abs=THRESHOLD_ABS, min_distance=MIN_DISTANCE, cell=cell,
    )
    print(f"  Maxima: {len(merged)} merged")

    # Use proper padding from the pipeline
    from mofography.compute import calculate_padding_thickness
    padding_thickness = int(calculate_padding_thickness(dgrid, mat_atoms, SPACING)) + int(np.ceil(5 / SPACING))
    print(f"  Padding thickness: {padding_thickness}")

    # Build padded arrays
    markers = np.zeros_like(dgrid, dtype=np.int32)
    coords = tuple(merged.T)
    markers[coords] = np.arange(1, len(merged) + 1)

    void_mask = dgrid > MASK_THICKNESS
    padded_grid = np.pad(dgrid, padding_thickness, mode="wrap")
    padded_markers = np.pad(markers, padding_thickness, mode="wrap")
    padded_void = np.pad(void_mask, padding_thickness, mode="wrap")

    # Watershed with skimage at given connectivity
    from skimage.segmentation import watershed as sk_watershed

    padded_regions = sk_watershed(
        -padded_grid.astype(float),
        padded_markers.astype(np.int32),
        mask=padded_void.astype(bool),
        compactness=COMPACTNESS,
        connectivity=conn,
    ).astype(np.uint32)

    # Extract central region → regions_pbc
    regions_pbc = padded_regions[
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
    ]

    # Now compute unmerged regions (PBC copies of maxima)
    # Build maxima_info (periodic copies)
    original_shape = np.array(dgrid.shape)
    padded_shape = np.array(padded_grid.shape)
    maxima_info_list = []
    counter = 1
    for m in merged:
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
                maxima_info_list.append([
                    counter,
                    image_coord[0], image_coord[1], image_coord[2],
                    offset[0], offset[1], offset[2],
                ])
                counter += 1
    maxima_info = np.array(maxima_info_list, dtype=int)

    # Build unmerged markers and watershed
    unmerged_markers = np.zeros_like(padded_grid, dtype=int)
    for info in maxima_info:
        pos = info[1:4] + padding_thickness + info[4:7] * original_shape
        unmerged_markers[tuple(pos)] = info[0]

    unmerged_padded = sk_watershed(
        -padded_grid.astype(float),
        unmerged_markers.astype(np.int32),
        mask=padded_void.astype(bool),
        compactness=COMPACTNESS,
        connectivity=conn,
    ).astype(np.uint32)

    unmerged_regions = unmerged_padded[
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
        padding_thickness:-padding_thickness,
    ]

    # ── Compute connections using Rust window_connections ────────
    A_unit = cellpar_to_cell(mat_atoms.cell.cellpar()).T
    min_win_sep = SPACING * 3

    connections_list = mgr_rust.window_connections(
        unmerged_padded.astype(np.uint32),
        padded_grid.astype(np.float32),
        padding_thickness,
        A_unit,
        float(min_win_sep),
    )
    connections_list = [
        np.array([c[0], c[1], c[2], c[3], c[4]], dtype="object")
        for c in connections_list
    ]

    # Classify internal vs periodic
    wall_region_thicknesses = compute_wall_cells(
        mat_atoms, SPACING, wall_thickness=0.05
    )
    maxima_indices = maxima_info[:, 1:4] + maxima_info[:, 4:7] * np.array(dgrid.shape)
    inside_close_outside = is_inside_or_close_to_boundary(
        maxima_indices, dgrid.shape, wall_region_thicknesses
    )
    valid_maxima_mask = np.array(
        [loc in ["inside", "close_outside"] for loc in inside_close_outside]
    )
    valid_maxima_ids = np.where(valid_maxima_mask)[0] + 1

    for c in connections_list:
        if np.any(~np.isin(c[0], valid_maxima_ids)):
            c[2] = "periodic"

    n_int = sum(1 for c in connections_list if c[2] == "internal")
    n_per = sum(1 for c in connections_list if c[2] == "periodic")
    n_conn = sum(1 for c in connections_list if c[1])
    print(f"  Connections: {n_conn} total, {n_int} internal, {n_per} periodic")

    # ── Build pore graph (internal only) ─────────────────────────
    G = nx.Graph()
    G.add_nodes_from(range(1, len(merged) + 1))

    # Add node coordinates (Cartesian)
    for node in G.nodes():
        m_idx = node - 1
        coord = find_coord_from_indices(
            np.array([merged[m_idx]]), shape=dgrid.shape, ase_atoms=mat_atoms
        )[0]
        G.nodes[node]["node_coords"] = coord

    internal_edges = [
        tuple(c[0]) for c in connections_list
        if c[1] and c[2] == "internal"
    ]

    # Only add edges where both nodes are in the 1..len(merged) range
    filtered_edges = []
    for a, b in internal_edges:
        ai, bi = int(a), int(b)
        if 1 <= ai <= len(merged) and 1 <= bi <= len(merged):
            filtered_edges.append((ai, bi))
    G.add_edges_from(filtered_edges)

    # Also add labels attribute
    nx.set_node_attributes(
        G,
        values={i + 1: [i + 1] for i in range(len(G.nodes))},
        name="labels"
    )

    print(f"  Pore graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} internal edges")

    # ── Build Plotly figure ──────────────────────────────────────
    print("  Building figure ...")

    # Window scatter — extract per-window centers & sizes
    windows_pts = []
    window_labels = []
    for c in connections_list:
        if not c[1]:
            continue
        centers = c[3]  # (n_windows, 3) float array
        radii = c[4]    # (n_windows,) float array
        reg_pair = f"{int(c[0][0])}-{int(c[0][1])}"
        for wi in range(len(centers)):
            windows_pts.append(centers[wi])
            try:
                r = float(radii[wi])
            except (TypeError, IndexError):
                r = float(radii) if np.ndim(radii) == 0 else 0.0
            window_labels.append(f"{reg_pair}<br>{c[2]}<br>{r:.3f} \u212b")
    windows_pts = np.array(windows_pts) if windows_pts else np.zeros((0, 3))

    # Maxima in Cartesian
    maxima_cart = find_coord_from_indices(merged, shape=dgrid.shape, ase_atoms=mat_atoms)
    cage_radii = dgrid[tuple(merged.T)]
    maxima_hover = [f"{i+1} {r:.2f} \u212b" for i, r in enumerate(cage_radii)]

    traces = []

    # 1) Regions surface (no PBC)
    traces.append(plot_surface_dgrid(
        dgrid, isovalue=ISOVALUE, opacity=0.6, stepsize=1,
        region_labels=unmerged_regions, select="all",
        ase_atoms=mat_atoms, stylize=True, cscale=None, random_seed=123,
        name=f"Regions (no PBC) [{label}]",
    ))

    # 2) Regions PBC surface
    traces.append(plot_surface_dgrid(
        dgrid, isovalue=ISOVALUE, opacity=0.8, stepsize=1,
        region_labels=regions_pbc, select="all",
        ase_atoms=mat_atoms, stylize=True, cscale=None, random_seed=456,
        name=f"Regions PBC [{label}]",
    ))

    # 3) Cell
    traces.append(plot_cell(mat_atoms, color="white", name="cell"))

    # 4) Maxima PBC
    if len(maxima_cart) > 0:
        traces.append(scatter3d(
            maxima_cart,
            name="maxima PBC",
            marker_color=cage_radii,
            marker_size=15,
            hovertext=maxima_hover,
        ))

    # 5) Windows (internal only)
    internal_windows = []
    internal_wlabels = []
    for c in connections_list:
        if not c[1] or c[2] != "internal":
            continue
        centers = c[3]
        radii = c[4]
        reg_pair = f"{int(c[0][0])}-{int(c[0][1])}"
        for wi in range(len(centers)):
            internal_windows.append(centers[wi])
            try:
                r = float(radii[wi])
            except (TypeError, IndexError):
                r = float(radii) if np.ndim(radii) == 0 else 0.0
            internal_wlabels.append(f"{reg_pair}<br>{r:.3f} \u212b")

    if internal_windows:
        iw_pts = np.array(internal_windows)
        traces.append(scatter3d(
            iw_pts,
            name="windows (internal)",
            marker_color="yellowgreen",
            marker_size=8,
            opacity=0.6,
            hovertext=internal_wlabels,
        ))

    # 6) Pore graph (internal)
    pg_nodes_trace, pg_edges_trace = plot_pore_graph(G)
    pg_nodes_trace.name = f"pore graph nodes [{label}]"
    pg_edges_trace.name = f"pore graph edges [{label}]"
    traces.append(pg_nodes_trace)
    traces.append(pg_edges_trace)

    # ── Layout ───────────────────────────────────────────────────
    dark = dict(template="plotly_dark",
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)")

    fig = go.Figure(data=traces, layout=dark)
    fig.update_layout(
        title=f"CuBTC — skimage backend connectivity={conn} — "
              f"mask={MASK_THICKNESS}, h={H}, spacing={SPACING}<br>"
              f"{G.number_of_edges()} internal edges, {n_int} internal connections",
        showlegend=True,
    )

    out_path = f"/Users/arungopalan/repos/mofography/test_cubtc_skimage_conn{conn}.html"
    fig.write_html(out_path)
    print(f"  Wrote: {out_path}")
    print(f"  File size: {os.path.getsize(out_path) / 1024 / 1024:.1f} MB")

import os
print(f"\n{'='*70}")
print("Done! Both HTML files generated.")

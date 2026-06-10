#!/usr/bin/env python3
"""Generate two Plotly HTML visualizations of CuBTC watershed segmentation
comparing skimage backend at connectivity=2 (18-NN) vs connectivity=3 (26-NN).

Uses full do_periodic_segmentation_rust pipeline for correct pore graph construction.
Shows: cell, regions surface, regions_pbc surface, maxima_pbc,
       pore graph (internal nodes + edges), internal windows.
"""

import sys, os, time, numpy as np
import plotly.graph_objects as go

sys.path.insert(0, "/Users/arungopalan/repos/mofography")
import mofography as mgr
from mofography.visualize import (
    plot_surface_dgrid, plot_cell, scatter3d, plot_pore_graph,
)
from mofography.helper import find_coord_from_indices
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
ISOVALUE = 2.0
PRESERVE_THRESHOLD = 3.0
SMOOTH_THRESHOLD = 6.0
MAX_SMOOTH_RADIUS = 2.0
MIN_REGION_SIZE = 1e-6  # don't remove any small regions

# ── Load ────────────────────────────────────────────────────────
print("Loading CuBTC ...")
mat_atoms = read(CIF_PATH)

# ── Compute dgrid once ──────────────────────────────────────────
print("Computing dgrid ...")
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
radii = covalent_radii[mat_updated.get_atomic_numbers()]
t0 = time.time()
grid = mgr.dgrid_from_atoms(mat_updated, radii=radii, spacing=SPACING, block_size=20000)
dgrid = grid.compute().astype(np.float32)
print(f"  shape={dgrid.shape}  range=[{dgrid.min():.3f}, {dgrid.max():.3f}]  time={time.time()-t0:.1f}s")

# ── Run pipeline for both connectivities ────────────────────────
for conn in [2, 3]:
    label = f"skimage_conn{conn}"
    print(f"\n{'='*70}")
    print(f"Running do_periodic_segmentation_rust with {label} ...")
    t0 = time.time()

    seg = mgr.do_periodic_segmentation_rust(
        dgrid,
        mat_updated,
        mask_thickness=MASK_THICKNESS,
        h=H,
        min_distance=MIN_DISTANCE,
        threshold_abs=THRESHOLD_ABS,
        preserve_threshold=PRESERVE_THRESHOLD,
        smooth_threshold=SMOOTH_THRESHOLD,
        max_smooth_radius=MAX_SMOOTH_RADIUS,
        compactness=COMPACTNESS,
        watershed_backend="skimage",
        watershed_connectivity=conn,
        min_region_size=MIN_REGION_SIZE,
        compute_connections=True,
        return_dgrid=False,
    )
    elapsed = time.time() - t0
    print(f"  Pipeline time: {elapsed:.1f}s")

    regions_pbc = seg["regions_pbc"]
    regions = seg["regions"]
    maxima_pbc = seg["maxima_pbc"]
    maxima_info = seg["maxima_info"]
    connections = seg["connections"]
    pbc_groups = seg["pbc_groups"]

    n_unique = len(np.unique(regions_pbc)) - 1
    n_regions = len(np.unique(regions)) - 1

    n_int_conn = sum(1 for c in connections if c[1] and c[2] == "internal")
    n_per_conn = sum(1 for c in connections if c[1] and c[2] == "periodic")
    print(f"  regions_pbc unique: {n_unique}, regions unique: {n_regions}")
    print(f"  Internal connections: {n_int_conn}, Periodic: {n_per_conn}")

    # ── Build pore graph (internal only) ─────────────────────────
    internal_edges = [
        tuple(int(x) for x in c[0])
        for c in connections
        if c[1] and c[2] == "internal"
    ]

    G = nx.Graph()
    # Add all nodes that appear in internal edges
    all_internal_nodes = set()
    for a, b in internal_edges:
        all_internal_nodes.add(a)
        all_internal_nodes.add(b)
    G.add_nodes_from(sorted(all_internal_nodes))
    G.add_edges_from(internal_edges)
    nx.set_node_attributes(
        G,
        values={n: [n] for n in G.nodes()},
        name="labels",
    )

    # Add node coordinates (Cartesian) from maxima_info
    # maxima_info columns: [label, i, j, k, ox, oy, oz]
    for info in maxima_info:
        lbl = info[0]
        if lbl in G.nodes():
            idx = info[1:4] + info[4:7] * np.array(dgrid.shape)
            cart = find_coord_from_indices(
                np.array([idx]), shape=dgrid.shape, ase_atoms=mat_updated
            )[0]
            G.nodes[lbl]["node_coords"] = cart

    print(f"  Pore graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} internal edges")

    # ── Build Plotly figure ──────────────────────────────────────
    print("  Building figure ...")

    # Window scatter (internal only)
    interna_windows = []
    interna_wlabels = []
    for c in connections:
        if not c[1] or c[2] != "internal":
            continue
        centers = c[3]
        radii = c[4]
        reg_pair = f"{int(c[0][0])}-{int(c[0][1])}"
        for wi in range(len(centers)):
            interna_windows.append(centers[wi])
            try:
                r = float(radii[wi])
            except (TypeError, IndexError):
                r = float(radii) if np.ndim(radii) == 0 else 0.0
            interna_wlabels.append(f"{reg_pair}<br>{r:.3f} \u212b")

    cage_radii = dgrid[tuple(maxima_pbc.T)]
    maxima_hover = [f"{i+1} {r:.2f} \u212b" for i, r in enumerate(cage_radii)]

    name_suffix = f" (conn={conn})"
    traces = []

    # 1) Regions surface (no PBC)
    traces.append(plot_surface_dgrid(
        dgrid, isovalue=ISOVALUE, opacity=0.5, stepsize=1,
        region_labels=regions, select="all",
        ase_atoms=mat_updated, stylize=True, cscale=None, random_seed=123,
        name=f"Regions{name_suffix}",
    ))

    # 2) Regions PBC surface
    traces.append(plot_surface_dgrid(
        dgrid, isovalue=ISOVALUE, opacity=0.8, stepsize=1,
        region_labels=regions_pbc, select="all",
        ase_atoms=mat_updated, stylize=True, cscale=None, random_seed=456,
        name=f"Regions PBC{name_suffix}",
    ))

    # 3) Cell
    traces.append(plot_cell(mat_updated, color="white", name="cell"))

    # 4) Maxima PBC
    maxima_cart = find_coord_from_indices(
        maxima_pbc, shape=dgrid.shape, ase_atoms=mat_updated
    )
    if len(maxima_cart) > 0:
        traces.append(scatter3d(
            maxima_cart,
            name=f"maxima PBC{name_suffix}",
            marker_color=cage_radii,
            marker_size=12,
            hovertext=maxima_hover,
        ))

    # 5) Internal windows
    if interna_windows:
        iw_pts = np.array(interna_windows)
        traces.append(scatter3d(
            iw_pts,
            name=f"windows internal{name_suffix}",
            marker_color="yellowgreen",
            marker_size=6,
            opacity=0.5,
            hovertext=interna_wlabels,
        ))

    # 6) Pore graph
    pg_nodes, pg_edges = plot_pore_graph(G)
    pg_nodes.name = f"pore graph nodes{name_suffix}"
    pg_edges.name = f"pore graph edges{name_suffix}"
    traces.append(pg_nodes)
    traces.append(pg_edges)

    # ── Layout ───────────────────────────────────────────────────
    dark = dict(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    fig = go.Figure(data=traces, layout=dark)
    fig.update_layout(
        title=(
            f"CuBTC — skimage backend connectivity={conn} — "
            f"mask={MASK_THICKNESS}, h={H}, spacing={SPACING}<br>"
            f"{G.number_of_edges()} internal edges, "
            f"{n_int_conn} internal / {n_per_conn} periodic connections"
        ),
        showlegend=True,
    )

    out_path = f"/Users/arungopalan/repos/mofography/test_cubtc_skimage_conn{conn}_v2.html"
    fig.write_html(out_path)
    mb = os.path.getsize(out_path) / (1024 * 1024)
    print(f"  Wrote: {out_path} ({mb:.1f} MB)")

print(f"\n{'='*70}")
print("Done! Both HTML files generated.")
print("  conn2: test_cubtc_skimage_conn2_v2.html")
print("  conn3: test_cubtc_skimage_conn3_v2.html")

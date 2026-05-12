#!/usr/bin/env python3
"""Generate CuBTC static HTML + LTA96 movie HTML using Figure3D API.
v3: unique per-region colors, label-0 fix, Laplacian mesh smoothing."""

import sys, time, numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mofography as mgr
from ase.io import read
from skimage.measure import regionprops


def make_surface_trace_smooth(
    atoms, dgrid, labels_array, isovalue=1.0, step_size=1,
    trace_id="surface", name="Surface",
    hover_prefix="Region",
    opacity=1.0, roughness=0.6, metalness=0.0,
    smooth_iters=5,
):
    """Create a surface trace with unique per-region colors, label-0 fix, and Laplacian smoothing."""
    import trimesh
    from mofography.movies import _make_trace as _mt

    # --- Step 1: raw marching cubes ---
    mesh_data = mgr.compute_surface_mesh(atoms, dgrid, labels_array, isovalue=isovalue, step_size=step_size)
    verts = np.array(mesh_data["positions"], dtype=float).reshape(-1, 3)
    faces = np.array(mesh_data["indices"], dtype=int).reshape(-1, 3)
    labels = np.array(mesh_data["labels"], dtype=int)
    n_v = len(verts)

    # --- Step 2: Laplacian smoothing ---
    if smooth_iters > 0 and n_v > 0:
        mesh_s = trimesh.Trimesh(vertices=verts.copy(), faces=faces.copy(), process=False)
        try:
            from trimesh.smoothing import filter_laplacian
            mesh_s = filter_laplacian(mesh_s, iterations=smooth_iters, lamb=0.5)
            verts = np.asarray(mesh_s.vertices, dtype=float)
        except Exception:
            pass  # smoothing is best-effort

    # --- Step 3: fix label-0 vertices (find nearest non-zero label) ---
    zero_mask = labels == 0
    if zero_mask.any():
        nonzero_idx = np.where(~zero_mask)[0]
        zero_idx = np.where(zero_mask)[0]
        from scipy.spatial import cKDTree
        tree = cKDTree(verts[nonzero_idx])
        _, nn = tree.query(verts[zero_idx], k=1)
        labels[zero_idx] = labels[nonzero_idx[nn]]

    # --- Step 4: unique per-region colors using HSV wheel ---
    uniq = sorted(int(x) for x in np.unique(labels) if x > 0)
    n_regions = len(uniq)

    # Generate n_regions distinct hues evenly spaced
    hues = np.linspace(0, 1, n_regions, endpoint=False)
    import colorsys
    palette = {}
    for i, rid in enumerate(uniq):
        # Use high saturation, medium-high value for vibrant distinct colors
        rgb = colorsys.hsv_to_rgb(hues[i], 0.7, 0.85)
        palette[rid] = list(rgb)

    # Fallback gray for label 0 (shouldn't happen after fix)
    colors = np.array([palette.get(int(l), [0.55, 0.55, 0.55]) for l in labels], dtype=float)

    # --- Step 5: hover text ---
    hover = [f"{hover_prefix}: {l}" if l > 0 else "Framework" for l in labels]

    # --- Step 6: trace ---
    data = {
        "positions": verts.reshape(-1).tolist(),
        "indices": faces.reshape(-1).tolist(),
        "colors": colors.reshape(-1).tolist(),
        "labels": labels.tolist(),
        "n_vertices": int(n_v),
        "n_triangles": int(len(faces)),
    }
    return _mt(
        trace_id=trace_id, trace_type="mesh", name=name,
        data=data,
        style={"opacity": opacity, "roughness": roughness, "metalness": metalness, "side": "DoubleSide"},
        hover=hover, visible=True,
    )


# =============================================================================
# PART 1: CuBTC static
# =============================================================================
print("=" * 60)
print("PART 1: CuBTC STATIC SCENE")
print("=" * 60)

CIF = "mofography/examples/TestCIFs/P1-Cu-BTC.cif"
SPACING = 0.5
ISOVAL = 1.0
MASK_THICKNESS = 0.5

print("\nLoading CuBTC...")
atoms_raw = read(CIF)
atoms = mgr.update_with_covalent_bonds(atoms_raw)
n_bonds = len(atoms.info.get("covalent_bonds", []))
print(f"  {atoms.get_chemical_formula()}, {len(atoms)} atoms, {n_bonds} bonds (covalent)")

print("Computing distance grid...")
t0 = time.time()
dgrid = mgr.dgrid_from_atoms_cpu(atoms, spacing=SPACING, block_size=20000).astype(np.float32)
print(f"  shape={dgrid.shape} ({time.time()-t0:.1f}s)")

print("Running periodic segmentation...")
t0 = time.time()
seg = mgr.do_periodic_segmentation_rust(
    dgrid, atoms, spacing=SPACING, mask_thickness=MASK_THICKNESS, h=0.2,
    threshold_abs=1.5, min_distance=1.5,
    compute_connections=True, return_dgrid=False, extra_padding=5,
)
elapsed = time.time() - t0
regions = seg["regions"]
regions_pbc = seg["regions_pbc"]
connections = seg["connections"]
all_maxima = seg["all_maxima"]
unique_regions = seg["unique_regions"]
pbc_groups = seg["pbc_groups"]

n_unique_pbc = len(np.unique(regions_pbc)) - 1
n_unique_raw = len(np.unique(regions)) - 1
n_int = sum(1 for c in connections if c[1] and c[2] == "internal")
n_per = sum(1 for c in connections if c[1] and c[2] == "periodic")
print(f"  {elapsed:.1f}s | regions_pbc={n_unique_pbc}, regions_raw={n_unique_raw}")
print(f"  connections: {n_int} int, {n_per} per, {len(connections)} total")

# Clustering
connections_trimmed = np.array([c for c in connections if np.any(np.isin(c[0], unique_regions))], dtype=object)
print(f"  trimmed: {len(connections_trimmed)}")

print("\nComputing pocket-type clustering...")
t0 = time.time()
maxima_in_cell = all_maxima[:regions.max()]
maxima_wrapped = np.mod(maxima_in_cell, np.asarray(dgrid.shape, dtype=int)).astype(int)
maxima_radii = dgrid[tuple(maxima_wrapped.T)]

n_max = len(maxima_in_cell)
conn_filtered = [c for c in connections_trimmed if 1 <= int(c[0][0]) <= n_max and 1 <= int(c[0][1]) <= n_max]
print(f"  connections in range: {len(conn_filtered)}")

rag = mgr.rag_from_connections(conn_filtered, maxima_in_cell, only_use_internal=False)
rag = mgr.add_maxima_to_rag(rag, maxima_in_cell, maxima_radii, dgrid.shape, atoms)
props = regionprops(regions)
rag = mgr.add_pbc_to_rag(rag, pbc_groups)
rag = mgr.add_centroids_to_rag_with_pbc(rag, props, pbc_groups, atoms, regions.shape)
rag = mgr.add_psd_to_rag_with_pbc(rag, props, pbc_groups, bins=50, ase_atoms=atoms, shape=regions.shape)
rag = mgr.add_clustering_to_rag_geometry_js(rag, props, pbc_groups, atoms, dgrid.shape, use_pbc=True, js_threshold=0.15)
regions_cluster = mgr.clustering_grid_from_regions(regions, rag, cluster_key='cluster_geom_pbc')

n_clusters = len(np.unique(regions_cluster)) - 1
print(f"  {time.time()-t0:.1f}s | pocket types: {n_clusters}")

# --- Build Figure3D with smoothed, uniquely-colored surfaces ---
print("\nBuilding Figure3D...")

fig = (
    mgr.Figure3D()
    .set_camera(eye=(1.5, 0.8, 1.2), center=(0, 0, 0), up=(0, 0, 1))
    .add_framework(atoms)
    .add_cell(atoms)
)

# Three surfaces — each with unique per-region colors + hover + smoothing
fig.add_trace(make_surface_trace_smooth(
    atoms, dgrid, regions, isovalue=ISOVAL, step_size=1,
    trace_id="surface_regions_raw", name="Regions (raw)", hover_prefix="Region",
    smooth_iters=5,
))
fig.add_trace(make_surface_trace_smooth(
    atoms, dgrid, regions_pbc, isovalue=ISOVAL, step_size=1,
    trace_id="surface_regions_pbc", name="Regions (PBC merged)", hover_prefix="PBC group",
    smooth_iters=5,
))
fig.add_trace(make_surface_trace_smooth(
    atoms, dgrid, regions_cluster, isovalue=ISOVAL, step_size=1,
    trace_id="surface_regions_cluster", name="Pocket types (cluster)", hover_prefix="Cluster",
    smooth_iters=5,
))

# Pore graph — ALL visible
fig.add_maxima(seg, atoms, hide_periodic=False)
fig.add_connections(seg, atoms, hide_periodic=False)
fig.add_windows(seg, hide_periodic=False)

fig.groups.update({
    "surface_regions_raw":    {"traces": ["surface_regions_raw"],    "label": "Regions (raw)"},
    "surface_regions_pbc":    {"traces": ["surface_regions_pbc"],    "label": "Regions (PBC merged)"},
    "surface_regions_cluster": {"traces": ["surface_regions_cluster"], "label": "Pocket types (cluster)"},
    "internal_pore_graph":   {"traces": ["maxima_internal", "windows_internal", "connections_internal", "connections_internal_nodes"], "label": "Internal Pore Graph"},
    "periodic_pore_graph":   {"traces": ["maxima_outside", "windows_periodic", "connections_periodic", "connections_periodic_nodes"], "label": "Periodic Pore Graph"},
})

fig.meta["description"] = "CuBTC — smoothed surfaces, unique per-region colors, hover labels"
fig.meta["n_regions_raw"] = n_unique_raw
fig.meta["n_regions_pbc"] = n_unique_pbc
fig.meta["n_pocket_types"] = n_clusters

out_static = "cubtc_static.html"
fig.to_html(out_static)
mb = Path(out_static).stat().st_size / (1024 * 1024)
print(f"\nSaved: {out_static} ({mb:.1f} MB)")

# =============================================================================
# PART 2: LTA96 movie
# =============================================================================
print("\n" + "=" * 60)
print("PART 2: LTA96 MOVIE")
print("=" * 60)

print("\nLoading LTA96 movie...")
snapshots = mgr.read_raspa_pdb_movie(
    "mofography/examples/Movie_LTA96_1.1.1_298.150000_0.000000_component_K_0.pdb",
    max_frames=50,
)
framework_raw = read("mofography/examples/TestCIFs/LTA96.cif")
framework = mgr.update_with_covalent_bonds(framework_raw)
n_frames = len(snapshots)
n_bonds_lta = len(framework.info.get("covalent_bonds", []))
print(f"  {n_frames} frames, {len(framework)} atoms, {n_bonds_lta} bonds (covalent)")

print("Computing segmentation...")
t0 = time.time()
seg_lta = mgr.compute_segmentation_for_snapshot({"atoms": framework}, spacing=0.5)
print(f"  {time.time()-t0:.1f}s")

from mofography.movies import _make_trace as _mt
from ase.data import covalent_radii

PINK = [1.0, 0.41, 0.71]

print("Building movie figure...")
movie_fig = mgr.Figure3D()
movie_fig.set_camera(eye=(2.0, 0.0, 1.35), center=(0, 0, 0), up=(0, 0, 1))
movie_fig.update_layout({"fps": 20.0, "n_frames": n_frames})

# Static traces — explicit material settings matching static figure
movie_fig.add_framework(framework, opacity=1.0, roughness=0.38, metalness=0.0, side="DoubleSide")
movie_fig.add_cell(framework)

# Smoothed, uniquely-colored surface
movie_fig.add_trace(make_surface_trace_smooth(
    framework, seg_lta["dgrid"], seg_lta["regions"],
    isovalue=1.0, step_size=1,
    trace_id="surface", name="Surface", hover_prefix="Region",
    smooth_iters=5,
))

maxima_lta = mgr.make_maxima_trace(seg_lta, framework)
movie_fig.add_trace(maxima_lta["internal"])
movie_fig.add_trace(maxima_lta["outside"])
windows_lta = mgr.make_windows_trace(seg_lta)
movie_fig.add_trace(windows_lta["internal"])
movie_fig.add_trace(windows_lta["periodic"])
conns_lta = mgr.make_connections_trace(seg_lta, framework)
movie_fig.add_trace(conns_lta["internal_lines"])
movie_fig.add_trace(conns_lta["internal_nodes"])
movie_fig.add_trace(conns_lta["periodic_lines"])
movie_fig.add_trace(conns_lta["periodic_nodes"])

movie_fig.groups = {
    "internal_pore_graph": {"traces": ["maxima_internal", "windows_internal", "connections_internal", "connections_internal_nodes"], "label": "Internal Pore Graph"},
    "periodic_pore_graph": {"traces": ["maxima_outside", "windows_periodic", "connections_periodic", "connections_periodic_nodes"], "label": "Periodic Pore Graph"},
    "guests": {"traces": ["guests"], "label": "K+ Adsorbates"},
}

for i, snap in enumerate(snapshots):
    a = snap["atoms"]
    n = int(len(a))
    r = covalent_radii[a.get_atomic_numbers()] * 0.8
    movie_fig.add_trace(_mt(
        trace_id="guests", trace_type="instanced_spheres", name="K+ Adsorbates",
        data={
            "positions": a.get_positions().reshape(-1).tolist(),
            "symbols": a.get_chemical_symbols(), "colors": [],
            "radii": r.astype(float).tolist(), "n_atoms": n,
        },
        style={"color": PINK, "opacity": 1.0, "shininess": 10, "specular": 0x333333},
        hover=None, visible=True,
    ), frame=i)

movie_fig.meta["description"] = "LTA96 — K+ movie with smoothed surface + unique region colors"

out_movie = "lta96_k_movie_v2.html"
movie_fig.to_html(out_movie)
mb = Path(out_movie).stat().st_size / (1024 * 1024)
print(f"\nSaved: {out_movie} ({mb:.1f} MB)")

print("\n" + "=" * 60)
print("DONE")
print(f"  Static: {out_static}")
print(f"  Movie:  {out_movie}")
print("=" * 60)

"""
Movie visualization module for Mofography.

Provides a Plotly-like API for building interactive Three.js HTML viewers
of framework + segmentation + MD snapshot animations.

Quick-start — fixed framework (framework + pores static, adsorbates animate)
--------------------------------------------------------------------------

```python
import mofography as mgr
from pathlib import Path

# 1. Read a RASPA PDB movie into canonical snapshots
snapshots = mgr.read_raspa_pdb_movie(
    "mofography/examples/Movie_LTA96_1.1.1_298.150000_0.000000_component_K_0.pdb",
    max_frames=200,
)

# 2. Load the framework CIF
framework = mgr.read("examples/TestCIFs/LTA96.cif")

# 3. Build the movie figure — segmentation computed once, guests animate
fig = mgr.build_fixed_framework_movie(
    snapshots=snapshots,
    framework_atoms=framework,
    spacing=0.25,          # dgrid spacing
    isovalue=1.0,          # surface isosurface level
    surf_step_size=2,      # marching-cubes decimation
    fps=20.0,
    eye=(2.0, 0.0, 1.35),
    center=(0.0, 0.0, 0.0),
    up=(0.0, 0.0, 1.0),
)

# 4. Export interactive HTML
fig.to_html("my_movie.html")
```

Or use the one-liner convenience wrapper:

```python
mgr.make_movie_from_raspa_pdb(
    pdb_path="movie.pdb",
    cif_path="framework.cif",
    out_html="output.html",
    max_frames=200,
    spacing=0.25,
)
```

The exported HTML is self-contained (single file) and works offline once
cached — no Python runtime required to view.

Supported input formats (v1)
----------------------------
- RASPA PDB movies → ``read_raspa_pdb_movie(...)`` or ``read_movie_snapshots(..., movie_format="raspa_pdb")``
- Any ASE-readable trajectory (XYZ, extxyz, DCD, TRR, XTC, LAMMPS dump, …) →
  ``read_movie_snapshots(path, movie_format="ase")`` or ``movie_format="auto"``
  (non-``.pdb`` paths default to ASE indexed reads, one frame at a time for memory).

All adapters produce the same canonical ``list[snapshot]`` dict so downstream
code (``build_fixed_framework_movie``, ``compute_segmentation_for_snapshot``)
works identically regardless of source format.

Use ``frame_stride`` and ``position_decimals`` on ``read_movie_snapshots`` to shrink
embedded JSON and speed up browser playback for long trajectories.

Flexible framework (recompute pores per frame)
----------------------------------------------
Planned for Phase 2.  You will pass ``recompute_indices`` to control which
snapshots trigger segmentation recomputation (e.g. every 100th frame, or an
explicit list).  Uncomputed frames reuse the last known surface / graph.

Parallel / high-throughput workflows
------------------------------------
The module is intentionally split into pure-computation helpers
(``compute_segmentation_for_snapshot``) and the ``MovieFigure`` builder so
that the expensive steps can be wrapped in ``dask.delayed`` or a process pool:

```python
from dask import delayed, compute

tasks = [delayed(mgr.compute_segmentation_for_snapshot)(snap)
         for snap in snapshots]
results = compute(*tasks)   # parallel via Rust + Dask

# Then assemble traces into a single MovieFigure
fig = mgr.MovieFigure()
for seg, snap in zip(results, snapshots):
    fig.add_trace(mgr.make_surface_trace(seg["dgrid"], seg["regions"], snap["atoms"]))
```

Sidecar / streaming mode (large trajectories)
---------------------------------------------
For movies with thousands of frames, ``fig.to_html(..., sidecar=True)``
writes frame data to a ``.frames.jsonl`` file next to the HTML and the
JavaScript viewer fetches chunks on demand rather than embedding everything.
"""

from __future__ import annotations

import json
import os
import warnings
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import trimesh
from ase import Atoms
from ase.data import atomic_masses as atomic_masses_table
from ase.data import covalent_radii
from ase.geometry import cellpar_to_cell
from ase.io import read
from matplotlib import colors as mcolors
from skimage.measure import marching_cubes

from .helper import get_fractional_coordinates, interpolate_labels

# ---------------------------------------------------------------------------
# Default styles registry (shared between Plotly and Three.js layers)
# ---------------------------------------------------------------------------

DEFAULT_STYLES: dict[str, dict[str, Any]] = {
    "framework": {
        "color": [0.7, 0.7, 0.7],
        "opacity": 0.95,
        "roughness": 0.38,
        "metalness": 0.0,
        "side": "DoubleSide",
    },
    "guests": {
        "color": [1.0, 0.41, 0.71],
        "opacity": 1.0,
        "shininess": 10,
        "specular": 0x151515,
    },
    "surface": {
        "opacity": 0.92,
        "roughness": 0.34,
        "metalness": 0.0,
        "side": "FrontSide",
    },
    "cell": {
        "color": 0xFFA500,
        "linewidth": 2,
    },
    "internal_maxima": {
        "color": [0.85, 0.85, 0.85],
        "size": 5.0,
    },
    "periodic_maxima": {
        "color": [0.85, 0.85, 0.85],
        "size": 5.0,
    },
    "internal_windows": {
        "color": 0xFFE600,
        "size": 2.5,
    },
    "periodic_windows": {
        "color": 0xFFE600,
        "size": 2.5,
    },
    "internal_connections": {
        "color": 0xFF0000,
        "linewidth": 16,
    },
    "periodic_connections": {
        "color": 0xFF0000,
        "dashSize": 0.9,
        "gapSize": 0.5,
        "linewidth": 16,
    },
    "internal_nodes": {
        "color": 0x66CCFF,
        "size": 1.10,
    },
    "periodic_nodes": {
        "color": 0xFF9966,
        "size": 1.10,
    },
}

# CPK colors for hover text generation
CPK_COLORS = {
    "H": "#FFFFFF",
    "C": "#909090",
    "N": "#3050F8",
    "O": "#FF0D0D",
    "S": "#FFFF30",
    "P": "#FF8000",
    "F": "#90E050",
    "CL": "#1FF01F",
    "BR": "#A62929",
    "I": "#940094",
    "NA": "#AB5CF2",
    "K": "#8F40D4",
    "CA": "#3DFF00",
}


def _hex_to_rgb01(hex_color: str) -> list[float]:
    c = hex_color.lstrip("#")
    if len(c) != 6:
        return [0.8, 0.8, 0.8]
    return [int(c[i : i + 2], 16) / 255.0 for i in (0, 2, 4)]


def _resolve_style(style_ref: str | dict | None) -> dict[str, Any]:
    if isinstance(style_ref, dict):
        return style_ref
    if style_ref is None:
        return {}
    return DEFAULT_STYLES.get(style_ref, {})


# ---------------------------------------------------------------------------
# Scene model helpers
# ---------------------------------------------------------------------------

def _make_trace(
    trace_id: str,
    trace_type: str,
    name: str,
    data: dict[str, Any],
    style: str | dict | None = None,
    hover: list[str] | None = None,
    visible: bool = True,
) -> dict[str, Any]:
    return {
        "id": trace_id,
        "type": trace_type,
        "name": name,
        "visible": visible,
        "style": _resolve_style(style),
        "data": data,
        "hover": hover,
    }


def _serialize_trimesh(mesh: trimesh.Trimesh) -> dict[str, Any]:
    v = np.asarray(mesh.vertices, dtype=float)
    f = np.asarray(mesh.faces, dtype=int)
    colors = None
    if hasattr(mesh, "visual") and hasattr(mesh.visual, "vertex_colors"):
        vc = np.asarray(mesh.visual.vertex_colors)
        if vc.ndim == 2 and vc.shape[0] == len(v):
            if vc.shape[1] >= 3:
                colors = vc[:, :3].astype(float) / 255.0
    if colors is None:
        colors = np.tile(np.asarray([[0.7, 0.7, 0.7]], dtype=float), (len(v), 1))
    return {
        "positions": v.reshape(-1).tolist(),
        "indices": f.reshape(-1).tolist(),
        "colors": colors.reshape(-1).tolist(),
        "n_vertices": int(len(v)),
        "n_triangles": int(len(f)),
    }


def _framework_vertex_hovertext(
    mesh: trimesh.Trimesh,
    ase_atoms: Atoms,
    bond_rgba: tuple[int, int, int, int] = (127, 128, 129, 255),
) -> list[str]:
    bonds = list(ase_atoms.info.get("covalent_bonds", []))
    n_atoms = int(ase_atoms.get_global_number_of_atoms())
    n_bonds = len(bonds)
    vc = np.asarray(mesh.visual.vertex_colors)
    if vc.size == 0 or vc.shape[0] != len(mesh.vertices):
        return [""] * len(mesh.vertices)
    if vc.shape[1] < 4:
        pad = np.full((vc.shape[0], 4 - vc.shape[1]), 255, dtype=vc.dtype)
        vc = np.hstack([vc, pad])
    vc_u8 = np.clip(
        np.rint(vc.astype(float) * 255.0) if vc.max() <= 1.0 else vc, 0, 255
    ).astype(np.uint8)
    colors = ["rgba(" + ",".join(c.astype(str)) + ")" for c in vc_u8]
    vertexcolors = np.asarray(colors, dtype=object)
    bond_rgba_string = "rgba(" + ",".join(str(x) for x in bond_rgba) + ")"
    n_verts = len(mesh.vertices)
    sphere_mask = vertexcolors != bond_rgba_string
    bond_mask = vertexcolors == bond_rgba_string
    n_sphere_verts = int(np.sum(sphere_mask))
    n_bond_verts = int(np.sum(bond_mask))
    if n_atoms <= 0 or n_sphere_verts % n_atoms != 0:
        return [""] * n_verts
    vertices_per_sphere = n_sphere_verts // n_atoms
    if n_bonds <= 0:
        vertices_per_bond = 0
    elif n_bond_verts % n_bonds != 0:
        return [""] * n_verts
    else:
        vertices_per_bond = n_bond_verts // n_bonds

    atomic_numbers = ase_atoms.get_atomic_numbers()
    chemical_symbols = ase_atoms.get_chemical_symbols()
    masses = atomic_masses_table[atomic_numbers]
    positions = np.round(ase_atoms.get_positions(), decimals=4)
    pos_str = [" ".join(map(str, row)) for row in positions]

    vertex_atom_ids = np.repeat(np.arange(n_atoms), vertices_per_sphere)
    vertex_symbols = np.repeat(chemical_symbols, vertices_per_sphere)
    vertex_masses = np.repeat(masses, vertices_per_sphere)
    vertex_pos = np.repeat(pos_str, vertices_per_sphere)

    hover_atoms = [
        f"Id: {i}\nElement: {s}\nMass: {m}\nPosition: {p}"
        for i, s, m, p in zip(
            vertex_atom_ids, vertex_symbols, vertex_masses, vertex_pos, strict=False
        )
    ]
    bond_atom_ids = np.repeat(np.asarray(bonds, dtype=int), vertices_per_bond, axis=0)
    hover_bonds = [
        "Bond between atoms " + " and ".join(v.astype(str)) for v in bond_atom_ids
    ]
    return hover_atoms + hover_bonds


def _serialize_framework_mesh(atoms: Atoms) -> dict[str, Any]:
    from .visualize import mesh_for_atoms

    mesh = mesh_for_atoms(
        atoms,
        kwargs_spheres=dict(scale=0.8),
        kwargs_cylinders=dict(radius=0.25),
    )
    payload = _serialize_trimesh(mesh)
    payload["vertex_hover"] = _framework_vertex_hovertext(mesh, atoms)
    if len(payload["vertex_hover"]) != payload["n_vertices"]:
        payload["vertex_hover"] = [""] * int(payload["n_vertices"])
    payload["n_atoms"] = int(len(atoms))
    return payload


def _serialize_cell_lines(atoms: Atoms) -> dict[str, Any]:
    A = cellpar_to_cell(atoms.cell.cellpar()).T
    a, b, c = A[:, 0], A[:, 1], A[:, 2]
    o = np.array([0.0, 0.0, 0.0])
    corners = [o, a, b, c, a + b, a + c, b + c, a + b + c]
    edges = [
        (0, 1), (0, 2), (0, 3),
        (1, 4), (2, 4),
        (2, 6), (3, 6),
        (3, 5), (1, 5),
        (4, 7), (5, 7), (6, 7),
    ]
    pts = []
    for i, j in edges:
        pts.extend([corners[i], corners[j]])
    arr = np.asarray(pts, dtype=float)
    return {"positions": arr.reshape(-1).tolist(), "n_segments": len(edges)}


def _serialize_segmentation_surface_mesh(
    atoms: Atoms,
    seg: dict[str, Any],
    isovalue: float,
    step_size: int,
) -> dict[str, Any]:
    dgrid = seg["dgrid"]
    regions = seg["regions"]
    verts, faces, _normals, _vals = marching_cubes(
        dgrid, level=isovalue, step_size=max(step_size, 1)
    )
    if len(faces) == 0:
        return {
            "positions": [],
            "indices": [],
            "colors": [],
            "labels": [],
            "n_vertices": 0,
            "n_triangles": 0,
        }
    cell = atoms.get_cell().array.T
    shape = np.asarray(dgrid.shape, dtype=float)
    verts_cart = np.dot(cell, (verts / shape).T).T
    mesh = trimesh.Trimesh(vertices=verts_cart, faces=faces, process=False)

    frac = get_fractional_coordinates(np.asarray(mesh.vertices, dtype=float), atoms)
    interp = interpolate_labels(regions)
    labels = interp(frac % 1).astype(int)
    cscale = [
        [102, 197, 204],
        [143, 200, 178],
        [184, 203, 152],
        [225, 206, 126],
        [246, 200, 113],
        [247, 171, 115],
        [240, 162, 152],
        [224, 173, 224],
    ]
    uniq = sorted([int(x) for x in np.unique(labels) if int(x) > 0])
    cmap = {rid: cscale[i % len(cscale)] for i, rid in enumerate(uniq)}
    colors = np.asarray(
        [cmap.get(int(r), [180, 180, 180]) for r in labels], dtype=float
    ) / 255.0
    return {
        "positions": np.asarray(mesh.vertices, dtype=float).reshape(-1).tolist(),
        "indices": np.asarray(mesh.faces, dtype=int).reshape(-1).tolist(),
        "colors": colors.reshape(-1).tolist(),
        "labels": labels.tolist(),
        "n_vertices": int(len(mesh.vertices)),
        "n_triangles": int(len(mesh.faces)),
    }


def _serialize_maxima_trace(seg: dict[str, Any], atoms: Atoms) -> dict[str, Any]:
    dgrid = np.asarray(seg["dgrid"])
    regions = np.asarray(seg.get("regions", np.zeros(dgrid.shape, dtype=int)))
    maxima = np.asarray(seg["all_maxima"], dtype=float)
    if maxima.size == 0:
        return {
            "internal_positions": [],
            "internal_sizes": [],
            "internal_colors": [],
            "internal_hover": [],
            "internal_n_points": 0,
            "outside_positions": [],
            "outside_sizes": [],
            "outside_colors": [],
            "outside_hover": [],
            "outside_n_points": 0,
            "n_points": 0,
        }
    maxima_wrapped = np.mod(maxima, np.asarray(dgrid.shape, dtype=int)).astype(int)
    cage_radii = dgrid[tuple(maxima_wrapped.T)]
    cell = atoms.get_cell().array.T
    frac = maxima / np.asarray(dgrid.shape, dtype=float)
    pos = np.dot(cell, frac.T).T

    # Segmentation palette (same 8 colors used for surface regions)
    cscale = [
        [102, 197, 204],
        [143, 200, 178],
        [184, 203, 152],
        [225, 206, 126],
        [246, 200, 113],
        [247, 171, 115],
        [240, 162, 152],
        [224, 173, 224],
    ]

    # Color maxima by their region label for consistency with the surface mesh.
    # Each maxima is at a grid point; we look up its region in the regions array.
    region_labels = regions[tuple(maxima_wrapped.T)]
    uniq_regions = sorted([int(x) for x in np.unique(region_labels) if int(x) > 0])
    cmap = {rid: cscale[i % len(cscale)] for i, rid in enumerate(uniq_regions)}
    colors = [list(np.asarray(cmap.get(int(r), [180, 180, 180]), dtype=float) / 255.0) for r in region_labels]

    # Marker scale (matches devtools/md_pdb_threejs_viewer: 6 + 5 * norm from cage radii)
    vals = np.asarray(cage_radii, dtype=float)
    vmin = float(vals.min()) if len(vals) else 0.0
    vmax = float(vals.max()) if len(vals) else 1.0
    denom = max(vmax - vmin, 1e-12)
    norm = (vals - vmin) / denom
    sizes = (6.0 + 5.0 * norm).astype(float).tolist()
    hover = []
    for i, (g, p, r) in enumerate(zip(maxima, pos, cage_radii, strict=False)):
        hover.append(
            f"trace: local maxima\n"
            f"index: {i + 1}\n"
            f"grid: [{int(g[0])}, {int(g[1])}, {int(g[2])}]\n"
            f"xyz: ({float(p[0]):.3f}, {float(p[1]):.3f}, {float(p[2]):.3f})\n"
            f"size: {float(r):.3f} A"
        )

    connections = list(seg.get("connections", []))
    periodic_nodes: set[int] = set()
    internal_nodes: set[int] = set()
    for c in connections:
        if len(c) < 3 or not c[1]:
            continue
        pair = {int(c[0][0]), int(c[0][1])}
        if str(c[2]) == "periodic":
            periodic_nodes.update(pair)
        elif str(c[2]) == "internal":
            internal_nodes.update(pair)
    outside_nodes = periodic_nodes - internal_nodes
    inside_nodes = {i + 1 for i in range(len(maxima)) if (i + 1) not in outside_nodes}

    def _pack(ids: set[int]):
        if not ids:
            return [], [], [], [], 0
        idx = np.asarray([i - 1 for i in ids], dtype=int)
        p = pos[idx]
        s = np.asarray(sizes, dtype=float)[idx]
        c = np.asarray(colors, dtype=float)[idx]
        h = [hover[i] for i in idx]
        return p.reshape(-1).tolist(), s.tolist(), c.reshape(-1).tolist(), h, int(len(ids))

    in_pos, in_sizes, in_cols, in_hover, in_n = _pack(inside_nodes)
    out_pos, out_sizes, out_cols, out_hover, out_n = _pack(outside_nodes)
    return {
        "internal_positions": in_pos,
        "internal_sizes": in_sizes,
        "internal_colors": in_cols,
        "internal_hover": in_hover,
        "internal_n_points": in_n,
        "outside_positions": out_pos,
        "outside_sizes": out_sizes,
        "outside_colors": out_cols,
        "outside_hover": out_hover,
        "outside_n_points": out_n,
        "n_points": int(len(pos)),
    }


def _serialize_windows_trace(seg: dict[str, Any]) -> dict[str, Any]:
    conns = list(seg.get("connections", []))
    internal_windows = []
    internal_hover = []
    periodic_windows = []
    periodic_hover = []
    for c in conns:
        if len(c) < 5 or not c[1]:
            continue
        ctype = str(c[2])
        pts = c[3] if c[3] is not None else []
        wr = c[4] if c[4] is not None else []
        n1, n2 = c[0]
        for p, r in zip(pts, wr, strict=False):
            wp = np.asarray(p, dtype=float)
            txt = (
                f"trace: {ctype} windows\n"
                f"connection maxima: {n1}-{n2}\n"
                f"type: {ctype}\n"
                f"size: {float(r):.3f} A"
            )
            if ctype == "internal":
                internal_windows.append(wp)
                internal_hover.append(txt)
            elif ctype == "periodic":
                periodic_windows.append(wp)
                periodic_hover.append(txt)
    out: dict[str, Any] = {
        "internal_positions": [],
        "internal_hover": internal_hover,
        "internal_n_points": 0,
        "periodic_positions": [],
        "periodic_hover": periodic_hover,
        "periodic_n_points": 0,
    }
    if internal_windows:
        w = np.vstack(internal_windows)
        out["internal_positions"] = w.reshape(-1).tolist()
        out["internal_n_points"] = int(len(w))
    if periodic_windows:
        w = np.vstack(periodic_windows)
        out["periodic_positions"] = w.reshape(-1).tolist()
        out["periodic_n_points"] = int(len(w))
    return out


def _serialize_connection_traces(
    seg: dict[str, Any], atoms: Atoms
) -> dict[str, Any]:
    conns = list(seg.get("connections", []))
    maxima = np.asarray(seg.get("all_maxima", []), dtype=float)
    if len(conns) == 0 or maxima.size == 0:
        return {
            "internal_segments": [],
            "internal_hover": [],
            "periodic_segments": [],
            "periodic_hover": [],
            "internal_points": [],
            "internal_points_hover": [],
            "periodic_points": [],
            "periodic_points_hover": [],
        }
    shape = np.asarray(seg["dgrid"].shape, dtype=float)
    cell = atoms.get_cell().array.T
    node_positions: dict[int, np.ndarray] = {}
    nodes: set[int] = set()
    for c in conns:
        nodes.update(c[0])
    for node in sorted(nodes):
        idx = int(node) - 1
        if idx < 0 or idx >= len(maxima):
            continue
        frac = maxima[idx] / shape
        node_positions[int(node)] = np.dot(cell, frac.T).T

    internal_segments: list[list[float]] = []
    internal_hover: list[str] = []
    periodic_segments: list[list[float]] = []
    periodic_hover: list[str] = []
    internal_nodes: set[int] = set()
    periodic_nodes: set[int] = set()
    for c in conns:
        if len(c) < 5 or not c[1]:
            continue
        n1, n2 = [int(x) for x in c[0]]
        p1 = node_positions.get(n1)
        p2 = node_positions.get(n2)
        if p1 is None or p2 is None:
            continue
        ctype = str(c[2])
        pts = c[3] if c[3] is not None else []
        wr = c[4] if c[4] is not None else []
        if ctype == "internal":
            internal_nodes.update([n1, n2])
            internal_segments.append(
                [
                    float(p1[0]), float(p1[1]), float(p1[2]),
                    float(p2[0]), float(p2[1]), float(p2[2]),
                ]
            )
            txt = f"trace: internal connections\nedge: {n1}-{n2}"
            if len(pts) and len(wr):
                txt += f"\nwindow radius: {float(wr[0]):.3f} A"
            internal_hover.append(txt)
        elif ctype == "periodic":
            periodic_nodes.update([n1, n2])
            periodic_segments.append(
                [
                    float(p1[0]), float(p1[1]), float(p1[2]),
                    float(p2[0]), float(p2[1]), float(p2[2]),
                ]
            )
            txt = f"trace: periodic connections\nedge: {n1}-{n2}"
            if len(pts) and len(wr):
                txt += f"\nwindow radius: {float(wr[0]):.3f} A"
            periodic_hover.append(txt)

    def _pack_nodes(node_ids: set[int], label: str) -> tuple[list[float], list[str]]:
        if not node_ids:
            return [], []
        arr: list[np.ndarray] = []
        hv: list[str] = []
        for n in sorted(node_ids):
            p = node_positions.get(n)
            if p is None:
                continue
            arr.append(p)
            hv.append(
                f"trace: {label}\n"
                f"node: {n}\n"
                f"xyz: ({float(p[0]):.3f}, {float(p[1]):.3f}, {float(p[2]):.3f})"
            )
        if not arr:
            return [], []
        a = np.vstack(arr)
        return a.reshape(-1).tolist(), hv

    internal_points, internal_points_hover = _pack_nodes(internal_nodes, "internal maxima")
    periodic_points, periodic_points_hover = _pack_nodes(periodic_nodes, "periodic maxima")
    return {
        "internal_segments": (
            np.asarray(internal_segments, dtype=float).reshape(-1).tolist()
            if internal_segments else []
        ),
        "internal_hover": internal_hover,
        "periodic_segments": (
            np.asarray(periodic_segments, dtype=float).reshape(-1).tolist()
            if periodic_segments else []
        ),
        "periodic_hover": periodic_hover,
        "internal_points": internal_points,
        "internal_points_hover": internal_points_hover,
        "periodic_points": periodic_points,
        "periodic_points_hover": periodic_points_hover,
    }


# ---------------------------------------------------------------------------
# I/O adapters
# ---------------------------------------------------------------------------

def _snapshot_dict(
    atoms: Atoms,
    *,
    path: str | Path,
    frame_index: int,
    source_format: str,
) -> dict[str, Any]:
    return {
        "atoms": atoms,
        "bonds": None,
        "cell": atoms.get_cell().array if atoms.cell else None,
        "dgrid": None,
        "regions": None,
        "connections": None,
        "meta": {
            "frame_index": frame_index,
            "source": str(path),
            "format": source_format,
        },
    }


def _normalize_element_symbol(symbol: str) -> str:
    s = str(symbol).strip()
    if not s:
        return "X"
    if len(s) == 1:
        return s.upper()
    return s[0].upper() + s[1:].lower()


def read_movie_snapshots(
    path: str | Path,
    *,
    movie_format: str = "auto",
    max_frames: int | None = None,
    frame_stride: int = 1,
    position_decimals: int = 4,
) -> list[dict[str, Any]]:
    """Load MD snapshots as canonical ``snapshot`` dicts (each with an ``ase.Atoms``).

    Parameters
    ----------
    path
        Trajectory or RASPA multi-model PDB file.
    movie_format
        ``"auto"`` — ``.pdb`` uses RASPA parser; everything else uses ASE indexed reads
        (one frame at a time, low peak memory).
        ``"raspa_pdb"`` / ``"ase"`` force a specific backend.
    max_frames
        Maximum number of *kept* frames (after stride), not counting skipped indices.
    frame_stride
        Keep every Nth source frame (reduces output size and work in the browser).
    position_decimals
        Round Cartesian positions before building ``Atoms`` (smaller JSON when embedded).
    """
    from .data import read_raspa_pdb

    p = Path(path).expanduser().resolve()
    stride = max(int(frame_stride), 1)
    decimals = max(int(position_decimals), 0)
    fmt = (movie_format or "auto").lower()
    if fmt == "auto":
        fmt = "raspa_pdb" if p.suffix.lower() == ".pdb" else "ase"

    snapshots: list[dict[str, Any]] = []

    if fmt == "raspa_pdb":
        raw = read_raspa_pdb(str(p))
        coords = raw["coords"]
        symbols = raw["symbols"]
        cells = raw["cells"]
        for i in range(len(coords)):
            if i % stride != 0:
                continue
            xyz = np.round(np.asarray(coords[i], dtype=np.float64), decimals=decimals)
            syms = [_normalize_element_symbol(s) for s in symbols[i]]
            cell = cells[i] if i < len(cells) else cells[-1] if len(cells) else None
            atoms = Atoms(symbols=syms, positions=xyz, cell=cell, pbc=True)
            snapshots.append(_snapshot_dict(atoms, path=p, frame_index=i, source_format="raspa_pdb"))
            if max_frames is not None and len(snapshots) >= max_frames:
                break
    elif fmt == "ase":
        i = 0
        while True:
            try:
                at = read(str(p), index=i)
            except Exception:
                break
            if i % stride == 0:
                pos = np.round(np.asarray(at.get_positions(), dtype=np.float64), decimals=decimals)
                pbc = bool(np.all(at.pbc)) if getattr(at, "pbc", None) is not None else True
                atoms = Atoms(
                    symbols=at.get_chemical_symbols(),
                    positions=pos,
                    cell=at.cell,
                    pbc=pbc,
                )
                snapshots.append(_snapshot_dict(atoms, path=p, frame_index=i, source_format="ase"))
                if max_frames is not None and len(snapshots) >= max_frames:
                    break
            i += 1
    else:
        raise ValueError(f"Unsupported movie_format: {movie_format!r}")

    return snapshots


def read_raspa_pdb_movie(
    path: str | Path,
    max_frames: int | None = None,
    frame_stride: int = 1,
    position_decimals: int = 4,
) -> list[dict[str, Any]]:
    """Parse a RASPA PDB movie into canonical snapshot dictionaries."""
    return read_movie_snapshots(
        path,
        movie_format="raspa_pdb",
        max_frames=max_frames,
        frame_stride=frame_stride,
        position_decimals=position_decimals,
    )


def embedded_movie_dict_from_snapshots(
    snapshots: list[dict[str, Any]],
    *,
    source: str,
    inferred_format: str,
    frame_stride: int,
    position_decimals: int,
) -> dict[str, Any]:
    """Compact per-frame dict for the legacy devtools Three.js HTML (positions + radii only)."""
    frames: list[dict[str, Any]] = []
    for snap in snapshots:
        at = snap["atoms"]
        syms = [_normalize_element_symbol(s) for s in at.get_chemical_symbols()]
        xyz = np.asarray(at.get_positions(), dtype=np.float32)
        atomic_numbers = at.get_atomic_numbers()
        radii = (covalent_radii[atomic_numbers] * 0.8).astype(np.float32)
        frames.append(
            {
                "positions": xyz.reshape(-1).tolist(),
                "symbols": syms,
                "radii": radii.tolist(),
            }
        )
    return {
        "source": source,
        "format": inferred_format,
        "frame_stride": int(frame_stride),
        "position_decimals": int(position_decimals),
        "n_frames": int(len(frames)),
        "n_atoms": int(len(frames[0]["symbols"])) if frames else 0,
        "frames": frames,
        "n_frames_scanned": int(len(snapshots)),
    }


def embedded_movie_dict_from_path(
    path: str | Path,
    *,
    movie_format: str = "auto",
    max_frames: int | None = None,
    frame_stride: int = 1,
    position_decimals: int = 4,
) -> dict[str, Any]:
    """Load trajectory and return the JSON ``movie`` block used by ``md_pdb_threejs_viewer``."""
    p = Path(path).expanduser().resolve()
    snaps = read_movie_snapshots(
        p,
        movie_format=movie_format,
        max_frames=max_frames,
        frame_stride=frame_stride,
        position_decimals=position_decimals,
    )
    inferred = str(snaps[0]["meta"].get("format", "unknown")) if snaps else "none"
    return embedded_movie_dict_from_snapshots(
        snaps,
        source=str(p),
        inferred_format=inferred,
        frame_stride=frame_stride,
        position_decimals=position_decimals,
    )


def _coerce_snapshot_indices(
    snapshot_indices: Sequence[int] | None,
    *,
    n_frames: int,
    max_frames: int | None = None,
) -> list[int]:
    """Return validated frame indices (unique, in-range, ordered)."""
    if n_frames <= 0:
        return []
    if snapshot_indices is None:
        idx = list(range(n_frames))
        if max_frames is not None:
            idx = idx[: max(int(max_frames), 0)]
        return idx
    seen: set[int] = set()
    out: list[int] = []
    for raw in snapshot_indices:
        i = int(raw)
        if i < 0 or i >= n_frames or i in seen:
            continue
        seen.add(i)
        out.append(i)
    out.sort()
    if max_frames is not None:
        out = out[: max(int(max_frames), 0)]
    return out


def _component_palette_rgb01(i: int) -> tuple[float, float, float]:
    """Simple deterministic palette for component traces."""
    palette = [
        (1.00, 0.41, 0.71),  # pink
        (0.35, 0.78, 0.98),  # cyan
        (1.00, 0.77, 0.20),  # amber
        (0.60, 0.93, 0.45),  # green
        (0.79, 0.64, 1.00),  # violet
        (1.00, 0.55, 0.40),  # salmon
    ]
    return palette[i % len(palette)]


def compute_sparse_segmentation_series(
    framework_snapshots: Sequence[dict[str, Any]],
    snapshot_indices: Sequence[int],
    *,
    spacing: float = 0.25,
    mask_thickness: float = 0.5,
    h: float = 0.2,
    threshold_abs: float = 1.5,
    min_distance: float = 1.5,
    extra_padding: int = 5,
) -> list[dict[str, Any]]:
    """Compute segmentation payloads for sparse framework frame indices."""
    results: list[dict[str, Any]] = []
    for frame_index in snapshot_indices:
        snap = framework_snapshots[frame_index]
        seg = compute_segmentation_for_snapshot(
            snap,
            spacing=spacing,
            mask_thickness=mask_thickness,
            h=h,
            threshold_abs=threshold_abs,
            min_distance=min_distance,
            extra_padding=extra_padding,
        )
        dgrid = np.asarray(seg["dgrid"])
        regions = np.asarray(seg["regions"])
        results.append(
            {
                "frame_index": int(frame_index),
                "segmentation": seg,
                "summary": {
                    "dgrid_shape": [int(v) for v in dgrid.shape],
                    "dgrid_max": float(np.nanmax(dgrid)),
                    "dgrid_min": float(np.nanmin(dgrid)),
                    "n_regions": int(np.unique(regions).size),
                    "n_connections": int(len(seg.get("connections", []))),
                    "n_maxima": int(len(seg.get("all_maxima", []))),
                },
            }
        )
    return results


# ---------------------------------------------------------------------------
# Segmentation helper
# ---------------------------------------------------------------------------

def compute_segmentation_for_snapshot(
    snapshot: dict[str, Any],
    spacing: float = 0.25,
    mask_thickness: float = 0.5,
    h: float = 0.2,
    threshold_abs: float = 1.5,
    min_distance: float = 1.5,
    extra_padding: int = 5,
    use_rust: bool = True,
    **kwargs: Any,
) -> dict[str, Any]:
    """Compute dgrid, regions, maxima, connections for a single snapshot."""
    from .rust_fun import do_periodic_segmentation_rust
    from .compute import dgrid_from_atoms_cpu

    atoms = snapshot["atoms"]
    dgrid = snapshot.get("dgrid")
    if dgrid is None:
        dgrid = dgrid_from_atoms_cpu(atoms, spacing=spacing, block_size=20000).astype(
            np.float32
        )

    if use_rust:
        seg = do_periodic_segmentation_rust(
            dgrid,
            atoms,
            spacing=spacing,
            mask_thickness=mask_thickness,
            h=h,
            threshold_abs=threshold_abs,
            min_distance=min_distance,
            return_dgrid=True,
            return_padded_regions=False,
            extra_padding=extra_padding,
            **kwargs,
        )
    else:
        raise NotImplementedError("Only Rust backend is supported for now.")

    return {
        "dgrid": seg.get("dgrid", dgrid),
        "regions": seg["regions"],
        "maxima_pbc": seg.get("maxima_pbc", []),
        "connections": seg.get("connections", []),
        "all_maxima": seg.get("all_maxima", []),
    }


# ---------------------------------------------------------------------------
# Trace builders
# ---------------------------------------------------------------------------

def make_framework_trace(
    atoms: Atoms,
    bonds: str | list | None = "auto",
    style: str | dict | None = "framework",
    name: str = "Framework",
    trace_id: str = "framework",
) -> dict[str, Any]:
    if bonds == "auto":
        from .analyze import update_with_covalent_bonds

        if "covalent_bonds" not in atoms.info:
            atoms = update_with_covalent_bonds(atoms)
        # bonds already stored in atoms.info
    payload = _serialize_framework_mesh(atoms)
    return _make_trace(
        trace_id=trace_id,
        trace_type="mesh",
        name=name,
        data=payload,
        style=style,
        hover=payload.get("vertex_hover"),
        visible=True,
    )


def make_surface_trace(
    dgrid: np.ndarray,
    regions: np.ndarray,
    atoms: Atoms,
    isovalue: float = 1.0,
    step_size: int = 2,
    style: str | dict | None = "surface",
    name: str = "Surface",
    trace_id: str = "surface",
) -> dict[str, Any]:
    seg = {"dgrid": dgrid, "regions": regions}
    payload = _serialize_segmentation_surface_mesh(
        atoms=atoms, seg=seg, isovalue=isovalue, step_size=step_size
    )
    return _make_trace(
        trace_id=trace_id,
        trace_type="mesh",
        name=name,
        data=payload,
        style=style,
        hover=None,
        visible=True,
    )


def make_maxima_trace(
    seg: dict[str, Any],
    atoms: Atoms,
    style: str | dict | None = "internal_maxima",
    name: str = "Maxima",
    trace_id: str = "maxima",
) -> dict[str, Any]:
    payload = _serialize_maxima_trace(seg, atoms)
    # Internal maxima points
    internal = _make_trace(
        trace_id=f"{trace_id}_internal",
        trace_type="points",
        name=f"{name} (internal)",
        data={
            "positions": payload["internal_positions"],
            "colors": payload["internal_colors"],
            "sizes": payload["internal_sizes"],
            "n_points": payload["internal_n_points"],
        },
        style=style,
        hover=payload["internal_hover"],
        visible=True,
    )
    # Outside maxima points
    outside = _make_trace(
        trace_id=f"{trace_id}_outside",
        trace_type="points",
        name=f"{name} (periodic)",
        data={
            "positions": payload["outside_positions"],
            "colors": payload["outside_colors"],
            "sizes": payload["outside_sizes"],
            "n_points": payload["outside_n_points"],
        },
        style="periodic_maxima",
        hover=payload["outside_hover"],
        visible=True,
    )
    return {"internal": internal, "outside": outside, "payload": payload}


def make_windows_trace(
    seg: dict[str, Any],
    style: str | dict | None = "internal_windows",
    name: str = "Windows",
    trace_id: str = "windows",
) -> dict[str, Any]:
    payload = _serialize_windows_trace(seg)
    internal = _make_trace(
        trace_id=f"{trace_id}_internal",
        trace_type="points",
        name=f"{name} (internal)",
        data={
            "positions": payload["internal_positions"],
            "n_points": payload["internal_n_points"],
        },
        style="internal_windows",
        hover=payload["internal_hover"],
        visible=True,
    )
    periodic = _make_trace(
        trace_id=f"{trace_id}_periodic",
        trace_type="points",
        name=f"{name} (periodic)",
        data={
            "positions": payload["periodic_positions"],
            "n_points": payload["periodic_n_points"],
        },
        style="periodic_windows",
        hover=payload["periodic_hover"],
        visible=True,
    )
    return {"internal": internal, "periodic": periodic, "payload": payload}


def make_connections_trace(
    seg: dict[str, Any],
    atoms: Atoms,
    style: str | dict | None = "internal_connections",
    name: str = "Connections",
    trace_id: str = "connections",
) -> dict[str, Any]:
    payload = _serialize_connection_traces(seg, atoms)
    internal_lines = _make_trace(
        trace_id=f"{trace_id}_internal",
        trace_type="line_segments",
        name=f"{name} (internal)",
        data={
            "positions": payload["internal_segments"],
            "n_segments": len(payload["internal_segments"]) // 6,
        },
        style="internal_connections",
        hover=payload["internal_hover"],
        visible=True,
    )
    periodic_lines = _make_trace(
        trace_id=f"{trace_id}_periodic",
        trace_type="line_segments",
        name=f"{name} (periodic)",
        data={
            "positions": payload["periodic_segments"],
            "n_segments": len(payload["periodic_segments"]) // 6,
        },
        style="periodic_connections",
        hover=payload["periodic_hover"],
        visible=True,
    )
    internal_nodes = _make_trace(
        trace_id=f"{trace_id}_internal_nodes",
        trace_type="points",
        name=f"{name} internal nodes",
        data={
            "positions": payload["internal_points"],
            "n_points": len(payload["internal_points"]) // 3,
        },
        style="internal_nodes",
        hover=payload["internal_points_hover"],
        visible=True,
    )
    periodic_nodes = _make_trace(
        trace_id=f"{trace_id}_periodic_nodes",
        trace_type="points",
        name=f"{name} periodic nodes",
        data={
            "positions": payload["periodic_points"],
            "n_points": len(payload["periodic_points"]) // 3,
        },
        style="periodic_nodes",
        hover=payload["periodic_points_hover"],
        visible=True,
    )
    return {
        "internal_lines": internal_lines,
        "periodic_lines": periodic_lines,
        "internal_nodes": internal_nodes,
        "periodic_nodes": periodic_nodes,
        "payload": payload,
    }


def make_guest_trace(
    atoms: Atoms,
    style: str | dict | None = "guests",
    name: str = "Adsorbates",
    trace_id: str = "guests",
) -> dict[str, Any]:
    syms = atoms.get_chemical_symbols()
    colors = np.asarray([_hex_to_rgb01(CPK_COLORS.get(s, "#CCCCCC")) for s in syms], dtype=float)
    atomic_numbers = atoms.get_atomic_numbers()
    radii = covalent_radii[atomic_numbers] * 0.8
    positions = atoms.get_positions().reshape(-1).tolist()
    return _make_trace(
        trace_id=trace_id,
        trace_type="instanced_spheres",
        name=name,
        data={
            "positions": positions,
            "symbols": syms,
            "colors": colors.reshape(-1).tolist(),
            "radii": radii.astype(float).tolist(),
            "n_atoms": int(len(atoms)),
        },
        style=style,
        hover=None,
        visible=True,
    )


def make_cell_trace(
    atoms: Atoms,
    style: str | dict | None = "cell",
    name: str = "Unit cell",
    trace_id: str = "cell",
) -> dict[str, Any]:
    payload = _serialize_cell_lines(atoms)
    return _make_trace(
        trace_id=trace_id,
        trace_type="unit_cell_lines",
        name=name,
        data=payload,
        style=style,
        hover=None,
        visible=True,
    )


# ---------------------------------------------------------------------------
# MovieFigure class
# ---------------------------------------------------------------------------

class MovieFigure:
    """Plotly-like API for building Three.js movie scenes."""

    def __init__(
        self,
        data: list[dict[str, Any]] | None = None,
        layout: dict[str, Any] | None = None,
        frames: list[dict[str, Any]] | None = None,
        groups: dict[str, Any] | None = None,
        styles: dict[str, Any] | None = None,
        meta: dict[str, Any] | None = None,
    ):
        self.data = list(data) if data is not None else []
        self.layout = dict(layout) if layout is not None else {}
        self.frames = list(frames) if frames is not None else []
        self.groups = dict(groups) if groups is not None else {}
        self.styles = dict(DEFAULT_STYLES)
        if styles is not None:
            self.styles.update(styles)
        self.meta = dict(meta) if meta is not None else {}

    # ------------------------------------------------------------------
    # Mutation API
    # ------------------------------------------------------------------

    def add_trace(self, trace: dict[str, Any], frame: int | None = None) -> "MovieFigure":
        if frame is None:
            self.data.append(trace)
        else:
            while len(self.frames) <= frame:
                self.frames.append({"frame_index": len(self.frames), "traces": []})
            self.frames[frame]["traces"].append(trace)
        return self

    def add_traces(self, traces: list[dict[str, Any]], frame: int | None = None) -> "MovieFigure":
        for t in traces:
            self.add_trace(t, frame=frame)
        return self

    def update_traces(
        self,
        patch: dict[str, Any],
        selector: dict[str, Any] | None = None,
        frame: int | None = None,
    ) -> "MovieFigure":
        targets = self.data if frame is None else self.frames[frame]["traces"]
        for t in targets:
            if selector is None or all(t.get(k) == v for k, v in selector.items()):
                t.update(patch)
        return self

    def update_layout(self, patch: dict[str, Any]) -> "MovieFigure":
        self.layout.update(patch)
        return self

    def set_camera(
        self,
        eye: tuple[float, float, float] | None = None,
        center: tuple[float, float, float] | None = None,
        up: tuple[float, float, float] | None = None,
        keyframes: list[dict[str, Any]] | None = None,
    ) -> "MovieFigure":
        if "scene" not in self.layout:
            self.layout["scene"] = {}
        cam = self.layout["scene"].setdefault("camera", {})
        if eye is not None:
            cam["eye"] = list(eye)
        if center is not None:
            cam["center"] = list(center)
        if up is not None:
            cam["up"] = list(up)
        if keyframes is not None:
            cam["keyframes"] = keyframes
        return self

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "version": "1.0",
            "frames": self.frames,
            "static": {"traces": self.data},
            "groups": self.groups,
            "styles": self.styles,
            "layout": self.layout,
            "meta": self.meta,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "MovieFigure":
        return cls(
            data=d.get("static", {}).get("traces", []),
            layout=d.get("layout", {}),
            frames=d.get("frames", []),
            groups=d.get("groups", {}),
            styles=d.get("styles", {}),
            meta=d.get("meta", {}),
        )

    # ------------------------------------------------------------------
    # HTML export
    # ------------------------------------------------------------------

    def to_html(
        self,
        path: str | Path,
        template: str = "threejs_default",
        sidecar: bool = False,
        fps: float = 20.0,
        point_size: float = 0.16,
    ) -> Path:
        out = Path(path).expanduser().resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        scene = self.to_dict()

        if sidecar:
            frames_data = scene.pop("frames")
            sidecar_path = out.with_suffix(".frames.jsonl")
            with open(sidecar_path, "w", encoding="utf-8") as fh:
                for fr in frames_data:
                    fh.write(json.dumps(fr) + "\n")
            scene["sidecar"] = str(sidecar_path.name)

        payload_json = json.dumps(scene)

        # Default camera
        cam = scene.get("layout", {}).get("scene", {}).get("camera", {})
        eye = tuple(cam.get("eye", [2.0, 0.0, 1.35]))
        center = tuple(cam.get("center", [0.0, 0.0, 0.0]))
        up = tuple(cam.get("up", [0.0, 0.0, 1.0]))

        html = _THREEJS_TEMPLATE.render(
            payload_json=payload_json,
            fps=fps,
            eye=eye,
            center=center,
            up=up,
            point_size=point_size,
        )
        out.write_text(html, encoding="utf-8")
        return out


# ---------------------------------------------------------------------------
# High-level builders
# ---------------------------------------------------------------------------

def build_fixed_framework_movie(
    snapshots: list[dict[str, Any]],
    framework_atoms: Atoms | None = None,
    spacing: float = 0.25,
    isovalue: float = 1.0,
    surf_step_size: int = 2,
    fps: float = 20.0,
    eye: tuple[float, float, float] = (2.0, 0.0, 1.35),
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    up: tuple[float, float, float] = (0.0, 0.0, 1.0),
    max_frames: int | None = None,
) -> MovieFigure:
    """Build a MovieFigure for fixed-framework + animated adsorbates workflow."""
    n_frames = len(snapshots) if max_frames is None else min(len(snapshots), max_frames)
    snapshots = snapshots[:n_frames]

    # Framework from first snapshot (or provided atoms)
    if framework_atoms is None:
        framework_atoms = snapshots[0]["atoms"]

    fig = MovieFigure()
    fig.set_camera(eye=eye, center=center, up=up)
    fig.update_layout({"fps": fps, "n_frames": n_frames})

    # Static traces
    fig.add_trace(make_framework_trace(framework_atoms))
    fig.add_trace(make_cell_trace(framework_atoms))

    # Compute segmentation once on the framework
    seg = compute_segmentation_for_snapshot(
        {"atoms": framework_atoms, "bonds": None, "cell": None, "dgrid": None, "regions": None, "connections": None, "meta": {}},
        spacing=spacing,
    )

    fig.add_trace(make_surface_trace(seg["dgrid"], seg["regions"], framework_atoms, isovalue=isovalue, step_size=surf_step_size))

    maxima = make_maxima_trace(seg, framework_atoms)
    fig.add_trace(maxima["internal"])
    fig.add_trace(maxima["outside"])

    windows = make_windows_trace(seg)
    fig.add_trace(windows["internal"])
    fig.add_trace(windows["periodic"])

    conns = make_connections_trace(seg, framework_atoms)
    fig.add_trace(conns["internal_lines"])
    fig.add_trace(conns["periodic_lines"])
    fig.add_trace(conns["internal_nodes"])
    fig.add_trace(conns["periodic_nodes"])

    # Groups for toggling
    fig.groups = {
        "internal_pore_graph": {
            "traces": ["maxima_internal", "windows_internal", "connections_internal", "connections_internal_nodes"],
            "label": "Internal Pore Graph",
        },
        "periodic_pore_graph": {
            "traces": ["maxima_outside", "windows_periodic", "connections_periodic", "connections_periodic_nodes"],
            "label": "Periodic Pore Graph",
        },
        "guests": {
            "traces": ["guests"],
            "label": "Adsorbates",
        },
    }

    # Per-frame guest traces
    for i, snap in enumerate(snapshots):
        guest = make_guest_trace(snap["atoms"])
        fig.add_trace(guest, frame=i)

    fig.meta = {
        "workflow": "fixed_framework",
        "n_frames": n_frames,
        "n_atoms_framework": int(len(framework_atoms)),
    }
    return fig


# ---------------------------------------------------------------------------
# Jinja2 template loader
# ---------------------------------------------------------------------------

class _TemplateLoader:
    """Minimal template loader that falls back to an embedded default."""

    def __init__(self) -> None:
        self._templates: dict[str, str] = {}

    def _default_threejs(self) -> str:
        return _DEFAULT_THREEJS_TEMPLATE

    def get(self, name: str) -> str:
        if name not in self._templates:
            if name == "threejs_default":
                self._templates[name] = self._default_threejs()
            else:
                raise ValueError(f"Unknown template: {name}")
        return self._templates[name]

    def render(self, **kwargs: Any) -> str:
        # Simple string substitution fallback
        template = self.get(kwargs.pop("template", "threejs_default"))
        for key, val in kwargs.items():
            placeholder = f"{{{{{key}}}}}"
            if isinstance(val, (tuple, list)):
                val = json.dumps(val)
            elif isinstance(val, (dict,)):
                val = json.dumps(val)
            template = template.replace(placeholder, str(val))
        return template


_THREEJS_TEMPLATE = _TemplateLoader()

# ---------------------------------------------------------------------------
# Embedded default Three.js template
# ---------------------------------------------------------------------------

_DEFAULT_THREEJS_TEMPLATE = r'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Framework + Segmentation + MD Snapshots</title>
  <style>
    html, body { margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden;
      background: radial-gradient(circle at center, #7a7a7a 0%, #111111 100%);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #efefef; }
    #app { position: fixed; inset: 0; }
    #hud { position: fixed; left: 12px; top: 12px; width: min(320px, calc(100vw - 24px));
      max-height: calc(100vh - 24px); overflow: auto;
      background: rgba(16, 16, 16, 0.12); border: 1px solid rgba(255, 255, 255, 0.15);
      border-radius: 10px; padding: 10px; display: grid; gap: 8px; backdrop-filter: blur(2px); }
    #row { display: flex; align-items: center; gap: 6px; flex-wrap: wrap; }
    #frameSlider { width: 100%; }
    .checkrow { display: grid; gap: 6px; align-items: start; }
    label { user-select: none; font-size: 13px; text-shadow: 0 1px 2px rgba(0,0,0,0.5); }
    #tooltip { position: fixed; pointer-events: none; display: none;
      background: rgba(0,0,0,0.85); border: 1px solid #666; border-radius: 8px;
      padding: 6px 8px; font-size: 12px; white-space: pre; }
    button { background: rgba(0,0,0,0.0); color: #f1f1f1; border: 1px solid rgba(255,255,255,0.4);
      border-radius: 6px; padding: 4px 8px; cursor: pointer; text-shadow: 0 1px 2px rgba(0,0,0,0.65); }
    button:hover { background: rgba(255,255,255,0.1); }
    #status, #meta { font-size: 12px; text-shadow: 0 1px 2px rgba(0,0,0,0.65); }
    #meta { display: block; opacity: 0.92; line-height: 1.35; }
    input[type="checkbox"] { accent-color: #b9d7ff; }
  </style>
</head>
<body>
  <div id="app"></div>
  <div id="tooltip"></div>
  <div id="hud">
    <div id="row">
      <button id="playPause">Play</button>
      <button id="prevFrame">Prev</button>
      <button id="nextFrame">Next</button>
      <span id="status">Frame 0</span>
    </div>
    <span id="meta"></span>
    <div class="checkrow">
      <button id="pivotScene">Pivot: Scene center</button>
      <button id="pivotOrigin">Pivot: Origin (0,0,0)</button>
    </div>
    <div class="checkrow">
      <label><input type="checkbox" id="showCell" checked /> Unit cell</label>
      <label><input type="checkbox" id="showFramework" checked /> Framework mesh</label>
      <label><input type="checkbox" id="showSurface" checked /> Segmentation surface (regions)</label>
      <label><input type="checkbox" id="showGuests" checked /> Snapshot adsorbates (mesh spheres)</label>
      <label><input type="checkbox" id="showWindows" checked /> Internal windows</label>
      <label><input type="checkbox" id="showPeriodicWindows" checked /> Periodic windows</label>
      <label><input type="checkbox" id="showInternalConnections" checked /> Internal connections</label>
      <label><input type="checkbox" id="showPeriodicConnections" checked /> Periodic connections</label>
    </div>
    <input id="frameSlider" type="range" min="0" max="0" step="1" value="0" />
  </div>

  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
  <script>
    const DATA = {{payload_json}};
    const MOV = DATA.static && DATA.static.traces ? DATA.static.traces : [];
    const FRAMES = DATA.frames || [];
    const STYLES = DATA.styles || {};
    let guestSphereSegments = 14;
    const _gtProbe = FRAMES.length && FRAMES[0].traces
      ? FRAMES[0].traces.find(t => t.id === "guests" || t.id.startsWith("guests_"))
      : null;
    if (_gtProbe && _gtProbe.data && _gtProbe.data.n_atoms > 2000) guestSphereSegments = 6;
    else if (_gtProbe && _gtProbe.data && _gtProbe.data.n_atoms > 500) guestSphereSegments = 8;
    const app = document.getElementById("app");
    const slider = document.getElementById("frameSlider");
    const statusEl = document.getElementById("status");
    const metaEl = document.getElementById("meta");
    const playBtn = document.getElementById("playPause");
    const prevBtn = document.getElementById("prevFrame");
    const nextBtn = document.getElementById("nextFrame");
    const tooltip = document.getElementById("tooltip");

    const showCell = document.getElementById("showCell");
    const showFramework = document.getElementById("showFramework");
    const showSurface = document.getElementById("showSurface");
    const showGuests = document.getElementById("showGuests");
    const showWindows = document.getElementById("showWindows");
    const showPeriodicWindows = document.getElementById("showPeriodicWindows");
    const showInternalConnections = document.getElementById("showInternalConnections");
    const showPeriodicConnections = document.getElementById("showPeriodicConnections");
    const pivotSceneBtn = document.getElementById("pivotScene");
    const pivotOriginBtn = document.getElementById("pivotOrigin");

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    renderer.setSize(window.innerWidth, window.innerHeight);
    app.appendChild(renderer.domElement);

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x111111);
    const camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.01, 2000);
    const eyeArr = {{eye}};
    const centerArr = {{center}};
    const upArr = {{up}};
    camera.position.set(eyeArr[0], eyeArr[1], eyeArr[2]);
    camera.up.set(upArr[0], upArr[1], upArr[2]);
    camera.lookAt(new THREE.Vector3(centerArr[0], centerArr[1], centerArr[2]));
    const controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.target.set(centerArr[0], centerArr[1], centerArr[2]);
    controls.update();

    const ambient = new THREE.AmbientLight(0xffffff, 0.4);
    scene.add(ambient);
    const dirKey = new THREE.DirectionalLight(0xffffff, 0.85);
    dirKey.position.set(5, 10, 7.5);
    scene.add(dirKey);
    const dirFill = new THREE.DirectionalLight(0xffffff, 0.55);
    dirFill.position.set(-5, -9, -7.5);
    scene.add(dirFill);

    const objectsById = {};

    function styleFor(trace) {
      const s = trace && trace.style;
      if (s && typeof s === "object" && !Array.isArray(s)) return s;
      if (typeof s === "string" && STYLES[s]) return STYLES[s];
      return {};
    }

    function threeColorFromStyle(c, fallbackHex) {
      if (c === undefined || c === null) return new THREE.Color(fallbackHex !== undefined ? fallbackHex : 0xffffff);
      if (typeof c === "number") return new THREE.Color(c);
      if (Array.isArray(c) && c.length >= 3) return new THREE.Color(c[0], c[1], c[2]);
      return new THREE.Color(fallbackHex !== undefined ? fallbackHex : 0xffffff);
    }

    function pointsBaseSize(trace) {
      const style = styleFor(trace);
      let size = style.size !== undefined ? style.size : 1.0;
      if (trace.data && trace.data.sizes && trace.data.sizes.length) {
        let mx = 0;
        for (let i = 0; i < trace.data.sizes.length; i++) {
          const v = Number(trace.data.sizes[i]);
          if (!isNaN(v)) mx = Math.max(mx, v);
        }
        if (mx > 0) size = mx;
      }
      return size;
    }

    function addMesh(trace) {
      if (!trace.data.positions || !trace.data.positions.length || !trace.data.indices || !trace.data.indices.length) {
        console.warn("Skipping empty mesh trace:", trace.id);
        return null;
      }
      const geo = new THREE.BufferGeometry();
      geo.setAttribute("position", new THREE.Float32BufferAttribute(trace.data.positions, 3));
      if (trace.data.colors && trace.data.colors.length) {
        geo.setAttribute("color", new THREE.Float32BufferAttribute(trace.data.colors, 3));
      }
      geo.setIndex(trace.data.indices);
      geo.computeVertexNormals();
      const style = styleFor(trace);
      const mat = new THREE.MeshStandardMaterial({
        vertexColors: !!(trace.data.colors && trace.data.colors.length),
        transparent: true,
        opacity: style.opacity !== undefined ? style.opacity : 0.95,
        roughness: style.roughness !== undefined ? style.roughness : 0.38,
        metalness: style.metalness !== undefined ? style.metalness : 0.0,
        side: style.side === "FrontSide" ? THREE.FrontSide : THREE.DoubleSide,
      });
      const mesh = new THREE.Mesh(geo, mat);
      mesh.name = trace.id;
      scene.add(mesh);
      objectsById[trace.id] = mesh;
      return mesh;
    }

    function addLineSegments(trace, colorOverride) {
      const pts = trace.data.positions || [];
      if (!pts.length) { console.warn("Skipping empty line trace:", trace.id); return null; }
      const geo = new THREE.BufferGeometry();
      geo.setAttribute("position", new THREE.Float32BufferAttribute(pts, 3));
      const style = styleFor(trace);
      const raw = colorOverride !== undefined ? colorOverride : (style.color !== undefined ? style.color : 0xffffff);
      const color = threeColorFromStyle(raw, 0xffffff);
      let mat;
      if (style.dashSize !== undefined) {
        mat = new THREE.LineDashedMaterial({
          color: color,
          dashSize: style.dashSize,
          gapSize: style.gapSize !== undefined ? style.gapSize : 0.5,
          linewidth: style.linewidth || 1
        });
      } else {
        mat = new THREE.LineBasicMaterial({ color: color, linewidth: style.linewidth || 1 });
      }
      const lines = new THREE.LineSegments(geo, mat);
      if (style.dashSize !== undefined) lines.computeLineDistances();
      lines.name = trace.id;
      scene.add(lines);
      objectsById[trace.id] = lines;
      return lines;
    }

    function addPoints(trace, colorOverride, sizeOverride) {
      const pts = trace.data.positions || [];
      if (!pts.length) { console.warn("Skipping empty points trace:", trace.id); return null; }
      const geo = new THREE.BufferGeometry();
      geo.setAttribute("position", new THREE.Float32BufferAttribute(pts, 3));
      if (trace.data.colors && trace.data.colors.length) {
        geo.setAttribute("color", new THREE.Float32BufferAttribute(trace.data.colors, 3));
      }
      const style = styleFor(trace);
      const size = sizeOverride !== undefined ? sizeOverride : pointsBaseSize(trace);
      const mat = new THREE.PointsMaterial({
        size: size,
        vertexColors: !!(trace.data.colors && trace.data.colors.length),
        sizeAttenuation: true,
        transparent: true,
        opacity: 0.95,
      });
      if (colorOverride !== undefined && (!trace.data.colors || !trace.data.colors.length)) {
        mat.color.copy(threeColorFromStyle(colorOverride, 0xffffff));
      } else if (!trace.data.colors || !trace.data.colors.length) {
        mat.color.copy(threeColorFromStyle(style.color, 0xffffff));
      }
      const points = new THREE.Points(geo, mat);
      points.name = trace.id;
      scene.add(points);
      objectsById[trace.id] = points;
      return points;
    }

    function addInstancedSpheres(trace, segLat, segLon) {
      const n = trace.data.n_atoms || 0;
      if (!n) { console.warn("Skipping empty instanced trace:", trace.id); return null; }
      const sw = segLat || 14;
      const sh = segLon || 14;
      const geo = new THREE.SphereBufferGeometry(1.0, sw, sh);
      const style = styleFor(trace);
      const c = style.color;
      const col = Array.isArray(c) && c.length >= 3
        ? new THREE.Color(c[0], c[1], c[2])
        : (typeof c === "number" ? new THREE.Color(c) : new THREE.Color(0xff69b4));
      const mat = new THREE.MeshPhongMaterial({
        color: col,
        shininess: style.shininess !== undefined ? style.shininess : 10,
        specular: style.specular !== undefined ? style.specular : 0x151515,
      });
      const mesh = new THREE.InstancedMesh(geo, mat, n);
      mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
      mesh.name = trace.id;
      scene.add(mesh);
      objectsById[trace.id] = mesh;
      return mesh;
    }

    function addUnitCellLines(trace) {
      const st = styleFor(trace);
      const c = st.color !== undefined ? st.color : 0xffa500;
      return addLineSegments(trace, c);
    }

    function replaceMeshGeometry(mesh, trace) {
      const oldGeo = mesh.geometry;
      const geo = new THREE.BufferGeometry();
      const pos = trace.data.positions || [];
      geo.setAttribute("position", new THREE.Float32BufferAttribute(new Float32Array(pos), 3));
      if (trace.data.colors && trace.data.colors.length) {
        geo.setAttribute("color", new THREE.Float32BufferAttribute(new Float32Array(trace.data.colors), 3));
      }
      geo.setIndex(trace.data.indices);
      geo.computeVertexNormals();
      mesh.geometry = geo;
      if (oldGeo) oldGeo.dispose();
    }

    function replaceLineGeometry(lines, trace) {
      const oldGeo = lines.geometry;
      const geo = new THREE.BufferGeometry();
      const pts = trace.data.positions || [];
      geo.setAttribute("position", new THREE.Float32BufferAttribute(new Float32Array(pts), 3));
      lines.geometry = geo;
      if (oldGeo) oldGeo.dispose();
      if (lines.material && lines.material.dashSize !== undefined) lines.computeLineDistances();
    }

    function replacePointsGeometry(points, trace) {
      const oldGeo = points.geometry;
      const geo = new THREE.BufferGeometry();
      const pts = trace.data.positions || [];
      geo.setAttribute("position", new THREE.Float32BufferAttribute(new Float32Array(pts), 3));
      if (trace.data.colors && trace.data.colors.length) {
        geo.setAttribute("color", new THREE.Float32BufferAttribute(new Float32Array(trace.data.colors), 3));
      }
      points.geometry = geo;
      if (oldGeo) oldGeo.dispose();
      const style = styleFor(trace);
      points.material.size = pointsBaseSize(trace);
      points.material.vertexColors = !!(trace.data.colors && trace.data.colors.length);
      if (!trace.data.colors || !trace.data.colors.length) {
        points.material.color.copy(threeColorFromStyle(style.color, 0xffffff));
      }
    }

    // Build static traces
    const staticTraces = DATA.static && DATA.static.traces ? DATA.static.traces : [];
    staticTraces.forEach(trace => {
      if (trace.type === "mesh") addMesh(trace);
      else if (trace.type === "line_segments") addLineSegments(trace);
      else if (trace.type === "points") addPoints(trace);
      else if (trace.type === "unit_cell_lines") addUnitCellLines(trace);
      else if (trace.type === "instanced_spheres") {
        if (trace.id === "guests" || trace.id.startsWith("guests_")) {
          return;
        }
        addInstancedSpheres(trace, 14, 14);
      }
    });

    // Traces that exist only on animation frames (e.g. flexible framework) — create once from frame 0
    const firstFrameTraces = (FRAMES.length > 0 && FRAMES[0].traces) ? FRAMES[0].traces : [];
    firstFrameTraces.forEach(trace => {
      if (objectsById[trace.id]) return;
      if (trace.type === "mesh") addMesh(trace);
      else if (trace.type === "line_segments") addLineSegments(trace);
      else if (trace.type === "points") addPoints(trace);
      else if (trace.type === "unit_cell_lines") addUnitCellLines(trace);
    });

    // Guest instanced meshes for animation (supports multi-component ids like guests_CO2)
    const guestMeshes = {};
    const lastRenderedGuestFrameById = {};
    const guestTraceIds = new Set();
    staticTraces.forEach(t => {
      if (t.type === "instanced_spheres" && (t.id === "guests" || t.id.startsWith("guests_"))) {
        guestTraceIds.add(t.id);
      }
    });
    firstFrameTraces.forEach(t => {
      if (t.type === "instanced_spheres" && (t.id === "guests" || t.id.startsWith("guests_"))) {
        guestTraceIds.add(t.id);
      }
    });
    guestTraceIds.forEach(id => {
      const staticTrace = staticTraces.find(t => t.id === id);
      const frameTrace = firstFrameTraces.find(t => t.id === id);
      const src = staticTrace || frameTrace;
      if (!src) return;
      guestMeshes[id] = addInstancedSpheres(src, guestSphereSegments, guestSphereSegments);
      guestMeshes[id].userData.maxInstances = src.data.n_atoms || 0;
      lastRenderedGuestFrameById[id] = -999;
    });

    // Fallback: if nothing expanded the bounding box (no objects created),
    // force a minimal scene center so camera doesn't look at NaN.
    const sceneBox = new THREE.Box3();
    Object.values(objectsById).forEach(obj => {
      if (obj && obj.geometry) {
        try { sceneBox.expandByObject(obj); } catch (e) {}
      }
    });
    // Also expand by unit cell corners if cell trace exists
    const cellObj = objectsById["cell"];
    if (cellObj && cellObj.geometry && cellObj.geometry.attributes && cellObj.geometry.attributes.position) {
      const pos = cellObj.geometry.attributes.position.array;
      for (let i = 0; i < pos.length; i += 3) {
        sceneBox.expandByPoint(new THREE.Vector3(pos[i], pos[i+1], pos[i+2]));
      }
    }
    const sceneCenter = new THREE.Vector3();
    if (sceneBox.isEmpty()) {
      sceneCenter.set(0, 0, 0);
    } else {
      sceneBox.getCenter(sceneCenter);
    }

    let frame = 0;
    let playing = false;
    let lastT = 0;
    const frameDtMs = 1000.0 / Math.max({{fps}}, 0.01);

    slider.max = String(Math.max(FRAMES.length - 1, 0));
    metaEl.textContent = `Frames: ${FRAMES.length}`;
    controls.target.copy(sceneCenter);
    controls.update();

    function setFrame(i) {
      if (FRAMES.length <= 0) {
        statusEl.textContent = "No frames";
        slider.value = "0";
        return;
      }
      frame = Math.max(0, Math.min(i, FRAMES.length - 1));
      const fr = FRAMES[frame];
      const traces = fr && fr.traces ? fr.traces : [];
      Object.keys(guestMeshes).forEach(guestId => {
        const guestTrace = traces.find(t => t.id === guestId);
        let guestMesh = guestMeshes[guestId];
        if (!guestTrace || !guestMesh) return;
        if (frame === lastRenderedGuestFrameById[guestId]) return;
        const positions = guestTrace.data.positions;
        const radii = guestTrace.data.radii;
        const n = guestTrace.data.n_atoms || 0;
        const cap = guestMesh.userData.maxInstances || guestMesh.count;
        if (n > cap) {
          scene.remove(guestMesh);
          if (guestMesh.geometry) guestMesh.geometry.dispose();
          guestMesh.material.dispose();
          guestMeshes[guestId] = addInstancedSpheres(guestTrace, guestSphereSegments, guestSphereSegments);
          guestMesh = guestMeshes[guestId];
          guestMesh.userData.maxInstances = n;
        }
        const dummy = new THREE.Object3D();
        for (let a = 0; a < n; a++) {
          const x = positions[a * 3], y = positions[a * 3 + 1], z = positions[a * 3 + 2];
          const r = (radii && radii[a]) || 0.6;
          dummy.position.set(x, y, z);
          dummy.scale.set(r, r, r);
          dummy.updateMatrix();
          guestMesh.setMatrixAt(a, dummy.matrix);
        }
        guestMesh.count = n;
        guestMesh.instanceMatrix.needsUpdate = true;
        lastRenderedGuestFrameById[guestId] = frame;
      });
      traces.forEach(trace => {
        if (trace.id === "guests" || (trace.id && trace.id.startsWith("guests_"))) return;
        const obj = objectsById[trace.id];
        if (!obj) return;
        if (trace.type === "mesh") replaceMeshGeometry(obj, trace);
        else if (trace.type === "line_segments" || trace.type === "unit_cell_lines") replaceLineGeometry(obj, trace);
        else if (trace.type === "points") replacePointsGeometry(obj, trace);
      });
      statusEl.textContent = `Frame ${frame}`;
      slider.value = String(frame);
    }

    slider.addEventListener("input", () => setFrame(parseInt(slider.value, 10)));
    prevBtn.addEventListener("click", () => setFrame(frame - 1));
    nextBtn.addEventListener("click", () => setFrame(frame + 1));
    playBtn.addEventListener("click", () => {
      playing = !playing;
      playBtn.textContent = playing ? "Pause" : "Play";
    });

    showCell.addEventListener("change", () => {
      const obj = objectsById["cell"];
      if (obj) obj.visible = showCell.checked;
    });
    showFramework.addEventListener("change", () => {
      const obj = objectsById["framework"];
      if (obj) obj.visible = showFramework.checked;
    });
    showSurface.addEventListener("change", () => {
      const obj = objectsById["surface"];
      if (obj) obj.visible = showSurface.checked;
    });
    showGuests.addEventListener("change", () => {
      Object.values(guestMeshes).forEach(m => { if (m) m.visible = showGuests.checked; });
    });
    showWindows.addEventListener("change", () => {
      const obj = objectsById["windows_internal"];
      if (obj) obj.visible = showWindows.checked;
    });
    showPeriodicWindows.addEventListener("change", () => {
      const obj = objectsById["windows_periodic"];
      if (obj) obj.visible = showPeriodicWindows.checked;
    });
    showInternalConnections.addEventListener("change", () => {
      const lineObj = objectsById["connections_internal"];
      const nodeObj = objectsById["connections_internal_nodes"];
      const maxObj = objectsById["maxima_internal"];
      if (lineObj) lineObj.visible = showInternalConnections.checked;
      if (nodeObj) nodeObj.visible = showInternalConnections.checked;
      if (maxObj) maxObj.visible = showInternalConnections.checked;
    });
    showPeriodicConnections.addEventListener("change", () => {
      const lineObj = objectsById["connections_periodic"];
      const nodeObj = objectsById["connections_periodic_nodes"];
      const maxObj = objectsById["maxima_outside"];
      if (lineObj) lineObj.visible = showPeriodicConnections.checked;
      if (nodeObj) nodeObj.visible = showPeriodicConnections.checked;
      if (maxObj) maxObj.visible = showPeriodicConnections.checked;
    });
    pivotSceneBtn.addEventListener("click", () => {
      controls.target.copy(sceneCenter);
      controls.update();
    });
    pivotOriginBtn.addEventListener("click", () => {
      controls.target.set(0, 0, 0);
      controls.update();
    });

    const raycaster = new THREE.Raycaster();
    raycaster.params.Points.threshold = 0.2;
    const mouse = new THREE.Vector2();
    let hoverRaf = null;
    let pendingHoverEv = null;
    function runHoverRaycast(ev) {
      const rect = renderer.domElement.getBoundingClientRect();
      mouse.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
      mouse.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;
      raycaster.setFromCamera(mouse, camera);
      const all = [];
      Object.values(objectsById).forEach(obj => {
        if (obj.visible) {
          const hits = raycaster.intersectObject(obj);
          if (hits.length) all.push(hits[0]);
        }
      });
      all.sort((a, b) => a.distance - b.distance);
      if (!all.length) {
        tooltip.style.display = "none";
        return;
      }
      const hit = all[0];
      let txt = "";
      if (hit.object && hit.object.name && (hit.object.name === "guests" || hit.object.name.startsWith("guests_"))) {
        const idx = hit.instanceId ?? -1;
        const fr = FRAMES[frame];
        const traces = fr && fr.traces ? fr.traces : [];
        const gt = traces.find(t => t.id === hit.object.name);
        if (gt) {
          const p = gt.data.positions.slice(idx * 3, idx * 3 + 3);
          const s = gt.data.symbols[idx] || "?";
          txt = `trace: ${hit.object.name}\nframe: ${frame}\natom: ${idx}\nelement: ${s}\nxyz: (${p[0].toFixed(3)}, ${p[1].toFixed(3)}, ${p[2].toFixed(3)})`;
        }
      } else if (hit.object.name && objectsById[hit.object.name]) {
        const frNow = FRAMES[frame];
        const tracesNow = frNow && frNow.traces ? frNow.traces : [];
        const trace = tracesNow.find(t => t.id === hit.object.name) || staticTraces.find(t => t.id === hit.object.name);
        if (trace && trace.hover && hit.index != null && hit.index >= 0 && hit.index < trace.hover.length) {
          txt = trace.hover[hit.index];
        } else if (hit.face) {
          txt = `trace: ${hit.object.name}\nxyz: (${hit.point.x.toFixed(3)}, ${hit.point.y.toFixed(3)}, ${hit.point.z.toFixed(3)})`;
        }
      }
      tooltip.style.display = "block";
      tooltip.style.left = `${ev.clientX + 12}px`;
      tooltip.style.top = `${ev.clientY + 12}px`;
      tooltip.textContent = txt;
    }
    function onPointerMove(ev) {
      pendingHoverEv = ev;
      if (hoverRaf != null) return;
      hoverRaf = requestAnimationFrame(() => {
        hoverRaf = null;
        if (pendingHoverEv) runHoverRaycast(pendingHoverEv);
        pendingHoverEv = null;
      });
    }
    renderer.domElement.addEventListener("pointermove", onPointerMove);
    renderer.domElement.addEventListener("pointerleave", () => tooltip.style.display = "none");

    function animate(t) {
      requestAnimationFrame(animate);
      if (playing && FRAMES.length > 1 && (t - lastT) >= frameDtMs) {
        setFrame((frame + 1) % FRAMES.length);
        lastT = t;
      }
      controls.update();
      renderer.render(scene, camera);
    }

    window.addEventListener("resize", () => {
      camera.aspect = window.innerWidth / window.innerHeight;
      camera.updateProjectionMatrix();
      renderer.setSize(window.innerWidth, window.innerHeight);
    });

    if (FRAMES.length <= 0) {
      slider.disabled = true;
      playBtn.disabled = true;
      prevBtn.disabled = true;
      nextBtn.disabled = true;
      statusEl.textContent = "No frames";
    }
    Object.values(guestMeshes).forEach(m => { if (m) m.visible = showGuests.checked; });
    setFrame(0);
    requestAnimationFrame(animate);
  </script>
</body>
</html>
'''


# ---------------------------------------------------------------------------
# Convenience CLI / script entry-point (mirror of md_pdb_threejs_viewer.py)
# ---------------------------------------------------------------------------

def make_sparse_multicomponent_movie(
    *,
    component_movie_paths: Mapping[str, str | Path],
    cif_path: str | Path,
    out_html: str | Path = "movie_sparse_multicomponent.html",
    framework_movie_path: str | Path | None = None,
    snapshot_indices: Sequence[int] | None = None,
    movie_format: str = "auto",
    framework_movie_format: str | None = None,
    max_frames: int | None = None,
    frame_stride: int = 1,
    position_decimals: int = 4,
    spacing: float = 0.25,
    isovalue: float = 1.0,
    surf_step_size: int = 2,
    fps: float = 20.0,
    eye: tuple[float, float, float] = (2.0, 0.0, 1.35),
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    up: tuple[float, float, float] = (0.0, 0.0, 1.0),
) -> tuple[Path, dict[str, Any]]:
    """Build sparse-index movie with separate guest traces per component.

    Notes
    -----
    - If ``framework_movie_path`` is set, **framework mesh, unit cell, surface, maxima,
      windows, and connection traces are stored per animation frame** so the viewer can
      show a flexible framework and recomputed segmentation for each sparse snapshot.
    - If ``framework_movie_path`` is omitted, the framework comes from ``cif_path`` and
      stays static (segmentation once); only guest traces animate.
    - Segmentation dictionaries are computed for every selected frame and returned in
      the second tuple element for downstream analysis.
    """
    if not component_movie_paths:
        raise ValueError("component_movie_paths must contain at least one component.")

    comp_snaps: dict[str, list[dict[str, Any]]] = {}
    n_frames_candidates: list[int] = []
    for label, path in component_movie_paths.items():
        snaps = read_movie_snapshots(
            path,
            movie_format=movie_format,
            max_frames=None,
            frame_stride=frame_stride,
            position_decimals=position_decimals,
        )
        if not snaps:
            raise ValueError(f"No snapshots parsed for component {label!r}: {path}")
        comp_snaps[str(label)] = snaps
        n_frames_candidates.append(len(snaps))

    if framework_movie_path is not None:
        fw_snaps = read_movie_snapshots(
            framework_movie_path,
            movie_format=framework_movie_format or movie_format,
            max_frames=None,
            frame_stride=frame_stride,
            position_decimals=position_decimals,
        )
        if not fw_snaps:
            raise ValueError(f"No framework snapshots parsed from: {framework_movie_path}")
        n_frames_candidates.append(len(fw_snaps))
    else:
        framework_atoms = read(str(cif_path))
        fw_snaps = [_snapshot_dict(framework_atoms, path=str(cif_path), frame_index=0, source_format="cif")]

    n_frames_common = int(min(n_frames_candidates))
    selected = _coerce_snapshot_indices(snapshot_indices, n_frames=n_frames_common, max_frames=max_frames)
    if not selected:
        raise ValueError("No valid snapshot indices selected.")

    seg_series = compute_sparse_segmentation_series(
        fw_snaps,
        selected,
        spacing=spacing,
    )
    seg_by_global_idx = {int(entry["frame_index"]): entry["segmentation"] for entry in seg_series}

    fig = MovieFigure()
    fig.set_camera(eye=eye, center=center, up=up)
    fig.update_layout({"fps": fps, "n_frames": len(selected)})

    component_labels = sorted(comp_snaps.keys())
    animate_framework = framework_movie_path is not None

    for local_frame, global_frame_idx in enumerate(selected):
        framework_atoms = fw_snaps[global_frame_idx]["atoms"]
        seg = seg_by_global_idx[int(global_frame_idx)]

        if animate_framework:
            fig.add_trace(make_framework_trace(framework_atoms), frame=local_frame)
            fig.add_trace(make_cell_trace(framework_atoms), frame=local_frame)
            fig.add_trace(
                make_surface_trace(
                    seg["dgrid"],
                    seg["regions"],
                    framework_atoms,
                    isovalue=isovalue,
                    step_size=surf_step_size,
                ),
                frame=local_frame,
            )
            maxima = make_maxima_trace(seg, framework_atoms)
            fig.add_trace(maxima["internal"], frame=local_frame)
            fig.add_trace(maxima["outside"], frame=local_frame)
            windows = make_windows_trace(seg)
            fig.add_trace(windows["internal"], frame=local_frame)
            fig.add_trace(windows["periodic"], frame=local_frame)
            conns = make_connections_trace(seg, framework_atoms)
            fig.add_trace(conns["internal_lines"], frame=local_frame)
            fig.add_trace(conns["periodic_lines"], frame=local_frame)
            fig.add_trace(conns["internal_nodes"], frame=local_frame)
            fig.add_trace(conns["periodic_nodes"], frame=local_frame)
        elif local_frame == 0:
            fig.add_trace(make_framework_trace(framework_atoms))
            fig.add_trace(make_cell_trace(framework_atoms))
            fig.add_trace(
                make_surface_trace(
                    seg["dgrid"],
                    seg["regions"],
                    framework_atoms,
                    isovalue=isovalue,
                    step_size=surf_step_size,
                )
            )
            maxima = make_maxima_trace(seg, framework_atoms)
            fig.add_trace(maxima["internal"])
            fig.add_trace(maxima["outside"])
            windows = make_windows_trace(seg)
            fig.add_trace(windows["internal"])
            fig.add_trace(windows["periodic"])
            conns = make_connections_trace(seg, framework_atoms)
            fig.add_trace(conns["internal_lines"])
            fig.add_trace(conns["periodic_lines"])
            fig.add_trace(conns["internal_nodes"])
            fig.add_trace(conns["periodic_nodes"])

        for ci, label in enumerate(component_labels):
            atoms = comp_snaps[label][global_frame_idx]["atoms"]
            trace = make_guest_trace(
                atoms,
                style={
                    "color": _component_palette_rgb01(ci),
                    "shininess": 10,
                    "specular": 0x151515,
                },
                name=f"Adsorbates ({label})",
                trace_id=f"guests_{label}",
            )
            fig.add_trace(trace, frame=local_frame)

    fig.groups = {
        "guests": {"traces": [f"guests_{lbl}" for lbl in component_labels], "label": "Adsorbates (all components)"}
    }
    fig.meta = {
        "workflow": "sparse_multicomponent",
        "component_labels": component_labels,
        "snapshot_indices": [int(i) for i in selected],
        "n_frames": int(len(selected)),
        "n_frames_common": int(n_frames_common),
        "flexible_framework_per_frame": bool(animate_framework),
    }

    out = fig.to_html(out_html, fps=fps)
    payload = {
        "html_path": str(out),
        "snapshot_indices": [int(i) for i in selected],
        "segmentation_series": seg_series,
    }
    return out, payload


def make_movie_from_trajectory(
    movie_path: str | Path,
    cif_path: str | Path,
    out_html: str | Path = "movie.html",
    *,
    movie_format: str = "auto",
    max_frames: int | None = 200,
    frame_stride: int = 1,
    position_decimals: int = 4,
    spacing: float = 0.25,
    isovalue: float = 1.0,
    surf_step_size: int = 2,
    fps: float = 20.0,
    eye: tuple[float, float, float] = (2.0, 0.0, 1.35),
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    up: tuple[float, float, float] = (0.0, 0.0, 1.0),
) -> Path:
    """Build a fixed-framework movie from any ASE-readable trajectory (or RASPA PDB) + CIF."""
    snapshots = read_movie_snapshots(
        movie_path,
        movie_format=movie_format,
        max_frames=max_frames,
        frame_stride=frame_stride,
        position_decimals=position_decimals,
    )
    framework_atoms = read(str(cif_path))
    fig = build_fixed_framework_movie(
        snapshots=snapshots,
        framework_atoms=framework_atoms,
        spacing=spacing,
        isovalue=isovalue,
        surf_step_size=surf_step_size,
        fps=fps,
        eye=eye,
        center=center,
        up=up,
        max_frames=None,
    )
    return fig.to_html(out_html, fps=fps)


def make_movie_from_raspa_pdb(
    pdb_path: str | Path,
    cif_path: str | Path,
    out_html: str | Path = "movie.html",
    max_frames: int | None = 200,
    frame_stride: int = 1,
    position_decimals: int = 4,
    spacing: float = 0.25,
    isovalue: float = 1.0,
    surf_step_size: int = 2,
    fps: float = 20.0,
    eye: tuple[float, float, float] = (2.0, 0.0, 1.35),
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    up: tuple[float, float, float] = (0.0, 0.0, 1.0),
) -> Path:
    """High-level convenience function to build a fixed-framework movie from RASPA PDB + CIF."""
    return make_movie_from_trajectory(
        movie_path=pdb_path,
        cif_path=cif_path,
        out_html=out_html,
        movie_format="raspa_pdb",
        max_frames=max_frames,
        frame_stride=frame_stride,
        position_decimals=position_decimals,
        spacing=spacing,
        isovalue=isovalue,
        surf_step_size=surf_step_size,
        fps=fps,
        eye=eye,
        center=center,
        up=up,
    )

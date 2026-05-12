#!/usr/bin/env python3
"""
Export trace stages from the org tutorial into GLB files for Blender.

Usage:
  python devtools/scripts/export_tutorial_traces_to_gltf.py \
      --cif examples/TestCIFs/P1-Cu-BTC.cif \
      --outdir outputs/blender_stages
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import plotly.express as px
import trimesh
from ase.io import read

import mofography as mgr


def _to_rgba(color: str | Sequence[float]) -> np.ndarray:
    from matplotlib.colors import to_rgba

    if isinstance(color, str):
        c = color.strip().lower()
        if c.startswith("rgba(") and c.endswith(")"):
            raw = [p.strip() for p in c[5:-1].split(",")]
            if len(raw) == 4:
                r, g, b, a = [float(v) for v in raw]
                if a > 1.0:
                    a = a / 255.0
                return np.array([int(r), int(g), int(b), int(round(a * 255.0))], dtype=np.uint8)
        if c.startswith("rgb(") and c.endswith(")"):
            raw = [p.strip() for p in c[4:-1].split(",")]
            if len(raw) == 3:
                r, g, b = [float(v) for v in raw]
                return np.array([int(r), int(g), int(b), 255], dtype=np.uint8)

    rgba = np.array(to_rgba(color), dtype=float)
    rgba[0:3] *= 255.0
    rgba[3] *= 255.0
    return rgba.astype(np.uint8)


def _mesh_from_plotly_mesh3d(trace, default_color: str = "lightgray") -> trimesh.Trimesh:
    verts = np.column_stack([np.asarray(trace.x), np.asarray(trace.y), np.asarray(trace.z)])
    faces = np.column_stack([np.asarray(trace.i), np.asarray(trace.j), np.asarray(trace.k)]).astype(int)
    finite_v = np.isfinite(verts).all(axis=1)
    if not finite_v.all():
        old_to_new = -np.ones(len(verts), dtype=int)
        kept_idx = np.where(finite_v)[0]
        old_to_new[kept_idx] = np.arange(len(kept_idx))
        keep_f = finite_v[faces].all(axis=1)
        faces = old_to_new[faces[keep_f]]
        verts = verts[finite_v]
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)

    vertexcolor = getattr(trace, "vertexcolor", None)
    if vertexcolor is not None:
        colors = np.array([_to_rgba(c) for c in vertexcolor], dtype=np.uint8)
        if len(colors) == len(verts):
            mesh.visual.vertex_colors = colors
            return mesh

    # Surface traces from plot_surface_dgrid store region-derived coloring via
    # intensity + colorscale; convert that to explicit vertex colors for GLB.
    intensity = getattr(trace, "intensity", None)
    colorscale = getattr(trace, "colorscale", None)
    if intensity is not None and len(intensity) == len(verts) and colorscale is not None:
        vals = np.asarray(intensity, dtype=float)
        cmin = getattr(trace, "cmin", None)
        cmax = getattr(trace, "cmax", None)
        lo = float(np.nanmin(vals) if cmin is None else cmin)
        hi = float(np.nanmax(vals) if cmax is None else cmax)
        if hi <= lo:
            norm = np.zeros_like(vals, dtype=float)
        else:
            norm = np.clip((vals - lo) / (hi - lo), 0.0, 1.0)
        # colorscale can be either a plain list of color strings or
        # [[stop, color], ...]. Convert robustly and interpolate ourselves.
        if isinstance(colorscale, (list, tuple)) and len(colorscale) > 0 and isinstance(colorscale[0], (list, tuple)):
            stops = np.array([float(s[0]) for s in colorscale], dtype=float)
            cols = np.array([_to_rgba(s[1]) for s in colorscale], dtype=float)
        else:
            # qualitative colors -> evenly spaced stops across [0, 1]
            cols = np.array([_to_rgba(c) for c in colorscale], dtype=float)
            stops = np.linspace(0.0, 1.0, len(cols))

        sampled_rgba = np.zeros((len(norm), 4), dtype=float)
        for idx, t in enumerate(norm):
            j = int(np.searchsorted(stops, t, side="right"))
            if j <= 0:
                sampled_rgba[idx] = cols[0]
            elif j >= len(stops):
                sampled_rgba[idx] = cols[-1]
            else:
                t0, t1 = stops[j - 1], stops[j]
                w = 0.0 if t1 <= t0 else (t - t0) / (t1 - t0)
                sampled_rgba[idx] = (1.0 - w) * cols[j - 1] + w * cols[j]

        sampled = sampled_rgba.astype(np.uint8)
        mesh.visual.vertex_colors = sampled
        return mesh

    mesh.visual.vertex_colors = np.tile(_to_rgba(default_color), (len(verts), 1))
    return mesh


def _sphere_cloud_from_scatter(trace, radius: float = 0.18, default_color: str = "gold") -> List[trimesh.Trimesh]:
    x = np.asarray(trace.x, dtype=float)
    y = np.asarray(trace.y, dtype=float)
    z = np.asarray(trace.z, dtype=float)
    xyz = np.column_stack([x, y, z])
    xyz = xyz[np.isfinite(xyz).all(axis=1)]

    marker = getattr(trace, "marker", None)
    marker_size = float(getattr(marker, "size", 8) or 8)
    sphere_radius = max(0.06, radius * marker_size / 8.0)

    color = default_color
    if marker is not None and getattr(marker, "color", None) is not None:
        c = marker.color
        if isinstance(c, str):
            color = c

    rgba = _to_rgba(color)
    meshes: List[trimesh.Trimesh] = []
    for p in xyz:
        s = trimesh.creation.icosphere(subdivisions=2, radius=sphere_radius)
        s.apply_translation(p)
        s.visual.vertex_colors = np.tile(rgba, (len(s.vertices), 1))
        meshes.append(s)
    return meshes


def _cylinders_from_line_trace(trace, radius: float = 0.06, default_color: str = "orange") -> List[trimesh.Trimesh]:
    x = np.asarray(trace.x, dtype=float)
    y = np.asarray(trace.y, dtype=float)
    z = np.asarray(trace.z, dtype=float)
    points = np.column_stack([x, y, z])
    finite = np.isfinite(points).all(axis=1)

    line = getattr(trace, "line", None)
    line_width = float(getattr(line, "width", 4) or 4)
    cyl_radius = max(0.02, radius * line_width / 4.0)

    color = default_color
    if line is not None and getattr(line, "color", None) is not None and isinstance(line.color, str):
        color = line.color
    rgba = _to_rgba(color)

    meshes: List[trimesh.Trimesh] = []
    start_idx = 0
    while start_idx < len(points):
        while start_idx < len(points) and not finite[start_idx]:
            start_idx += 1
        end_idx = start_idx
        while end_idx < len(points) and finite[end_idx]:
            end_idx += 1
        seg = points[start_idx:end_idx]
        for p0, p1 in zip(seg[:-1], seg[1:]):
            if np.linalg.norm(p1 - p0) < 1e-8:
                continue
            c = trimesh.creation.cylinder(radius=cyl_radius, segment=np.vstack([p0, p1]), sections=12)
            c.visual.vertex_colors = np.tile(rgba, (len(c.vertices), 1))
            meshes.append(c)
        start_idx = end_idx + 1
    return meshes


def _trace_to_meshes(trace, name: str) -> List[trimesh.Trimesh]:
    t = trace.type
    if t == "mesh3d":
        return [_mesh_from_plotly_mesh3d(trace)]
    if t == "scatter3d":
        mode = getattr(trace, "mode", "markers")
        if "lines" in mode:
            return _cylinders_from_line_trace(trace)
        return _sphere_cloud_from_scatter(trace)
    raise ValueError(f"Unsupported trace type for '{name}': {t}")


def _build_tutorial_traces(cif_path: Path):
    mat_atoms = read(str(cif_path))

    spacing = 0.25
    h = 0.2
    mask_thickness = 0.5
    threshold_abs = 1.5
    extra_padding = 5
    min_distance = 1.5

    dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=20000).astype(np.float32)
    results = mgr.do_periodic_segmentation_rust(
        dgrid,
        mat_atoms,
        mask_thickness=mask_thickness,
        h=h,
        spacing=spacing,
        return_dgrid=True,
        threshold_abs=threshold_abs,
        min_distance=min_distance,
        return_padded_regions=True,
        extra_padding=extra_padding,
    )

    dgrid = results["dgrid"]
    regions = results["regions"]
    regions_pbc = results["regions_pbc"]
    all_maxima = results["all_maxima"]
    maxima_pbc = results["maxima_pbc"]
    unique_regions = results["unique_regions"]
    maxima = all_maxima[: regions.max()]
    connections = results["connections"]
    connections_trimmed = [c for c in connections if np.any(np.isin(c[0], unique_regions))]

    rag = mgr.rag_from_connections(connections_trimmed, maxima)
    rag = mgr.add_maxima_to_rag(rag, maxima, dgrid[tuple(results["maxima_info"][:, 1:4].T)][: regions.max()], dgrid.shape, mat_atoms)
    rag = mgr.add_pbc_to_rag(rag, results["pbc_groups"])
    from skimage.measure import regionprops

    props = regionprops(regions)
    bins = np.linspace(0, 50, 30)
    rag = mgr.add_psd_to_rag_with_pbc(rag, props, results["pbc_groups"], bins, mat_atoms, dgrid.shape)
    rag = mgr.add_clustering_to_rag_geometry_js(rag, props, results["pbc_groups"], mat_atoms, dgrid.shape)
    regions_cluster = mgr.clustering_grid_from_regions(regions, rag, cluster_key="cluster_geom_pbc")

    cscale = px.colors.qualitative.Pastel
    mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
    isovalue = 1.0

    frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)
    cell_trace = mgr.plot_cell(mat_atoms, color="orange")
    isosurface_segmentation = mgr.plot_surface_dgrid(
        dgrid, stylize=True, isovalue=isovalue, region_labels=regions, ase_atoms=mat_atoms, name="Regions", opacity=1, cscale=cscale
    )
    isosurface_pbc = mgr.plot_surface_dgrid(
        dgrid, stylize=True, isovalue=isovalue, region_labels=regions_pbc, ase_atoms=mat_atoms, name="Regions PBC", opacity=1, cscale=cscale
    )
    isosurface_clusters = mgr.plot_surface_dgrid(
        dgrid, stylize=True, isovalue=isovalue, region_labels=regions_cluster, name="Pocket Types by geometry", ase_atoms=mat_atoms, opacity=1, cscale=cscale
    )
    maxima_trace = mgr.plot_maxima(
        dgrid, maxima, mat_atoms, marker_size=10, show_spheres=False, marker_colorscale=cscale, name="Maxima (distinct)"
    )
    maxima_trace_pbc = mgr.plot_maxima(
        dgrid, maxima_pbc, mat_atoms, show_spheres=False, marker_size=7.5, marker_colorscale=cscale, name="Maxima (pbc)"
    )
    connections_traces = mgr.plot_connections_as_it_is(
        connections_trimmed, all_maxima, regions.shape, mat_atoms, cscale=cscale, plot_disconnected=False, plot_periodic=True
    )
    windows_trace = mgr.plot_windows(
        connections_trimmed, name="Windows", marker_color="yellow", marker_size=10, marker_opacity=0.8
    )

    stage_traces = [
        ("01_framework", list(frame_trace)),
        ("02_cell", [cell_trace]),
        ("03_regions", [isosurface_segmentation]),
        ("04_regions_pbc", [isosurface_pbc]),
        ("05_pocket_types", [isosurface_clusters]),
        ("06_connections", list(connections_traces)),
        ("07_maxima_distinct", [maxima_trace]),
        ("08_windows", [windows_trace]),
        ("09_maxima_pbc", [maxima_trace_pbc]),
    ]
    return stage_traces


def _export_stages(stage_traces: Iterable[Tuple[str, Sequence]], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for stage_name, traces in stage_traces:
        new_meshes: List[Tuple[str, trimesh.Trimesh]] = []
        for i, trace in enumerate(traces):
            meshes = _trace_to_meshes(trace, f"{stage_name}_{i}")
            for j, mesh in enumerate(meshes):
                new_meshes.append((f"{stage_name}_trace{i}_mesh{j}", mesh))

        scene = trimesh.Scene()
        for node_name, mesh in new_meshes:
            scene.add_geometry(mesh, node_name=node_name, geom_name=node_name)
        out_file = outdir / f"{stage_name}.glb"
        scene.export(out_file)
        print(f"Exported {out_file}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Export tutorial traces to staged GLB files for Blender.")
    parser.add_argument("--cif", type=Path, default=Path("examples/TestCIFs/P1-Cu-BTC.cif"))
    parser.add_argument("--outdir", type=Path, default=Path("outputs/blender_stages"))
    args = parser.parse_args()

    stage_traces = _build_tutorial_traces(args.cif)
    _export_stages(stage_traces, args.outdir)
    print("Done.")


if __name__ == "__main__":
    main()

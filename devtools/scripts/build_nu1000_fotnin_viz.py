#!/usr/bin/env python3
"""Granular Plotly + Three.js viz (upstream_features §7) for NU-1000 and FOTNIN.

Defaults: spacing=0.5 Å, mask_thickness=0.05 Å, segmentation otherwise default.
Outputs: examples/viz_outputs/{name}_plotly.html and {name}_threejs.html
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
TEST_CIFS = ROOT / "examples" / "TestCIFs"
OUT_DIR = ROOT / "examples" / "viz_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SPACING = 0.5
MASK_THICKNESS = 0.05
ISOVALUE = 1.0

JOBS = [
    ("NU-1000", TEST_CIFS / "NU-1000.cif"),
    ("FOTNIN_clean", TEST_CIFS / "FOTNIN_clean.cif"),
    ("TOVKOG_clean", TEST_CIFS / "TOVKOG_clean.cif"),
]


def run_pipeline(cif: Path):
    import mofography as mgr
    from mofography.cif_io import read_cif_to_ase

    print(f"\n{'=' * 60}\n{cif.name}\n{'=' * 60}")
    t0 = time.perf_counter()
    atoms = read_cif_to_ase(cif)
    print(f"  load: {len(atoms)} atoms ({time.perf_counter() - t0:.1f}s)")

    t0 = time.perf_counter()
    atoms = mgr.update_with_covalent_bonds(atoms)
    mgr.compute_aromatic_cycles(atoms)
    print(f"  bonds+rings: {time.perf_counter() - t0:.1f}s")

    t0 = time.perf_counter()
    dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=SPACING).astype(np.float32)
    print(f"  dgrid {dgrid.shape}: {time.perf_counter() - t0:.1f}s")

    t0 = time.perf_counter()
    seg = mgr.do_periodic_segmentation_rust(
        dgrid,
        atoms,
        spacing=SPACING,
        mask_thickness=MASK_THICKNESS,
        compute_connections=True,
        return_dgrid=True,
        extra_padding=5,
    )
    seg["dgrid"] = dgrid
    n_reg = len(seg.get("unique_regions", []))
    n_conn = len(seg.get("connections", []))
    print(f"  segmentation: {time.perf_counter() - t0:.1f}s — {n_reg} regions, {n_conn} connections")

    from mofography.helper import trim_connections_to_unique_regions

    connections_trimmed = trim_connections_to_unique_regions(
        seg.get("connections", []),
        seg.get("unique_regions", []),
    )

    return atoms, dgrid, seg, connections_trimmed


def build_plotly(name: str, atoms, seg, connections_trimmed) -> Path:
    import plotly.graph_objects as go
    from mofography.visualize import (
        plot_cell,
        plot_connections_as_it_is,
        plot_framework_as_mesh,
        plot_maxima,
        plot_surface_dgrid,
        plot_windows,
        white_layout,
    )

    out = OUT_DIR / f"{name}_plotly.html"
    dgrid = np.asarray(seg["dgrid"])
    regions = np.asarray(seg["regions"])
    regions_pbc = seg.get("regions_pbc")
    all_maxima = np.asarray(seg.get("all_maxima", []))
    maxima_pbc = seg.get("maxima_pbc")
    cscale = __import__("plotly.express", fromlist=["px"]).colors.qualitative.Pastel

    traces = []
    traces.extend(plot_framework_as_mesh(atoms, stylize=True, show_rings=False))
    traces.append(plot_cell(atoms, color="orange"))

    traces.append(
        plot_surface_dgrid(
            dgrid,
            isovalue=ISOVALUE,
            region_labels=regions,
            ase_atoms=atoms,
            name="Regions",
            opacity=1.0,
            cscale=cscale,
            stepsize=1,
            stylize=True,
        )
    )
    if regions_pbc is not None:
        traces.append(
            plot_surface_dgrid(
                dgrid,
                isovalue=ISOVALUE,
                region_labels=regions_pbc,
                ase_atoms=atoms,
                name="Regions PBC",
                opacity=0.35,
                cscale=cscale,
                stepsize=1,
                stylize=True,
            )
        )

    mx_local = plot_maxima(
        dgrid,
        all_maxima,
        atoms,
        name="Local maxima",
        marker_size=10,
        show_spheres=False,
        marker_colorscale=cscale,
    )
    traces.append(mx_local[0] if isinstance(mx_local, tuple) else mx_local)

    if maxima_pbc is not None and len(maxima_pbc):
        mp = np.asarray(maxima_pbc, dtype=float)
        mx_pbc = plot_maxima(
            dgrid,
            mp,
            atoms,
            name="Maxima (PBC reps)",
            marker_size=12,
            marker_color="cyan",
            show_spheres=False,
        )
        traces.append(mx_pbc[0] if isinstance(mx_pbc, tuple) else mx_pbc)

    traces.extend(
        plot_connections_as_it_is(
            connections_trimmed,
            all_maxima,
            dgrid.shape,
            atoms,
            plot_periodic=True,
            plot_disconnected=False,
        )
    )

    internal_w = [c for c in connections_trimmed if len(c) > 2 and c[2] == "internal" and c[1]]
    periodic_w = [
        c
        for c in connections_trimmed
        if len(c) > 2 and c[2] == "periodic" and c[1] and len(c) > 3 and c[3] is not None
    ]
    if internal_w:
        traces.append(
            plot_windows(
                internal_w,
                name="Windows (internal)",
                marker_color="yellow",
                marker_size=10,
                marker_opacity=0.85,
            )
        )
    if periodic_w:
        traces.append(
            plot_windows(
                periodic_w,
                name="Windows (periodic)",
                marker_color="orange",
                marker_size=10,
                marker_opacity=0.65,
            )
        )

    fig = go.Figure(data=traces, layout=white_layout)
    fig.update_layout(title=dict(text=name))
    fig.write_html(str(out), include_plotlyjs="cdn")
    print(f"  plotly: {out} ({out.stat().st_size / 1024 / 1024:.1f} MB)")
    return out


def build_threejs(name: str, atoms, seg, connections_trimmed) -> Path:
    from mofography import threejs_viz as tjs
    from mofography.movies import make_framework_instanced_trace

    out = OUT_DIR / f"{name}_threejs.html"
    dgrid = np.asarray(seg["dgrid"])
    regions = np.asarray(seg["regions"])
    regions_pbc = seg.get("regions_pbc")
    all_maxima = np.asarray(seg.get("all_maxima", []))
    cscale = tjs.qualitative_colorscale("Pastel")

    seg_dict = {
        "dgrid": dgrid,
        "regions": regions,
        "all_maxima": all_maxima,
        "maxima_pbc": seg.get("maxima_pbc"),
        "maxima_info": seg.get("maxima_info"),
        "maxima_pbc_dict": seg.get("maxima_pbc_dict", {}),
        "connections": connections_trimmed,
    }

    print("  building Three.js figure...")
    t0 = time.perf_counter()
    fig = tjs.Figure3JS(layout="dark", stylize=True)
    fig.set_camera(eye=(2.0, 0.0, 1.35), center=(0, 0, 0), up=(0, 0, 1))

    inst = make_framework_instanced_trace(atoms, representation="ball_and_stick", atom_scale=0.65)
    fig.add_trace(inst["atoms"])
    fig.add_trace(inst["bonds"])

    fig.plot_cell(atoms, color="orange")

    fig.plot_surface_dgrid(
        dgrid,
        isovalue=ISOVALUE,
        ase_atoms=atoms,
        region_labels=regions,
        name="Regions",
        opacity=1.0,
        cscale=cscale,
        stepsize=1,
        stylize=True,
    )
    if regions_pbc is not None:
        fig.plot_surface_dgrid(
            dgrid,
            isovalue=ISOVALUE,
            ase_atoms=atoms,
            region_labels=regions_pbc,
            name="Regions PBC",
            opacity=0.35,
            cscale=cscale,
            stepsize=1,
            stylize=True,
        )

    fig.plot_maxima(
        seg_dict,
        atoms,
        name="Local maxima",
        trace_id="maxima",
        cscale=cscale,
        hide_periodic=False,
    )

    fig.plot_connections_as_it_is(
        connections_trimmed,
        all_maxima,
        dgrid.shape,
        atoms,
        plot_periodic=True,
        plot_disconnected=False,
        name="Connections",
    )

    fig.figure.add_windows(
        seg_dict,
        hide_periodic=False,
        atoms=atoms,
        spacing=SPACING,
        reclassify_far_outside=True,
    )

    fig.to_html(str(out))
    print(f"  threejs: {out} ({out.stat().st_size / 1024 / 1024:.1f} MB, {time.perf_counter() - t0:.1f}s)")
    return out


def main() -> None:
    print(f"Output dir: {OUT_DIR}")
    print(f"spacing={SPACING}  mask_thickness={MASK_THICKNESS}")

    for name, cif in JOBS:
        if not cif.is_file():
            raise FileNotFoundError(cif)
        atoms, _dgrid, seg, conn = run_pipeline(cif)
        build_plotly(name, atoms, seg, conn)
        build_threejs(name, atoms, seg, conn)

    print(f"\nDone. Open files under {OUT_DIR}/")


if __name__ == "__main__":
    main()

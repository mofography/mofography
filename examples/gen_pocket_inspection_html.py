#!/usr/bin/env python3
"""Full pocket / segmentation inspection HTML (Plotly + Three.js).

Layers (toggle in legend): framework mesh, unit cell, regions surface,
regions_pbc surface, pocket-type surface (env_cluster), all_maxima[:regions.max()],
maxima_pbc representatives.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import mofography as mgr
from mofography import Figure3js
from mofography.pocket_typing import (
    DEFAULT_MASK_THICKNESS_A,
    DEFAULT_PROBE_RADIUS_PSD_A,
    DEFAULT_SPACING_A,
    assign_pocket_types_from_segmentation,
)
from mofography.workflow import enrich_segmentation_with_pocket_psds

ISOVALUE = 1.0
_TAUBIN = dict(
    taubin_smooth=True,
    taubin_iterations=3,
    taubin_lambda=0.2175,
    taubin_nu=-0.2175,
    stepsize=1,
)
_SEG = dict(
    mask_thickness=DEFAULT_MASK_THICKNESS_A,
    h=0.2,
    min_distance=1.5,
    threshold_abs=1.5,
    extra_padding=5.0,
    return_dgrid=True,
    compute_connections=True,
)


def build_dgrid(atoms, spacing: float, backend: str) -> np.ndarray:
    if backend == "rtree":
        return mgr.dgrid_from_atoms_rtree_rust(atoms, spacing=spacing).astype(np.float32)
    return mgr.dgrid_from_atoms_rust(atoms, spacing=spacing).astype(np.float32)


def run_segmentation_and_typing(cif: Path, spacing: float, dgrid_backend: str):
    atoms = read(str(cif))
    mat_updated = mgr.update_with_covalent_bonds(atoms)
    t0 = time.time()
    dgrid = build_dgrid(atoms, spacing, dgrid_backend)
    results = mgr.do_periodic_segmentation_rust(dgrid, atoms, spacing=spacing, **_SEG)
    enrich_segmentation_with_pocket_psds(
        results, atoms, probe_radius_psd=DEFAULT_PROBE_RADIUS_PSD_A
    )
    typed = assign_pocket_types_from_segmentation(
        results,
        atoms,
        method="env_cluster",
        probe_radius_psd=DEFAULT_PROBE_RADIUS_PSD_A,
        cif_path=cif,
    )
    if not typed.ok:
        raise RuntimeError(typed.error)
    print(
        f"  seg+typing {time.time() - t0:.1f}s — {typed.n_types} types, "
        f"{typed.n_pockets} pockets"
    )
    return atoms, mat_updated, results, typed


def build_threejs(
    atoms,
    mat_updated,
    results,
    typed,
    *,
    cscale,
    use_mesh_framework: bool = True,
):
    dgrid = results["dgrid"]
    regions = results["regions"]
    regions_pbc = results["regions_pbc"]
    all_maxima = results["all_maxima"]
    maxima_pbc = results["maxima_pbc"]
    n_cap = int(regions.max()) if regions.size else 0
    maxima_local = all_maxima[:n_cap] if n_cap > 0 else all_maxima[:0]

    if use_mesh_framework:
        fig3 = (
            Figure3js(layout="white", stylize=True)
            .plot_framework_as_mesh_3js(mat_updated, stylize=True, show_rings=False)
        )
    else:
        fig3 = (
            Figure3js(layout="white", stylize=True)
            .plot_framework_instanced_3js(mat_updated, atom_scale=0.5)
        )

    fig3.plot_cell_3js(atoms, color="orange", name="Unit cell")

    if typed.regions_cluster is not None:
        fig3.plot_surface_dgrid_3js(
            dgrid,
            ISOVALUE,
            atoms,
            region_labels=typed.regions_cluster,
            cscale=cscale,
            opacity=0.7,
            name="Pocket types (env_cluster)",
            stylize=True,
            **_TAUBIN,
        )

    fig3.plot_surface_dgrid_3js(
        dgrid,
        ISOVALUE,
        atoms,
        region_labels=regions,
        cscale=cscale,
        opacity=0.35,
        name="Regions (unmerged)",
        stylize=True,
        **_TAUBIN,
    )
    fig3.plot_surface_dgrid_3js(
        dgrid,
        ISOVALUE,
        atoms,
        region_labels=regions_pbc,
        cscale=cscale,
        opacity=0.35,
        name="Regions PBC (merged)",
        stylize=True,
        **_TAUBIN,
    )
    fig3.plot_maxima_from_grid_3js(
        dgrid,
        maxima_local,
        atoms,
        regions=regions,
        cscale=cscale,
        name=f"Maxima [:regions.max()={n_cap}]",
        maxima_info=results.get("maxima_info"),
        hide_periodic=True,
    )
    fig3.plot_maxima_from_grid_3js(
        dgrid,
        maxima_pbc,
        atoms,
        regions=regions_pbc,
        cscale=cscale,
        name="Maxima PBC (representatives)",
        maxima_pbc=maxima_pbc,
        maxima_info=results.get("maxima_info"),
        maxima_pbc_dict=results.get("maxima_pbc_dict", {}),
        hide_periodic=False,
    )
    return fig3


def build_plotly(atoms, mat_updated, results, typed, *, cscale):
    dgrid = results["dgrid"]
    regions = results["regions"]
    regions_pbc = results["regions_pbc"]
    all_maxima = results["all_maxima"]
    maxima_pbc = results["maxima_pbc"]
    n_cap = int(regions.max()) if regions.size else 0
    maxima_local = all_maxima[:n_cap] if n_cap > 0 else all_maxima[:0]
    mi = np.asarray(results["maxima_info"], dtype=int)
    local_labels = mi[: len(maxima_local), 0] if len(maxima_local) else np.array([], dtype=int)

    traces = [
        *mgr.plot_framework_as_mesh(mat_updated, stylize=True, show_rings=False),
        mgr.plot_cell(mat_atoms, color="orange", name="Unit cell"),
    ]
    if typed.regions_cluster is not None:
        traces.append(
            mgr.plot_surface_dgrid(
                dgrid,
                ISOVALUE,
                mat_atoms,
                region_labels=typed.regions_cluster,
                opacity=0.7,
                cscale=cscale,
                stylize=True,
                name="Pocket types (env_cluster)",
                **_TAUBIN,
            )
        )
    traces.extend(
        [
            mgr.plot_surface_dgrid(
                dgrid,
                ISOVALUE,
                mat_atoms,
                region_labels=regions,
                opacity=0.35,
                cscale=cscale,
                stylize=True,
                name="Regions (unmerged)",
                **_TAUBIN,
            ),
            mgr.plot_surface_dgrid(
                dgrid,
                ISOVALUE,
                mat_atoms,
                region_labels=regions_pbc,
                opacity=0.35,
                cscale=cscale,
                stylize=True,
                name="Regions PBC (merged)",
                **_TAUBIN,
            ),
            mgr.plot_maxima(
                dgrid,
                maxima_local,
                mat_atoms,
                show_spheres=False,
                marker_size=6,
                marker_color=local_labels,
                marker_colorscale=cscale,
                name=f"Maxima [:regions.max()={n_cap}]",
            ),
            mgr.plot_maxima(
                dgrid,
                maxima_pbc,
                mat_atoms,
                show_spheres=False,
                marker_size=9,
                marker_color=list(range(1, len(maxima_pbc) + 1)),
                marker_colorscale=cscale,
                name="Maxima PBC (representatives)",
            ),
        ]
    )
    return go.Figure(data=traces, layout=mgr.white_layout)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cif", type=Path, required=True)
    p.add_argument("--name", default=None, help="Title prefix for HTML")
    p.add_argument("--out-dir", type=Path, default=ROOT / "examples/viz_outputs")
    p.add_argument("--spacing", type=float, default=DEFAULT_SPACING_A)
    p.add_argument("--dgrid-backend", choices=("mic", "rtree"), default="mic")
    p.add_argument("--skip-plotly", action="store_true", help="Large cells: Three.js only")
    p.add_argument("--instanced-framework", action="store_true", help="Lighter framework in Three.js")
    args = p.parse_args()

    cif = args.cif.expanduser().resolve()
    title = args.name or cif.stem
    stem = title.replace(" ", "_").replace("(", "").replace(")", "")
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_3js = out_dir / f"{stem}_pocket_inspection_threejs.html"
    out_plotly = out_dir / f"{stem}_pocket_inspection_plotly.html"
    cscale = px.colors.qualitative.Pastel

    print(f"Loading {cif.name} (dgrid={args.dgrid_backend})...")
    atoms, mat_updated, results, typed = run_segmentation_and_typing(
        cif, args.spacing, args.dgrid_backend
    )

    print("Three.js scene...")
    fig3 = build_threejs(
        atoms,
        mat_updated,
        results,
        typed,
        cscale=cscale,
        use_mesh_framework=not args.instanced_framework,
    )
    fig3.to_html(str(out_3js))
    print(f"  {out_3js} ({out_3js.stat().st_size / 1024 / 1024:.1f} MB)")

    if not args.skip_plotly:
        print("Plotly scene (may be large)...")
        fig_p = build_plotly(atoms, mat_updated, results, typed, cscale=cscale)
        fig_p.update_layout(title=f"{title} — pocket inspection")
        fig_p.write_html(str(out_plotly))
        print(f"  {out_plotly} ({out_plotly.stat().st_size / 1024 / 1024:.1f} MB)")

    print("Done. Toggle layers in legend.")


if __name__ == "__main__":
    main()

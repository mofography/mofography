#!/usr/bin/env python3
"""Generate the small curated interactive HTML gallery assets (each < 50 MB).

Output goes to a staging dir (default /tmp/mofography-gallery-assets, folder per material);
that dir is a separate private repo (mofography/mofography-gallery-assets). Nothing is committed here.

Usage (from repo root):
    python scripts/dev/gen_gallery_assets.py                    # all materials
    python scripts/dev/gen_gallery_assets.py --material cubtc   # one of cubtc nu1000 lta zif8 irmof1
    python scripts/dev/gen_gallery_assets.py --material cubtc --out /tmp/x --step-size 3

Size levers (see MANIFEST.md): dgrid spacing 0.5 A, marching-cubes ``step_size`` per layer,
``to_html(pack_binary=True)`` (base64 float32/uint16 arrays instead of JSON lists), instanced
atoms/bonds (client-side spheres) instead of per-atom meshes, per-vertex hover text off.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import mofography as mgr  # noqa: E402
from ase.io import read  # noqa: E402
from mofography import Figure3js  # noqa: E402
from mofography.helper import trim_connections_to_unique_regions  # noqa: E402
from mofography.movies import make_aromatic_rings_trace, make_percolation_clusters_trace  # noqa: E402
from mofography.threejs_viz import qualitative_colorscale  # noqa: E402

CIFS = ROOT / "examples/TestCIFs"
MATERIALS = {
    "cubtc": ("Cu-BTC", CIFS / "P1-Cu-BTC.cif"),
    "nu1000": ("NU-1000", CIFS / "NU-1000_S1.cif"),
    "lta": ("LTA", CIFS / "LTA-96.cif"),
    "zif8": ("ZIF-8", CIFS / "ZIF-8-RASPA.cif"),
    "irmof1": ("IRMOF-1", CIFS / "P1-IRMOF-1.cif"),
}
# default viewer background per material (see MovieFigure.to_html(background=...))
BACKGROUNDS = {
    "cubtc": ("radial", "#1d3b73", "#040a1c"),   # dark navy radial
    "nu1000": ("linear", "#0d4f55", "#25282d"),  # deep teal -> charcoal, left to right
    "lta": "#2a3340",                            # solid dark slate
    "zif8": ("radial", "#ffffff", "#d3dce8"),    # soft light gradient
    "irmof1": "#ffffff",                         # plain white (deliberate)
}
SPACING = 0.5
ISO = 1.0
PASTEL = qualitative_colorscale("Pastel")
MB = 1024 * 1024


# --------------------------------------------------------------------------- pipeline
def run_pipeline(cif: Path, spacing: float = SPACING, typing: bool = True):
    atoms = read(cif)
    fw_atoms = mgr.update_with_covalent_bonds(atoms)
    dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=spacing).astype(np.float32)
    seg = mgr.do_periodic_segmentation_rust(
        dgrid, atoms, spacing=spacing, mask_thickness=0.05, h=0.2, min_distance=1.5,
        threshold_abs=1.5, extra_padding=5.0, compute_connections=True,
        watershed_backend="rust_det", return_dgrid=True, disable_tqdm=True)
    typed = None
    if typing:
        typed = mgr.assign_pocket_types(seg, atoms, method="env_cluster")
        if not typed.ok:
            print("  pocket typing failed:", typed.error)
            typed = None
    return atoms, fw_atoms, seg, typed


# --------------------------------------------------------------------------- layers
def add_framework(fig, fw_atoms, atom_scale=0.45, rings=True, extra_styles=True):
    """Instanced atoms + bonds (+ aromatic rings): client-side spheres, tiny payload.

    ``extra_styles`` also embeds instanced stick-only and vdW-sphere variants (ids
    ``framework_stick_*`` / ``framework_spacefill_atoms``); the viewer's "Framework style" dropdown
    (ball-and-stick / stick / space-filling) switches between them and only the selected one is shown.
    """
    fig.plot_framework_instanced_3js(fw_atoms, atom_scale=atom_scale, name="Framework")
    fig.update_traces({"visible": True}, selector={"id": "framework_atoms"})
    if extra_styles:
        fig.plot_framework_instanced_3js(fw_atoms, trace_id="framework_stick", representation="stick",
                                         name="Framework (stick)")
        fig.plot_framework_instanced_3js(fw_atoms, trace_id="framework_spacefill", representation="space_filling",
                                         name="Framework (vdW spheres)")
    if rings:
        try:
            fig.add_traces(make_aromatic_rings_trace(fw_atoms))
        except Exception as exc:  # materials without aromatic rings (e.g. LTA)
            print("  no aromatic rings:", exc)


def add_surface(fig, seg, atoms, labels, name, trace_id, step, opacity=1.0, cscale=PASTEL, visible=True):
    fig.plot_surface(seg["dgrid"], labels, atoms, isovalue=ISO, name=name, trace_id=trace_id,
                     opacity=opacity, cscale=cscale, step_size=step, stylize=True,
                     taubin_smooth=True, taubin_iterations=3, taubin_lambda=0.2175,
                     taubin_nu=-0.2175, hover_label="")  # hover_label="" -> no per-vertex hover strings
    fig.update_traces({"visible": visible}, selector={"id": trace_id})


def set_legend_label(fig, trace_id, label):
    """Override the generic schema legend label of the entry containing ``trace_id`` (viewer ``legend_label``)."""
    for t in fig.data:
        if t.get("id") == trace_id:
            t["legend_label"] = label
            t["name"] = label


def add_percolation(fig, seg, hist, atoms, probe, step, visible=True, opacity=0.9):
    """Percolated-cluster mesh(es) at one probe radius, as a separately toggleable layer."""
    traces = make_percolation_clusters_trace(hist, seg["dgrid"], seg["regions"], atoms, probe_radius=probe,
                                             isovalue=ISO, step_size=step, cluster_opacity=opacity)
    for t in traces:
        t["id"] = f"{t['id']}_probe{probe:g}"
        t["name"] = f"Percolated cluster, probe r={probe:g} A"
        t["visible"] = visible
    fig.add_traces(traces)
    for t in traces:
        fig.groups[t["id"]] = {"traces": [t["id"]], "label": t["name"]}


def add_pore_graph(fig, seg, atoms, with_periodic=True):
    conns = trim_connections_to_unique_regions(seg["connections"], seg["unique_regions"])
    fig.plot_connections_as_it_is_3js(conns, seg["all_maxima"], seg["regions"].shape, atoms,
                                      plot_periodic=with_periodic, name="Connections", trace_id="connections")
    spacing = float(np.mean(np.linalg.norm(atoms.get_cell().array.T, axis=0) / np.asarray(seg["regions"].shape)))
    # plot_windows_from_results (not plot_windows) reclassifies windows far outside the cell as "periodic",
    # so the default view shows only windows inside the central cell; periodic ones are an off-by-default toggle.
    fig.plot_windows_from_results(conns, dgrid=seg["dgrid"], atoms=atoms, spacing=spacing, hide_periodic=True,
                                  name="Windows", marker_color="red", marker_size=8)


def finish(fig, atoms, path: Path, hide: list[str] = (), camera=None, background=None):
    for tid in hide:
        fig.update_traces({"visible": False}, selector={"id": tid})
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.to_html(path, pack_binary=True, background=background)
    size = path.stat().st_size / MB
    print(f"  {path.relative_to(path.parents[1])}: {size:.1f} MB")
    return size


def layer_names(fig):
    return [t.get("name") for t in fig.data]


def trace_report(fig):
    rows = []
    for t in fig.data:
        n = len(json.dumps(t.get("data", {}))) / MB
        rows.append((t.get("id"), round(n, 2)))
    return rows


# --------------------------------------------------------------------------- scenes
def build_cubtc(out: Path, step: int = 1, perc_step: int = 1, probe: float = 1.2):
    name, cif = MATERIALS["cubtc"]
    atoms, fw, seg, typed = run_pipeline(cif)
    regions, regions_pbc = seg["regions"], seg["regions_pbc"]
    fig = Figure3js(layout="white")
    add_framework(fig, fw)
    fig.plot_cell(atoms, color="orange", name="Unit cell")
    fig.plot_maxima(seg, atoms, name="Local maxima (PBC)", trace_id="maxima", cscale=PASTEL, hide_periodic=False)
    add_surface(fig, seg, atoms, np.ones_like(regions), "Pore surface (distance grid, isovalue 1.0 A)", "surface", step,
                opacity=0.5, cscale=[[0.6, 0.75, 0.95]])
    set_legend_label(fig, "surface", "Pore surface (distance grid, isovalue 1.0 A)")
    add_surface(fig, seg, atoms, regions, "Regions (no PBC)", "regions_nopbc", step, visible=False)
    add_surface(fig, seg, atoms, regions_pbc, "Regions (PBC-merged)", "surface_regions_pbc", step, visible=False)
    set_legend_label(fig, "surface_regions_pbc", "Regions (PBC-merged)")
    if typed is not None:
        add_surface(fig, seg, atoms, typed.regions_cluster, "Pocket types (env_cluster)", "surface_regions_cluster", step,
                    visible=False)
    add_pore_graph(fig, seg, atoms)
    # percolated cluster at a probe radius
    hist, _, _ = mgr.compute_percolation_history_winding(seg, step=0.1)
    add_percolation(fig, seg, hist, atoms, probe, perc_step, visible=False)
    set_legend_label(fig, "surface_regions_cluster", "Pocket types (env_cluster)")
    keep = {"framework_atoms", "framework_bonds", "framework_stick_atoms", "framework_stick_bonds",
            "framework_spacefill_atoms", "cell", "surface", "aromatic_rings", "aromatic_highlights",
            "aromatic_highlights_hex", "aromatic_highlights_pent"}
    hide = [t["id"] for t in fig.data if t["id"] not in keep and t.get("visible", True)]
    path = out / "cubtc" / "cubtc_all_in_one_threejs.html"
    size = finish(fig, atoms, path, hide=hide, background=BACKGROUNDS["cubtc"])
    print("  trace payload MB (json):", trace_report(fig))
    return fig, seg, atoms, typed, size


def build_simple(key: str, out: Path, step: int = 1, region_step: int | None = None):
    """Pockets + windows + framework (+ pore graph) toggles."""
    region_step = step if region_step is None else region_step
    name, cif = MATERIALS[key]
    atoms, fw, seg, typed = run_pipeline(cif)
    fig = Figure3js(layout="white")
    add_framework(fig, fw)
    fig.plot_cell(atoms, color="orange", name="Unit cell")
    labels = typed.regions_cluster if typed is not None else seg["regions_pbc"]
    add_surface(fig, seg, atoms, labels, "Pocket types (env_cluster)" if typed is not None else "Pockets (PBC)",
                "surface_regions_cluster", step, opacity=0.9)
    set_legend_label(fig, "surface_regions_cluster", "Pocket types (env_cluster)" if typed is not None else "Pockets (PBC)")
    # extra region colourings, off by default (separate legend toggles)
    add_surface(fig, seg, atoms, seg["regions"], "Regions (no PBC)", "regions_nopbc", region_step, visible=False)
    add_surface(fig, seg, atoms, seg["regions_pbc"], "Regions (PBC-merged)", "surface_regions_pbc", region_step,
                visible=False)
    set_legend_label(fig, "surface_regions_pbc", "Regions (PBC-merged)")
    fig.plot_maxima(seg, atoms, name="Local maxima", trace_id="maxima", cscale=PASTEL, hide_periodic=True)
    add_pore_graph(fig, seg, atoms, with_periodic=False)
    hide = ["maxima_internal", "maxima_outside", "maxima_pbc"]
    path = out / key / f"{key}_pockets_windows_threejs.html"
    size = finish(fig, atoms, path, hide=[h for h in hide if any(t["id"] == h for t in fig.data)],
                  background=BACKGROUNDS[key])
    print("  trace payload MB (json):", trace_report(fig))
    return fig, seg, atoms, typed, size


def build_nu1000_percolation(out: Path, step: int = 1, probes=(1.0, 2.0, 3.0, 4.0)):
    """Percolation-history demo: one toggleable cluster layer per probe radius (3.0 A shown first)."""
    name, cif = MATERIALS["nu1000"]
    atoms, fw, seg, _ = run_pipeline(cif, typing=False)
    hist, _, _ = mgr.compute_percolation_history_winding(seg, step=0.1)
    analysis = dict(seg, percolation_history=hist, isovalue=ISO)
    fig = Figure3js(layout="white")
    add_framework(fig, fw)
    fig.plot_cell(atoms, color="orange", name="Unit cell")
    for r in probes:
        add_percolation(fig, seg, hist, atoms, r, step, visible=(r == 3.0))
    try:
        fig.plot_percolation_windows(analysis, atoms, probe_radius=3.0)
    except Exception as exc:
        print("  percolation windows skipped:", exc)
    path = out / "nu1000" / "nu1000_percolation_history_threejs.html"
    finish(fig, atoms, path)
    return path


def build_cubtc_sa_vv(out: Path, probes=(0.5, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0), n_samples=20000):
    """Standalone Plotly: surface area & void volume (total vs accessible) vs probe radius, + per-pocket SA."""
    import plotly.graph_objects as go
    name, cif = MATERIALS["cubtc"]
    atoms, fw, seg, _ = run_pipeline(cif, typing=False)
    sa, vv = [], []
    for r in probes:
        sa.append(mgr.compute_surface_area_mc_winding(seg, atoms, probe_radius=r, n_samples=n_samples,
                                                      shell_thickness=0.5))
        vv.append(mgr.compute_void_volume_from_results_winding(seg, atoms, probe_radius=r))
    files = []
    for fname, title, key_t, key_a, unit, rows in (
        ("cubtc_surface_area_plotly.html", "Cu-BTC surface area vs probe radius (MC shell)",
         "total_surface_area", "accessible_surface_area", "m<sup>2</sup>/g", sa),
        ("cubtc_void_volume_plotly.html", "Cu-BTC void volume vs probe radius",
         "total_void_volume", "accessible_void_volume", "cm<sup>3</sup>/g", vv),
    ):
        f = go.Figure()
        f.add_scatter(x=list(probes), y=[r[key_t] for r in rows], mode="lines+markers", name="total")
        f.add_scatter(x=list(probes), y=[r[key_a] for r in rows], mode="lines+markers", name="accessible")
        f.update_layout(title=title, xaxis_title="probe radius (A)", yaxis_title=unit, template="plotly_white")
        path = out / "cubtc" / fname
        f.write_html(path, include_plotlyjs="cdn")
        print(f"  {fname}: {path.stat().st_size / MB:.2f} MB")
        files.append(path)
    return files


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--material", nargs="*", default=list(MATERIALS))
    ap.add_argument("--out", type=Path, default=Path("/tmp/mofography-gallery-assets"))
    ap.add_argument("--step-size", type=int, default=1, help="marching-cubes step_size for surface layers")
    ap.add_argument("--sa-vv", action="store_true", help="also build Cu-BTC Plotly SA / void-volume plots")
    ap.add_argument("--perc", action="store_true", help="also build NU-1000 percolation demo")
    a = ap.parse_args()
    for m in a.material:
        t0 = time.time()
        print(f"== {m}")
        if m == "cubtc":
            build_cubtc(a.out, step=a.step_size, perc_step=a.step_size)
        else:
            build_simple(m, a.out, step=a.step_size)
        print(f"  ({time.time() - t0:.0f}s)")
    if a.perc:
        build_nu1000_percolation(a.out, step=a.step_size)
    if a.sa_vv:
        build_cubtc_sa_vv(a.out)

#!/usr/bin/env python3
"""Pocket-type SI report (4-panel PDF + Three.js HTML). Batch-friendly."""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
from matplotlib.colors import to_hex, to_rgb
from ase.io import read
from matplotlib.gridspec import GridSpec

import mofography as mgr
from mofography import Figure3JS
from mofography.helper import trim_connections_to_unique_regions
from mofography.pocket_framework_env import (
    PocketFingerprint,
    _cosine,
    fingerprint_cosine_similarity,
)
from mofography.pocket_typing import (
    DEFAULT_MASK_THICKNESS_A,
    DEFAULT_PROBE_RADIUS_PSD_A,
    DEFAULT_SPACING_A,
    assign_pocket_types_from_segmentation,
)

REPO = Path(__file__).resolve().parents[1]
RASPA_CIF = Path.home() / "repos/RASPA2/structures/mofs/cif"
TEST_CIF = REPO / "examples/TestCIFs"

TYPE_COLORS = ["#4C72B0", "#55A868", "#DD8452", "#C44E52", "#8172B2", "#937860"]  # legacy fallback
CSCALE_NAME = "Pastel"  # must match build_pocket_types_figure3js / plot_surface_dgrid


def pocket_type_palette(name: str = CSCALE_NAME) -> List[str]:
    """Plotly qualitative palette used for typed surfaces in panel (a)."""
    return list(getattr(px.colors.qualitative, name))


def _mpl_color(color: str) -> str:
    """Plotly Pastel entries are often ``rgb(r,g,b)``; matplotlib wants hex."""
    s = str(color).strip()
    if s.startswith("rgb"):
        nums = [float(x) for x in re.findall(r"[\d.]+", s)[:3]]
        if nums and max(nums) > 1.0:
            nums = [n / 255.0 for n in nums]
        return to_hex(nums[:3])
    return to_hex(to_rgb(s))


def color_for_pocket_type(type_id: int, palette: Optional[Sequence[str]] = None) -> str:
    """Color for env_cluster type_id, matching Three.js surface coloring.

    ``regions_cluster`` paints label ``type_id + 1``; surface code assigns
    ``palette[i]`` to the i-th sorted positive label → type *k* → ``palette[k]``.
    """
    pal = list(palette) if palette is not None else pocket_type_palette()
    return _mpl_color(pal[int(type_id) % len(pal)])
SEG_KW = dict(
    mask_thickness=DEFAULT_MASK_THICKNESS_A,
    h=0.2,
    min_distance=1.5,
    threshold_abs=1.5,
    extra_padding=5.0,
    return_dgrid=True,
)
TAUBIN = dict(
    taubin_smooth=True,
    taubin_iterations=3,
    taubin_lambda=0.2175,
    taubin_nu=-0.2175,
    stepsize=1,
)
ISOVALUE = 1.0
VIEW_DIR_C = (0.0, 0.0, 1.0)
# [001] orthographic: offset along c by factor × max(la, lb) so ab-plane fits frustum.
PANEL_A_CAMERA_Z_FACTOR = 1.8
# Panel (a) overlay visibility — cell + local maxima must read clearly in PDF PNG export.
PANEL_A_CELL_KW = dict(color=0xFF6600, linewidth=5)
PANEL_A_MAXIMA_SIZE_SCALE = 2.8
PANEL_A_MAXIMA_STYLE_SIZE = 18.0

# slug -> (cif_path, display_name, dgrid_backend)
BATCH_MATERIALS: Dict[str, tuple] = {
    "cubtc": (TEST_CIF / "P1-Cu-BTC.cif", "Cu-BTC", "mic"),
    "irmof1": (TEST_CIF / "P1-IRMOF-1.cif", "IRMOF-1", "mic"),
    "uio66": (TEST_CIF / "UiO-66-RASPA.cif", "UiO-66", "mic"),
    "dut75": (TEST_CIF / "DUT-75_uncapped.cif", "DUT-75", "mic"),
    "lta96": (TEST_CIF / "LTA-96.cif", "LTA (2x2x2)", "mic"),
    "zif20": (RASPA_CIF / "ZIF-20.cif", "ZIF-20", "rtree"),
    "nu1000": (TEST_CIF / "NU-1000_S1.cif", "NU-1000", "mic"),
    "irmof10": (RASPA_CIF / "IRMOF-10.cif", "IRMOF-10", "mic"),
    "umcm1": (RASPA_CIF / "UMCM-1.cif", "UMCM-1", "mic"),
}


def _expand(path) -> Path:
    return Path(path).expanduser().resolve()


def display_name_from_stem(stem: str) -> str:
    """Pretty title: drop P1-, _RASPA, uncapped; LTA96 -> LTA (2x2x2)."""
    s = stem.replace("_", "-")
    if re.match(r"(?i)^P1-", s):
        s = s[3:]
    s = re.sub(r"(?i)-RASPA$", "", s)
    s = re.sub(r"(?i)_?uncapped$", "", s)
    s = re.sub(r"(?i)_?experimental.*$", "", s)
    if re.match(r"(?i)^LTA-?96$", s) or s.upper() == "LTA96":
        return "LTA (2x2x2)"
    if re.match(r"(?i)^NU-1000[_-]?S1$", s):
        return "NU-1000"
    return s


def _pair_key(ti, tj):
    a, b = min(ti, tj), max(ti, tj)
    return f"{a}_vs_{b}"


def _lookup_pairwise(pairwise, ti, tj, default=0.0):
    return pairwise.get(_pair_key(ti, tj), default)


def channel_cosine(a: PocketFingerprint, b: PocketFingerprint):
    out = {}
    if a.psd_vec is not None and b.psd_vec is not None and len(a.psd_vec) and len(b.psd_vec):
        out["psd"] = float(_cosine(a.psd_vec, b.psd_vec))
    if len(a.global_sh_power) and len(b.global_sh_power):
        out["sh_global"] = float(_cosine(a.global_sh_power, b.global_sh_power))
    if len(a.global_radial) and len(b.global_radial):
        out["radial"] = float(_cosine(a.global_radial, b.global_radial))
    if (
        a.conn_sh_power is not None
        and b.conn_sh_power is not None
        and a.n_connections >= 3
        and b.n_connections >= 3
    ):
        out["connectivity"] = float(_cosine(a.conn_sh_power, b.conn_sh_power))
    return out


def build_dgrid(atoms, spacing: float, backend: str = "mic") -> np.ndarray:
    if backend == "rtree":
        return mgr.dgrid_from_atoms_rtree_rust(atoms, spacing=spacing).astype(np.float32)
    return mgr.dgrid_from_atoms_rust(atoms, spacing=spacing).astype(np.float32)


def run_pipeline(
    cif_path,
    spacing=DEFAULT_SPACING_A,
    probe_radius_psd=DEFAULT_PROBE_RADIUS_PSD_A,
    dgrid_backend: str = "mic",
):
    atoms = read(cif_path)
    mat_updated = mgr.update_with_covalent_bonds(atoms)
    dgrid = build_dgrid(atoms, spacing, dgrid_backend)
    results = mgr.do_periodic_segmentation_rust(dgrid, atoms, spacing=spacing, **SEG_KW)
    if results.get("pocket_psds") is None:
        from mofography.pocket_typing import psd_bins

        regions_pbc = results.get("regions_pbc", results["regions"])
        results["pocket_psds"] = mgr.compute_psds_pocketwise_rust(
            results["dgrid"],
            regions_pbc,
            atoms,
            num_samples_per_region=1000,
            bins=psd_bins(results["dgrid"]),
            probe_center_min_radius=probe_radius_psd,
        )
    typed = assign_pocket_types_from_segmentation(
        results,
        atoms,
        method="env_cluster",
        probe_radius_psd=probe_radius_psd,
        cif_path=cif_path,
    )
    if not typed.ok:
        raise RuntimeError(typed.error)
    return atoms, mat_updated, results, typed, {
        "cif": str(cif_path),
        "n_types": typed.n_types,
        "n_pockets": typed.n_pockets,
        "audit": typed.audit,
        "dgrid_backend": dgrid_backend,
    }


def type_groups(rep_to_type):
    g = defaultdict(list)
    for rep, tid in rep_to_type.items():
        g[int(tid)].append(int(rep))
    return {k: sorted(v) for k, v in sorted(g.items())}


def per_type_reps(rep_to_type, fps):
    groups = type_groups(rep_to_type)
    return {tid: min(reps, key=lambda r: fps[r].max_probe_radius_a) for tid, reps in groups.items()}


def fingerprint_tables(rep_to_type, fps, pockets):
    groups = type_groups(rep_to_type)
    reps_by_type = per_type_reps(rep_to_type, fps)
    pocket_by_rep = {int(p["id"]): p for p in pockets}
    per_type = {}
    for tid, rep in reps_by_type.items():
        fp, p = fps[rep], pocket_by_rep[rep]
        per_type[tid] = {
            "representative_rep": rep,
            "n_pockets": len(groups[tid]),
            "d1_A": round(float(p.get("d1_A", 0)), 2),
            "max_probe_radius_A": round(float(fp.max_probe_radius_a), 3),
            "n_connections": int(fp.n_connections),
        }
    types = sorted(per_type)
    pairwise, channel = {}, {}
    for i, ti in enumerate(types):
        for tj in types[i:]:
            key = _pair_key(ti, tj)
            pairwise[key] = float(
                fingerprint_cosine_similarity(fps[reps_by_type[ti]], fps[reps_by_type[tj]])
            )
            if ti != tj:
                channel[key] = channel_cosine(fps[reps_by_type[ti]], fps[reps_by_type[tj]])
    return {
        "n_types": len(types),
        "type_groups": {str(k): groups[k] for k in groups},
        "per_type_info": {str(k): per_type[k] for k in per_type},
        "pairwise_cosine": pairwise,
        "channel_cosine": channel,
    }


def plot_four_panel(typed, tables, out_pdf, panel_a_image=None, palette: Optional[Sequence[str]] = None):
    if palette is None:
        palette = pocket_type_palette()
    plt.rcParams["font.size"] = 15
    fps, pockets, rep_to_type = typed.fingerprints, typed.pockets, typed.rep_to_type
    reps_by_type = per_type_reps(rep_to_type, fps)
    types = sorted(reps_by_type)
    n_types = len(types)
    fig = plt.figure(figsize=(11, 9))
    gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.28)
    ax_a = fig.add_subplot(gs[0, 0])
    if panel_a_image and Path(panel_a_image).is_file():
        ax_a.imshow(plt.imread(panel_a_image))
        ax_a.axis("off")
    else:
        ax_a.text(
            0.5,
            0.5,
            f"{n_types} types\nExport PNG from HTML",
            ha="center",
            va="center",
            transform=ax_a.transAxes,
        )
        ax_a.axis("off")
    ax_a.set_title("(a) Pocket types (3D)", fontweight="bold")
    ax_b = fig.add_subplot(gs[0, 1])
    for i, tid in enumerate(types):
        rep = reps_by_type[tid]
        p = next(x for x in pockets if int(x["id"]) == rep)
        x_psd, y_psd = np.asarray(p.get("psd_x", [])), np.asarray(p.get("psd_vec", []))
        if x_psd.size and y_psd.size:
            ax_b.plot(
                x_psd,
                y_psd,
                color=color_for_pocket_type(tid, palette),
                lw=2,
                label=f"Type {tid}",
            )
    ax_b.set_xlabel("Pore diameter (Å)")
    ax_b.set_ylabel("Normalized PSD")
    ax_b.legend(fontsize=12)
    ax_b.set_title("(b) Pocket PSDs", fontweight="bold")
    ax_c = fig.add_subplot(gs[1, 0])
    bar_w, x_pos = 0.35, np.arange(7)
    for i, tid in enumerate(types):
        fp = fps[reps_by_type[tid]]
        fw = np.asarray(fp.global_sh_power, dtype=float)
        conn = np.asarray(
            fp.conn_sh_power if fp.conn_sh_power is not None else np.zeros(7), dtype=float
        )
        if fw.max() > 0:
            fw /= fw.max()
        if conn.max() > 0:
            conn /= conn.max()
        off = (i - (n_types - 1) / 2) * bar_w
        tc = color_for_pocket_type(tid, palette)
        ax_c.bar(x_pos + off, fw, bar_w * 0.9, color=tc, alpha=0.85)
        ax_c.bar(
            x_pos + off,
            conn,
            bar_w * 0.9,
            color=tc,
            alpha=0.35,
            hatch="///",
            edgecolor="gray",
        )
    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels([f"ℓ={l}" for l in range(7)], fontsize=12)
    ax_c.set_ylabel("Norm. SH power")
    ax_c.set_title("(c) Framework (solid) & connectivity (hatched) SH", fontsize=14, fontweight="bold")
    ax_d = fig.add_subplot(gs[1, 1])
    mat = np.eye(n_types)
    for i, ti in enumerate(types):
        for j, tj in enumerate(types):
            if i != j:
                mat[i, j] = _lookup_pairwise(tables["pairwise_cosine"], ti, tj)
    im = ax_d.imshow(mat, vmin=0, vmax=1, cmap="RdYlGn")
    ax_d.set_xticks(range(n_types))
    ax_d.set_yticks(range(n_types))
    tick_colors = [color_for_pocket_type(t, palette) for t in types]
    ax_d.set_xticklabels([f"T{t}" for t in types])
    ax_d.set_yticklabels([f"T{t}" for t in types])
    for lbl, c in zip(ax_d.get_xticklabels(), tick_colors):
        lbl.set_color(c)
    for lbl, c in zip(ax_d.get_yticklabels(), tick_colors):
        lbl.set_color(c)
    for i in range(n_types):
        for j in range(n_types):
            ax_d.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=13)
    plt.colorbar(im, ax=ax_d, fraction=0.046)
    ax_d.set_title("(d) Cosine similarity", fontweight="bold")
    fig.suptitle(
        f"Pocket types — {tables.get('material_name', '')}",
        fontweight="bold",
        y=0.98,
    )
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _camera_eye_center(atoms):
    cell = atoms.get_cell().array  # rows = a, b, c (ASE convention)
    lengths = np.linalg.norm(cell, axis=1)
    center = 0.5 * cell.sum(axis=0)  # Cartesian center of unit cell (frac 0.5, 0.5, 0.5)
    in_plane = float(max(lengths[0], lengths[1]))
    c_offset = in_plane * PANEL_A_CAMERA_Z_FACTOR
    eye = center + np.asarray(VIEW_DIR_C, dtype=float) * c_offset
    return tuple(eye), tuple(center)


def _apply_panel_a_scene(fig3, atoms):
    """Orthographic [001] camera + explicit lighting for panel (a) export."""
    eye, center = _camera_eye_center(atoms)
    fig3.set_camera(eye=eye, center=center, up=VIEW_DIR_C, projection="orthographic")
    fig3.set_lights(
        ambient=0.015,
        hemisphere_intensity=0.100,
        hemisphere_sky="#ddeeff",
        hemisphere_ground="#3b3b3b",
        dir_key_intensity=0.576,
        dir_key_position=(5, 10, 7.5),
        dir_fill_intensity=0.070,
        dir_fill_position=(-5, -9, -7.5),
        ssao=False,
        explicit=True,
    )
    return eye, center


def _boost_panel_a_trace_visibility(fig3) -> None:
    """Thicken unit cell + enlarge local-maxima markers for panel (a) PNG."""
    for trace in getattr(fig3._fig, "data", []):
        tid = str(trace.get("id", ""))
        if trace.get("trace_type") == "unit_cell_lines" or tid == "cell":
            st = trace.setdefault("style", {})
            st.update(PANEL_A_CELL_KW)
        if trace.get("trace_type") == "points" and "maxima" in tid:
            data = trace.setdefault("data", {})
            sizes = data.get("sizes") or []
            if sizes:
                data["sizes"] = [float(s) * PANEL_A_MAXIMA_SIZE_SCALE for s in sizes]
            st = trace.setdefault("style", {})
            st["size"] = PANEL_A_MAXIMA_STYLE_SIZE


def _add_panel_a_cell_and_maxima(fig3, atoms, results, cscale):
    dgrid, regions = results["dgrid"], results["regions"]
    maxima_pbc = results["maxima_pbc"]
    fig3.plot_cell(atoms, name="cell", **PANEL_A_CELL_KW)
    fig3.plot_maxima_from_grid_3js(
        dgrid,
        maxima_pbc,
        atoms,
        regions=regions,
        cscale=cscale,
        name="Local maxima",
        hide_periodic=False,
    )
    _boost_panel_a_trace_visibility(fig3)


def should_use_lite_panel_a(results, typed) -> bool:
    """Large cells crash Playwright on full ball-and-stick + graph overlays."""
    dgrid = np.asarray(results["dgrid"])
    n_vox = int(np.prod(dgrid.shape))
    return int(typed.n_pockets) > 40 or n_vox > 150_000 or len(results.get("connections", [])) > 500


def build_panel_a_figure3js(atoms, mat_updated, results, typed, *, lite: bool = False):
    """3D scene for panel (a) PNG export — lite omits graph overlays for huge cells."""
    dgrid, regions = results["dgrid"], results["regions"]
    cscale = pocket_type_palette(CSCALE_NAME)
    if lite:
        fig3 = (
            Figure3JS(layout="white", stylize=True)
            .plot_framework_instanced_3js(mat_updated, atom_scale=0.45, name="Framework")
        )
        _add_panel_a_cell_and_maxima(fig3, atoms, results, cscale)
        stepsize = 2
    else:
        connections, all_maxima = results["connections"], results["all_maxima"]
        maxima_pbc, unique_regions = results["maxima_pbc"], results["unique_regions"]
        connections_trimmed = trim_connections_to_unique_regions(connections, unique_regions)
        spacing_win = float(
            np.mean(np.linalg.norm(atoms.get_cell().array.T, axis=0)
                    / np.maximum(np.asarray(dgrid.shape, dtype=float), 1.0))
        )
        fig3 = (
            Figure3JS(layout="white", stylize=True)
            .plot_framework(mat_updated, representation="ball_and_stick", name="Framework")
        )
        _add_panel_a_cell_and_maxima(fig3, atoms, results, cscale)
        stepsize = 1
    if typed.regions_cluster is not None:
        surf_kw = dict(TAUBIN)
        surf_kw["stepsize"] = stepsize
        fig3.plot_surface_dgrid(
            dgrid,
            ISOVALUE,
            atoms,
            region_labels=typed.regions_cluster,
            cscale=cscale,
            name="Pocket types",
            stylize=True,
            **surf_kw,
        )
    if not lite:
        fig3.plot_connections_as_it_is_3js(
            connections_trimmed,
            all_maxima,
            regions.shape,
            atoms,
            plot_periodic=False,
            name="Internal connections",
        )
        fig3.plot_windows(
            connections_trimmed,
            atoms=atoms,
            dgrid=dgrid,
            maxima=all_maxima,
            grid_shape=regions.shape,
            spacing=spacing_win,
            hide_periodic=True,
            name="Internal windows",
            marker_color="red",
            marker_size=8,
        )
    _apply_panel_a_scene(fig3, atoms)
    return fig3


def build_pocket_types_figure3js(atoms, mat_updated, results, typed):
    dgrid, regions = results["dgrid"], results["regions"]
    connections, all_maxima = results["connections"], results["all_maxima"]
    maxima_pbc, unique_regions = results["maxima_pbc"], results["unique_regions"]
    connections_trimmed = trim_connections_to_unique_regions(connections, unique_regions)
    cscale = pocket_type_palette(CSCALE_NAME)
    spacing_win = float(
        np.mean(np.linalg.norm(atoms.get_cell().array.T, axis=0)
                / np.maximum(np.asarray(dgrid.shape, dtype=float), 1.0))
    )
    fig3 = (
        Figure3JS(layout="white", stylize=True)
        .plot_framework(mat_updated, representation="ball_and_stick", name="Framework")
    )
    _add_panel_a_cell_and_maxima(fig3, atoms, results, cscale)
    if typed.regions_cluster is not None:
        fig3.plot_surface_dgrid(
            dgrid,
            ISOVALUE,
            atoms,
            region_labels=typed.regions_cluster,
            cscale=cscale,
            name="Pocket types",
            stylize=True,
            **TAUBIN,
        )
    fig3.plot_connections_as_it_is_3js(
        connections_trimmed,
        all_maxima,
        regions.shape,
        atoms,
        plot_periodic=False,
        name="Internal connections",
    )
    fig3.plot_windows(
        connections_trimmed,
        atoms=atoms,
        dgrid=dgrid,
        maxima=all_maxima,
        grid_shape=regions.shape,
        spacing=spacing_win,
        hide_periodic=True,
        name="Internal windows",
        marker_color="red",
        marker_size=8,
    )
    _apply_panel_a_scene(fig3, atoms)
    return fig3


_SET_PANEL_A_CAMERA_JS = """
({eye, center, up}) => {
  camera.position.set(eye[0], eye[1], eye[2]);
  controls.target.set(center[0], center[1], center[2]);
  camera.up.set(up[0], up[1], up[2]);
  camera.lookAt(controls.target);
  controls.update();
  [perspectiveCam, orthoCam].forEach(function(cam) {
    cam.position.copy(camera.position);
    cam.up.copy(camera.up);
    cam.lookAt(controls.target);
  });
  if (typeof updateOrthoFrustum === 'function') updateOrthoFrustum();
}
"""


def export_panel_a_png(fig3, out_png, *, eye, center, wait_ms=2500, up=VIEW_DIR_C):
    import base64
    import tempfile

    out = _expand(out_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix.lower() != ".png":
        out = out.with_suffix(".png")
    tmp_dir = Path(tempfile.mkdtemp(prefix="mofography_panel_a_"))
    html_path = tmp_dir / "scene.html"
    fig3.to_html(html_path)
    from playwright.sync_api import sync_playwright

    w, h = 1280, 960
    cam = {
        "eye": [float(eye[0]), float(eye[1]), float(eye[2])],
        "center": [float(center[0]), float(center[1]), float(center[2])],
        "up": [float(up[0]), float(up[1]), float(up[2])],
    }
    with sync_playwright() as p:
        browser = p.chromium.launch(
            args=["--disable-dev-shm-usage", "--disable-gpu", "--no-sandbox"],
        )
        page = browser.new_page(viewport={"width": w, "height": h})
        page.goto(html_path.as_uri(), wait_until="domcontentloaded", timeout=120_000)
        page.wait_for_timeout(wait_ms)
        page.evaluate(_SET_PANEL_A_CAMERA_JS, cam)
        page.wait_for_timeout(200)
        data_url = page.evaluate("(mult) => window.__mofFigureExportPng(mult)", 2)
        browser.close()
    b64 = str(data_url).split(",", 1)[1]
    out.write_bytes(base64.b64decode(b64))
    return out


def emit_org_si(name, meta, tables, out_org, pdf_name, html_name):
    lines = [
        f"* Pocket type report — {name}",
        "",
        f"CIF: ={meta['cif']}= · Types: *{tables['n_types']}* · Pockets: *{meta['n_pockets']}*",
        f"PDF: [[file:{pdf_name}]] · HTML: [[file:{html_name}]]",
        "",
        "** Type summary",
        "",
        "| Type | Rep | # | d₁ (Å) | max_r (Å) | n_conn |",
        "|------+-----+---+--------+-----------+--------+",
    ]
    for tid in sorted(tables["per_type_info"], key=int):
        r = tables["per_type_info"][tid]
        lines.append(
            f"| {tid} | {r['representative_rep']} | {r['n_pockets']} | {r['d1_A']:.2f} | {r['max_probe_radius_A']:.3f} | {r['n_connections']} |"
        )
    if tables["channel_cosine"]:
        lines += [
            "",
            "** Channel cosine",
            "",
            "| Pair | PSD | SH | Radial | Conn | Overall |",
            "|------+-----+----+--------+------+---------|",
        ]
        types = sorted({int(k.split("_")[0]) for k in tables["pairwise_cosine"] if "_vs_" in k})
        for i, ti in enumerate(types):
            for tj in types[i + 1 :]:
                ch = tables["channel_cosine"].get(_pair_key(ti, tj), {})
                ov = _lookup_pairwise(tables["pairwise_cosine"], ti, tj)
                lines.append(
                    f"| {ti} vs {tj} | {ch.get('psd', float('nan')):.2f} | {ch.get('sh_global', float('nan')):.2f} | {ch.get('radial', float('nan')):.2f} | {ch.get('connectivity', float('nan')):.2f} | {ov:.2f} |"
                )
    lines += [
        "",
        "** Panel (c)",
        "",
        "Bar height = SH power P_ℓ (ℓ=0…6). Solid = framework; hatched = connectivity.",
        "",
    ]
    out_org.write_text("\n".join(lines), encoding="utf-8")


def file_stem(display_name: str) -> str:
    """Filesystem-safe stem from display name."""
    s = display_name.replace(" ", "_").replace("(", "").replace(")", "").replace("x", "x")
    s = s.replace("/", "-")
    return re.sub(r"[^\w.-]+", "_", s)


def generate_report(
    cif_path: Path,
    out_dir: Path,
    *,
    display_name: Optional[str] = None,
    spacing: float = DEFAULT_SPACING_A,
    probe_radius_psd: float = DEFAULT_PROBE_RADIUS_PSD_A,
    dgrid_backend: str = "mic",
    no_html: bool = False,
    no_auto_panel_a: bool = False,
    panel_a_wait_ms: int = 2500,
    write_json: bool = False,
    panel_a_png_only: bool = False,
) -> Dict[str, Any]:
    cif_path = _expand(cif_path)
    if not cif_path.is_file():
        raise FileNotFoundError(cif_path)
    title = display_name or display_name_from_stem(cif_path.stem)
    stem = file_stem(title)
    out_dir = _expand(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Run {title} ({cif_path.name}) dgrid={dgrid_backend} …")
    atoms, mat_updated, results, typed, meta = run_pipeline(
        cif_path,
        spacing=spacing,
        probe_radius_psd=probe_radius_psd,
        dgrid_backend=dgrid_backend,
    )
    meta["n_pockets"] = typed.n_pockets
    tables = fingerprint_tables(typed.rep_to_type, typed.fingerprints, typed.pockets)
    tables["material_name"] = title

    pdf_path = out_dir / f"{stem}_pocket_types_report.pdf"
    html_path = out_dir / f"{stem}_pocket_types_panel_a.html"
    png_path = out_dir / f"{stem}_pocket_types_panel_a.png"
    org_path = out_dir / f"{stem}_pocket_types_si.org"

    fig3 = build_pocket_types_figure3js(atoms, mat_updated, results, typed)
    panel_a_image = None
    if not no_auto_panel_a:
        lite = should_use_lite_panel_a(results, typed)
        panel_fig3 = (
            build_panel_a_figure3js(atoms, mat_updated, results, typed, lite=lite)
            if lite
            else fig3
        )
        try:
            wait = panel_a_wait_ms * (3 if lite else 1)
            eye, center = _camera_eye_center(atoms)
            export_panel_a_png(panel_fig3, png_path, eye=eye, center=center, wait_ms=wait)
            panel_a_image = png_path
            mode = "lite" if lite else "full"
            print(
                f"  Panel (a) PNG {png_path.name} ({mode}) "
                f"eye=({eye[0]:.2f},{eye[1]:.2f},{eye[2]:.2f}) "
                f"target=({center[0]:.2f},{center[1]:.2f},{center[2]:.2f})"
            )
        except Exception as exc:
            print(f"  Warning: panel (a) PNG failed ({exc})", file=sys.stderr)

    if panel_a_png_only:
        print(f"Done: {title} — panel (a) PNG only")
        return {"title": title, "stem": stem, "png": str(png_path)}

    plot_four_panel(typed, tables, pdf_path, panel_a_image=panel_a_image, palette=pocket_type_palette())
    print(f"  PDF {pdf_path.name}")
    if not no_html:
        fig3.to_html(str(html_path))
        print(f"  HTML {html_path.name}")
    emit_org_si(title, meta, tables, org_path, pdf_path.name, html_path.name)
    print(f"  Org {org_path.name}")
    if write_json:
        (out_dir / f"{stem}_fingerprint_data.json").write_text(json.dumps(tables, indent=2))
    print(f"Done: {title} — {typed.n_types} types, {typed.n_pockets} pockets")
    return {"title": title, "stem": stem, "n_types": typed.n_types, "n_pockets": typed.n_pockets, "pdf": str(pdf_path)}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cif", type=Path)
    p.add_argument("--display-name", help="Report title (overrides auto from CIF stem)")
    p.add_argument("--material", choices=sorted(BATCH_MATERIALS), nargs="*")
    p.add_argument("--batch-all", action="store_true", help="All materials in BATCH_MATERIALS")
    p.add_argument("--out-dir", type=Path, default=REPO / "examples/pocket_reports")
    p.add_argument("--spacing", type=float, default=DEFAULT_SPACING_A)
    p.add_argument("--probe-radius-psd", type=float, default=DEFAULT_PROBE_RADIUS_PSD_A)
    p.add_argument("--dgrid-backend", choices=("mic", "rtree"), default="mic")
    p.add_argument("--no-auto-panel-a", action="store_true")
    p.add_argument("--panel-a-png-only", action="store_true", help="Export panel (a) PNG only (skip PDF/org/HTML)")
    p.add_argument("--panel-a-wait-ms", type=int, default=2500)
    p.add_argument("--no-html", action="store_true")
    p.add_argument("--json", action="store_true")
    args = p.parse_args(argv)

    jobs: List[tuple] = []
    if args.batch_all:
        jobs = list(BATCH_MATERIALS.items())
    elif args.material:
        for key in args.material:
            jobs.append((key, BATCH_MATERIALS[key]))
    elif args.cif:
        jobs = [
            (
                "custom",
                (
                    _expand(args.cif),
                    args.display_name or display_name_from_stem(args.cif.stem),
                    args.dgrid_backend,
                ),
            )
        ]
    else:
        p.error("Need --cif, --material, or --batch-all")

    failures = []
    for key, spec in jobs:
        cif_path, title, backend = spec[0], spec[1], spec[2]
        if not Path(cif_path).is_file():
            failures.append((key, f"CIF not found: {cif_path}"))
            continue
        try:
            generate_report(
                cif_path,
                args.out_dir,
                display_name=title,
                spacing=args.spacing,
                probe_radius_psd=args.probe_radius_psd,
                dgrid_backend=backend,
                no_html=args.no_html,
                no_auto_panel_a=args.no_auto_panel_a,
                panel_a_wait_ms=args.panel_a_wait_ms,
                panel_a_png_only=args.panel_a_png_only,
                write_json=args.json,
            )
        except Exception as exc:
            failures.append((key, str(exc)))
            print(f"FAILED {key}: {exc}", file=sys.stderr)

    if failures:
        print("\n=== FAILURES ===", file=sys.stderr)
        for key, msg in failures:
            print(f"  {key}: {msg}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

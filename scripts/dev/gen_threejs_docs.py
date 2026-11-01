"""Build IRMOF-1 Three.js scenes for docs/threejs_plotting.qmd: HTML + Playwright PNGs.

Run: PLAYWRIGHT_BROWSERS_PATH=~/Library/Caches/ms-playwright python scripts/dev/gen_threejs_docs.py
"""
import os
import sys
import warnings
from pathlib import Path

import numpy as np
from ase.io import read

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import mofography as mgr  # noqa: E402
from mofography import Figure3js  # noqa: E402
from mofography.threejs_viz import qualitative_colorscale  # noqa: E402

HTML_DIR = ROOT / "docs/assets/threejs"
IMG_DIR = ROOT / "docs/tutorial_images"
HTML_DIR.mkdir(parents=True, exist_ok=True)
IMG_DIR.mkdir(parents=True, exist_ok=True)

atoms = read(ROOT / "examples/TestCIFs/P1-IRMOF-1.cif")
atoms = mgr.update_with_covalent_bonds(atoms)
sp = 0.5
dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=sp).astype(np.float32)
seg = mgr.do_periodic_segmentation_rust(
    dgrid, atoms, spacing=sp, mask_thickness=0.05, h=0.2, min_distance=1.5,
    threshold_abs=1.5, compute_connections=True, return_dgrid=True, disable_tqdm=True)
regions, conns, maxima = seg["regions"], seg["connections"], seg["all_maxima"]
cs = qualitative_colorscale("Pastel")
hist, _, _ = mgr.compute_percolation_history_winding(seg, step=0.1)
analysis = dict(seg, percolation_history=hist, isovalue=1.6)


def fw(f=None):
    f = f or Figure3js(layout="white")
    return f.plot_framework_as_mesh_3js(atoms)


SCENES = {
    "framework": lambda: fw().plot_cell_3js(atoms, color="orange"),
    "maxima": lambda: fw().plot_cell_3js(atoms, color="orange").plot_maxima_3js(seg, atoms, cscale=cs),
    "windows": lambda: fw().plot_cell_3js(atoms, color="orange").plot_connections_3js(
        conns, maxima, regions.shape, atoms).plot_windows_3js(
        conns, name="Windows", marker_color="red", grid_shape=regions.shape, maxima=maxima),
    "surface": lambda: Figure3js(layout="white").plot_cell_3js(atoms, color="orange").plot_surface_dgrid_3js(
        seg["dgrid"], 1.6, atoms, opacity=0.6),
    "pockets": lambda: Figure3js(layout="white").plot_cell_3js(atoms, color="orange").plot_surface_dgrid_3js(
        seg["dgrid"], 1.6, atoms, region_labels=regions, cscale=cs, name="Pockets"),
    "percolation": lambda: fw().plot_cell_3js(atoms, color="orange").plot_percolated_cluster_meshes_3js(
        hist, seg["dgrid"], regions, atoms, isovalue=1.6, probe_radius=hist[-1][0]),
    "pore_graph": lambda: fw().plot_cell_3js(atoms, color="orange").plot_maxima_3js(
        seg, atoms, cscale=cs).plot_connections_3js(conns, maxima, regions.shape, atoms, cscale=cs),
}

paths = {}
for name, make in SCENES.items():
    p = HTML_DIR / f"irmof1_{name}.html"
    if "--shots-only" not in sys.argv:
        make().save_html(p)
    paths[name] = p
    print(name, f"{p.stat().st_size / 1e6:.2f} MB")

from playwright.sync_api import sync_playwright  # noqa: E402

with sync_playwright() as pw:
    b = pw.chromium.launch(args=["--use-gl=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"])
    for name, p in paths.items():
        pg = b.new_page(viewport={"width": 1280, "height": 720})
        pg.goto(p.as_uri(), wait_until="domcontentloaded", timeout=120000)
        pg.wait_for_selector("canvas", timeout=60000)
        pg.wait_for_timeout(8000)
        pg.screenshot(path=str(IMG_DIR / f"threejs_{name}.png"))
        print("shot", name)
    b.close()

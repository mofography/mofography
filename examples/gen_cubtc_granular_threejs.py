#!/usr/bin/env python3
"""Cu-BTC granular Three.js scene (Plotly-style trace composition)."""
import sys, time
from pathlib import Path
import numpy as np
from ase.io import read
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import mofography as mgr
from mofography.movies import Figure3D
from mofography import threejs_viz as tjs
CIF = ROOT / "mofography/examples/TestCIFs/P1-Cu-BTC.cif"
OUT = ROOT / "cubtc_granular_threejs.html"
SPACING = 0.5
ISOVALS = [1.0]
print("Loading Cu-BTC...")
atoms_raw = read(str(CIF))
atoms = mgr.update_with_covalent_bonds(atoms_raw)
print(f"  {len(atoms)} atoms, {len(atoms.info.get('covalent_bonds', []))} bonds")
print("Distance grid...")
t0 = time.time()
dgrid = mgr.dgrid_from_atoms_cpu(atoms, spacing=SPACING, block_size=20000).astype(np.float32)
print(f"  shape={dgrid.shape} ({time.time()-t0:.1f}s)")
print("Segmentation...")
t0 = time.time()
seg = mgr.do_periodic_segmentation_rust(dgrid, atoms, spacing=SPACING, mask_thickness=0.5, h=0.2,
    threshold_abs=1.5, min_distance=1.5, compute_connections=True, return_dgrid=False, extra_padding=5)
print(f"  done in {time.time()-t0:.1f}s")
regions, regions_pbc = seg["regions"], seg["regions_pbc"]
connections, all_maxima = seg["connections"], seg["all_maxima"]
unique_regions = seg["unique_regions"]
connections_trimmed = np.array([c for c in connections if np.any(np.isin(c[0], unique_regions))], dtype=object)
cscale = tjs.qualitative_colorscale("Pastel")
print("Building Figure3D...")
fig = Figure3D()
fig.set_camera(eye=(2.0, 0.0, 1.35), center=(0, 0, 0), up=(0, 0, 1))
tjs.add_trace(fig, tjs.plot_cell(atoms, color="orange", name="Unit cell"))
tjs.add_traces(fig, tjs.draw_framework(atoms))
for isoval in ISOVALS:
    tjs.add_trace(fig, tjs.plot_surface_dgrid(dgrid, isovalue=isoval, ase_atoms=atoms, region_labels=regions,
        name="Regions", opacity=1.0, cscale=cscale, stepsize=2, stylize=True))
    tjs.add_trace(fig, tjs.plot_surface_dgrid(dgrid, isovalue=isoval, ase_atoms=atoms, region_labels=regions_pbc,
        name="Regions PBC", opacity=0.35, cscale=cscale, stepsize=2, stylize=True))
tjs.add_traces(fig, tjs.plot_maxima(seg, atoms, name="Maxima (distinct)", cscale=cscale))
tjs.add_traces(fig, tjs.plot_connections_as_it_is(connections_trimmed, all_maxima, dgrid.shape, atoms,
    cscale=cscale, plot_periodic=True, plot_disconnected=False))
tjs.add_traces(fig, tjs.plot_windows(connections_trimmed, name="Windows", marker_color="yellow",
    marker_size=10, marker_opacity=0.8, grid_shape=dgrid.shape, maxima=all_maxima))
fig.meta["description"] = "Cu-BTC granular Three.js"
fig.to_html(str(OUT))
print(f"Saved {OUT} ({OUT.stat().st_size/1024/1024:.1f} MB)")

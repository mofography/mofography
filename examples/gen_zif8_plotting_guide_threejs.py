#!/usr/bin/env python3
"""Generate ZIF-8 Three.js HTML for manual pore-graph / window-filter testing."""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import mofography as mgr
from mofography.workflow import assign_pocket_types, show_mof_threejs

CIF = ROOT / "examples/TestCIFs/ZIF-8-RASPA.cif"
OUT = ROOT / "examples/viz_outputs/ZIF-8-RASPA_spacing0.25_threejs.html"
OUT.parent.mkdir(parents=True, exist_ok=True)

SPACING = 0.25
MASK_THICKNESS = 0.05
H_VALUE = 0.2
MIN_DISTANCE = 1.5
THRESHOLD_ABS = 1.5
EXTRA_PADDING = 5.0

print(f"Loading {CIF.name}...")
from ase.io import read

mat_atoms = read(str(CIF))
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

print(f"Distance grid (spacing={SPACING} Å)...")
t0 = time.time()
dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=SPACING).astype(np.float32)
print(f"  shape={dgrid.shape} ({time.time() - t0:.1f}s)")

print("Segmentation...")
t0 = time.time()
results = mgr.do_periodic_segmentation_rust(
    dgrid,
    mat_atoms,
    spacing=SPACING,
    mask_thickness=MASK_THICKNESS,
    h=H_VALUE,
    min_distance=MIN_DISTANCE,
    threshold_abs=THRESHOLD_ABS,
    extra_padding=EXTRA_PADDING,
    return_dgrid=True,
)
print(f"  done in {time.time() - t0:.1f}s")

typed = assign_pocket_types(results, mat_atoms)
print("Building Three.js scene...")
path = show_mof_threejs(
    results,
    mat_atoms,
    typed_result=typed,
    framework_atoms=mat_updated,
    output_html=str(OUT),
    plot_framework_representation="instanced",
    plot_layout="white",
    plot_spacing=SPACING,
    plot_windows_style={
        "color": 0xFF0000,
        "size": 8.0,
        "symbol": "circle",
    },
)
print(f"Saved {path} ({OUT.stat().st_size / 1024 / 1024:.1f} MB)")

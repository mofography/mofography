#!/usr/bin/env python3
"""Parity check: Python dgrid MIC helpers vs Rust dgrid_from_atoms."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from ase.io import read

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

import mofography as mgr
from mofography.compute import dgrid_from_atoms_cpu_no_aabb
from mofography.rust_fun import dgrid_from_atoms_rust
from mofography.vdw_radii import default_vdw_radii_angstrom

TOL = 1e-5
SPACING = 0.5

CASES = [
    REPO / "examples/TestCIFs/FOTNIN_clean.cif",
    REPO / "examples/TestCIFs/WUPNIF_clean.cif",
    REPO / "examples/TestCIFs/TOVKOG_clean.cif",
    REPO / "examples/TestCIFs/P1-Cu-BTC.cif",
    REPO / "examples/TestCIFs/NU-1000_S1.cif",
]


def main() -> int:
    failed = False
    for cif in CASES:
        if not cif.exists():
            print(f"SKIP missing: {cif.name}")
            continue
        atoms = read(cif)
        radii = default_vdw_radii_angstrom(atoms)
        rust = dgrid_from_atoms_rust(atoms, spacing=SPACING, radii=radii).astype(np.float64)
        py = dgrid_from_atoms_cpu_no_aabb(
            atoms, spacing=SPACING, radii=radii, block_size=5000
        ).astype(np.float64)
        if rust.shape != py.shape:
            print(f"FAIL {cif.name}: shape {rust.shape} vs {py.shape}")
            failed = True
            continue
        max_diff = float(np.max(np.abs(rust - py)))
        mean_diff = float(np.mean(np.abs(rust - py)))
        status = "PASS" if max_diff < TOL else "FAIL"
        print(f"{status} {cif.name}: max={max_diff:.3e} mean={mean_diff:.3e}")
        if max_diff >= TOL:
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())

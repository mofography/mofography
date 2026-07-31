#!/usr/bin/env python3
"""Minimal ZIF-20 Three.js HTML: framework + unit cell + aromatic rings."""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import mofography as mgr
from mofography.movies import (
    MovieFigure,
    make_aromatic_rings_trace,
    make_cell_trace,
    make_framework_instanced_trace,
)

CIF = Path.home() / "repos/RASPA2/structures/mofs/cif/ZIF-20.cif"
OUT = ROOT / "examples/viz_outputs/ZIF-20_rings_check_threejs.html"
OUT.parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    from ase.io import read

    if not CIF.is_file():
        raise FileNotFoundError(f"ZIF-20 CIF not found: {CIF}")

    print(f"Loading {CIF.name}...")
    atoms = read(str(CIF))
    print(f"  {len(atoms):,} atoms")

    t0 = time.time()
    atoms = mgr.update_with_covalent_bonds(atoms)
    print(f"  bonds: {len(atoms.info['covalent_bonds']):,} ({time.time() - t0:.1f}s)")

    t0 = time.time()
    cycles = mgr.compute_aromatic_cycles(atoms, force_recompute=True)
    sizes = sorted({len(c) for c in cycles})
    print(f"  aromatic rings: {len(cycles):,}  sizes={sizes} ({time.time() - t0:.1f}s)")

    print("Building Three.js scene (framework + cell + rings)...")
    fw = make_framework_instanced_trace(
        atoms,
        representation="ball_and_stick",
        trace_id="framework",
    )
    rings = make_aromatic_rings_trace(atoms, cycles=cycles)

    fig = MovieFigure()
    fig.add_trace(fw["atoms"])
    fig.add_trace(fw["bonds"])
    fig.add_trace(make_cell_trace(atoms))
    fig.add_trace(rings["circles"])
    fig.add_trace(rings["highlights"])
    if "fills_hex" in rings:
        fig.add_trace(rings["fills_hex"])
    if "fills_pent" in rings:
        fig.add_trace(rings["fills_pent"])

    fig.to_html(OUT)
    print(f"Saved {OUT} ({OUT.stat().st_size / 1024 / 1024:.1f} MB)")


if __name__ == "__main__":
    main()

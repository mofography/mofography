#!/usr/bin/env python3
"""NU-1000 granular Three.js test with cell-matrix fix (C = M.T) for Rust dgrid."""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import mofography as mgr
from mofography import threejs_viz as tjs
from mofography.cif_io import read_cif_to_ase

CIF = ROOT / "mofography/examples/TestCIFs/NU-1000_S1.cif"
OUT = Path("/tmp/nu1000_granular_threejs_cellfix.html")
SPACING = 0.5


def dgrid_from_atoms_rust_fixed(ase_atoms, spacing, radii=None):
    import mgr_rust
    from ase.data import vdw_radii as _ase_vdw
    from ase.geometry import cellpar_to_cell

    if radii is None:
        radii = _ase_vdw[ase_atoms.get_atomic_numbers()]
    fpoints = ase_atoms.get_scaled_positions() - 0.5
    cell = cellpar_to_cell(ase_atoms.cell.cellpar()).T
    nx, ny, nz = (ase_atoms.cell.cellpar()[:3] / spacing).astype(int) + 1
    return mgr_rust.dgrid_from_atoms(fpoints, cell, radii, int(nx), int(ny), int(nz))


def compare_dgrid_wrappers(atoms, spacing: float) -> None:
    from mofography.compute import dgrid_from_atoms_cpu

    print("\n=== Dgrid wrapper comparison (NU-1000) ===")
    t0 = time.perf_counter()
    d_legacy = mgr.dgrid_from_atoms_rust(atoms, spacing=spacing).astype(np.float32)
    t_legacy = time.perf_counter() - t0

    t0 = time.perf_counter()
    d_fixed = dgrid_from_atoms_rust_fixed(atoms, spacing=spacing).astype(np.float32)
    t_fixed = time.perf_counter() - t0

    t0 = time.perf_counter()
    d_cpu = dgrid_from_atoms_cpu(atoms, spacing=spacing, block_size=20000).astype(np.float32)
    t_cpu = time.perf_counter() - t0

    def stats(a, b, label):
        d = np.abs(a - b)
        print(
            f"  {label}: mean={d.mean():.4f} max={d.max():.4f} A  "
            f"pct>0.01={100 * np.mean(d > 0.01):.1f}%"
        )

    stats(d_legacy, d_cpu, "legacy rust vs cpu_aabb")
    stats(d_fixed, d_cpu, "fixed rust   vs cpu_aabb")
    stats(d_legacy, d_fixed, "legacy vs fixed")
    print(f"  timings: legacy={t_legacy:.1f}s fixed={t_fixed:.1f}s cpu={t_cpu:.1f}s\n")


def main() -> None:
    print(f"Loading {CIF.name}...")
    atoms = read_cif_to_ase(CIF)
    print(f"  {len(atoms):,} atoms  gamma={atoms.cell.cellpar()[5]:.1f} deg")

    compare_dgrid_wrappers(atoms, SPACING)

    print("Covalent bonds...")
    t0 = time.perf_counter()
    atoms = mgr.update_with_covalent_bonds(atoms)
    print(f"  {time.perf_counter() - t0:.1f}s - {len(atoms.info['covalent_bonds']):,} bonds")

    print("Aromatic rings...")
    t0 = time.perf_counter()
    mgr.compute_aromatic_cycles(atoms)
    print(f"  {time.perf_counter() - t0:.1f}s - {len(atoms.info.get('aromatic_cycles', [])):,} cycles")

    print("Distance grid (C=M.T via rust_fun)...")
    t0 = time.perf_counter()
    dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=SPACING).astype(np.float32)
    print(f"  {time.perf_counter() - t0:.1f}s - shape={dgrid.shape}")

    print("Segmentation...")
    t0 = time.perf_counter()
    seg = mgr.do_periodic_segmentation_rust(
        dgrid,
        atoms,
        spacing=SPACING,
        mask_thickness=0.5,
        h=0.2,
        threshold_abs=1.5,
        min_distance=1.5,
        extra_padding=5,
        compute_connections=True,
        return_dgrid=False,
    )
    print(f"  {time.perf_counter() - t0:.1f}s")

    regions = seg["regions"]
    regions_pbc = seg["regions_pbc"]
    connections = seg["connections"]
    all_maxima = seg["all_maxima"]
    unique_regions = seg["unique_regions"]
    connections_trimmed = [c for c in connections if np.any(np.isin(c[0], unique_regions))]

    print("Building Three.js figure...")
    t0 = time.perf_counter()
    cscale = tjs.qualitative_colorscale("Pastel")

    fig = tjs.Figure3JS(layout="dark", stylize=True)
    fig.set_camera(eye=(2.0, 0.0, 1.35), center=(0, 0, 0), up=(0, 0, 1))

    fig.draw_framework(atoms, representation="scatter")
    fig.plot_cell(atoms, color="orange")

    fig.plot_surface_dgrid(
        dgrid,
        isovalue=1.0,
        ase_atoms=atoms,
        region_labels=regions,
        name="Regions",
        opacity=1.0,
        cscale=cscale,
        stepsize=1,
        stylize=True,
    )
    fig.plot_surface_dgrid(
        dgrid,
        isovalue=1.0,
        ase_atoms=atoms,
        region_labels=regions_pbc,
        name="Regions PBC",
        opacity=0.35,
        cscale=cscale,
        stepsize=1,
        stylize=True,
    )

    fig.plot_maxima(seg, atoms, name="Maxima", cscale=cscale)
    fig.plot_connections_as_it_is(
        connections_trimmed,
        all_maxima,
        dgrid.shape,
        atoms,
        cscale=cscale,
        plot_periodic=True,
        plot_disconnected=False,
    )
    fig.plot_windows(
        connections_trimmed,
        name="Windows",
        marker_color="yellow",
        marker_size=10,
        marker_opacity=0.8,
        grid_shape=dgrid.shape,
        maxima=all_maxima,
    )

    fig.to_html(str(OUT))
    mb = OUT.stat().st_size / 1024 / 1024
    print(f"  {time.perf_counter() - t0:.1f}s - {mb:.0f} MB")
    print(f"Saved to {OUT}")


if __name__ == "__main__":
    main()

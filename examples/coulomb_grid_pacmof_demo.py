#!/usr/bin/env python3
"""Tier 1 demo: PACMOF → LJ + Coulomb grids → guest_energy_kelvin.

Examples::

    python examples/coulomb_grid_pacmof_demo.py --cif cubtc
    python examples/coulomb_grid_pacmof_demo.py --cif irmof --spacing 0.5
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

CIFS = {
    "cubtc": ROOT / "mofography/examples/TestCIFs/P1-Cu-BTC.cif",
    "irmof": ROOT / "mofography/examples/TestCIFs/P1-IRMOF-1.cif",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cif", choices=sorted(CIFS), default="cubtc")
    parser.add_argument("--spacing", type=float, default=0.5)
    parser.add_argument("--cutoff", type=float, default=12.5)
    parser.add_argument("--probe", default="He", help="LJ probe symbol")
    parser.add_argument("--q-guest", type=float, default=-0.35, help="Guest charge (e)")
    parser.add_argument("--validate", action="store_true", help="Run interpolation/FD checks")
    args = parser.parse_args()

    from mofography.charges import get_pacmof_charges_single, resolve_pacmof_model_path
    from mofography.compute import build_guest_interpolators, guest_energy_kelvin
    from mofography.coulomb_validate import (
        validate_derivative_channels_fd,
        validate_grid_channel0_vs_direct,
        validate_interpolation_vs_direct,
    )

    cif = CIFS[args.cif]
    print(f"PACMOF model: {resolve_pacmof_model_path()}")
    print(f"CIF: {cif.name}")

    atoms = get_pacmof_charges_single(cif)
    ch = np.asarray(atoms.info["_atom_site_charge"])
    print(f"  atoms={len(atoms)}  |sum q|={abs(ch.sum()):.2e}  q∈[{ch.min():.3f}, {ch.max():.3f}]")

    print(f"Building grids (spacing={args.spacing} Å, cutoff={args.cutoff} Å, probe={args.probe})...")
    lj_ip, coul_ip = build_guest_interpolators(
        atoms, spacing=args.spacing, cutoff=args.cutoff, probe_symbol=args.probe
    )
    if coul_ip is None:
        raise RuntimeError("Coulomb interpolator missing — charges not assigned?")

    r_ase = np.array([[0.65, 0.65, 0.65]])  # ASE frac [0,1); away from singularities
    e_lj = lj_ip(r_ase)
    phi_k = coul_ip(r_ase)
    e_tot = guest_energy_kelvin(r_ase, args.q_guest, lj_ip, coul_ip)
    print(f"\nAt ASE frac {r_ase[0]} (guest q={args.q_guest} e):")
    print(f"  E_LJ ({args.probe})  = {float(np.asarray(e_lj).ravel()[0]):.2f} K")
    print(f"  φ_K                  = {float(np.asarray(phi_k).ravel()[0]):.2f} K/e")
    print(f"  E_coul               = {args.q_guest * float(np.asarray(phi_k).ravel()[0]):.2f} K")
    print(f"  E_total              = {float(np.asarray(e_tot).ravel()[0]):.2f} K")

    if args.validate:
        from mofography.coulomb_grid import coulomb_tricubic_from_atoms_rust

        coulomb8 = coulomb_tricubic_from_atoms_rust(atoms, spacing=args.spacing, cutoff=args.cutoff)
        g_max, g_rms, _ = validate_grid_channel0_vs_direct(atoms, coulomb8, cutoff=args.cutoff)
        i_max, i_rms, _ = validate_interpolation_vs_direct(atoms, coul_ip, cutoff=args.cutoff)
        d_max, d_rms = validate_derivative_channels_fd(coulomb8)
        print("\nValidation:")
        print(f"  grid ch0 vs direct:  max_rel={g_max:.2e}  rms={g_rms:.2e}")
        print(f"  interp vs direct:    max_rel={i_max:.2e}  rms={i_rms:.2e}")
        print(f"  deriv ch1-3 vs FD:   max_rel={d_max:.2e}  rms={d_rms:.2e}")


if __name__ == "__main__":
    main()

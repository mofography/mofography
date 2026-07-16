#!/usr/bin/env python3
"""Tier-1 validation: PACMOF charges → Coulomb grid → interpolation checks."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

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
    parser.add_argument("--spacing", type=float, default=1.0, help="Å (use 0.5 for production)")
    parser.add_argument("--cutoff", type=float, default=12.5)
    parser.add_argument("--n-interp", type=int, default=12, help="Random interpolation samples")
    args = parser.parse_args()

    from mofography.charges import get_pacmof_charges_single, resolve_pacmof_model_path
    from mofography.coulomb_grid import coulomb_tricubic_from_atoms_rust
    from mofography.coulomb_validate import (
        validate_derivative_channels_fd,
        validate_grid_channel0_vs_direct,
        validate_interpolation_vs_direct,
    )
    from mofography.compute import get_energy_interpolator_tricubic

    cif = CIFS[args.cif]
    print(f"PACMOF model: {resolve_pacmof_model_path()}")
    print(f"CIF: {cif.name}  spacing={args.spacing} Å  cutoff={args.cutoff} Å")

    atoms = get_pacmof_charges_single(cif)
    import numpy as np

    ch = np.asarray(atoms.info["_atom_site_charge"])
    print(f"  atoms={len(atoms)}  |sum q|={abs(ch.sum()):.2e}")

    coulomb8 = coulomb_tricubic_from_atoms_rust(atoms, spacing=args.spacing, cutoff=args.cutoff)
    print(f"  Coulomb grid shape {coulomb8.shape}")

    coul_ip = get_energy_interpolator_tricubic(coulomb8[0], coulomb8)

    g_max, g_rms, g_n = validate_grid_channel0_vs_direct(
        atoms, coulomb8, cutoff=args.cutoff, margin=2
    )
    print(f"Grid ch0 vs direct Wolf ({g_n} nodes): max_rel={g_max:.2e}  rms_rel={g_rms:.2e}")

    i_max, i_rms, i_n = validate_interpolation_vs_direct(
        atoms, coul_ip, cutoff=args.cutoff, n_samples=args.n_interp
    )
    print(f"Tricubic vs direct Wolf ({i_n} samples): max_rel={i_max:.2e}  rms_rel={i_rms:.2e}")

    d_max, d_rms = validate_derivative_channels_fd(coulomb8, margin=2)
    print(f"Derivative ch1-3 vs FD on ch0: max_rel={d_max:.2e}  rms_rel={d_rms:.2e}")

    ok = g_max < 1e-8 and i_max < 0.12
    print("PASS (grid + interpolation)" if ok else "FAIL")
    if d_max > 0.2:
        print("  note: derivative FD noisy at coarse spacing — retry with --spacing 0.5")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

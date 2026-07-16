"""Coulomb electrostatic potential grid builder — Wolf DSF (φ_K in Kelvin)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)


def coulomb_tricubic_from_cif_rust(
    path_to_cif: Union[str, Path],
    spacing: float,
    cutoff: float = 12.5,
    alpha: Optional[float] = None,
    return_ncells: bool = False,
    assign_charges: str = "pacmof",
    **pacmof_kwargs,
):
    """Build Coulomb grid from a CIF, assigning charges first if needed.

    Parameters
    ----------
    assign_charges
        ``'pacmof'`` — run PACMOF DDEC charges via
        :func:`mofography.charges.get_pacmof_charges_single`.
        ``'none'`` — CIF must already carry ``_atom_site_charge``.
    **pacmof_kwargs
        Passed to :func:`get_pacmof_charges_single` when ``assign_charges='pacmof'``.
    """
    path_to_cif = Path(path_to_cif)
    if assign_charges == "pacmof":
        from mofography.charges import get_pacmof_charges_single

        atoms = get_pacmof_charges_single(path_to_cif, **pacmof_kwargs)
    elif assign_charges == "none":
        from ase.io import read

        atoms = read(str(path_to_cif))
    else:
        raise ValueError(f"Unknown assign_charges={assign_charges!r}; use 'pacmof' or 'none'")

    return coulomb_tricubic_from_atoms_rust(
        atoms,
        spacing=spacing,
        cutoff=cutoff,
        alpha=alpha,
        return_ncells=return_ncells,
    )


def _coulomb_unit_cell_grid_inputs(ase_atoms, mat_atoms_repeat, charges_repeat, spacing, cutoff):
    """Shared unit-cell grid geometry for Coulomb grids (mirrors LJ path)."""
    from ase.geometry import cellpar_to_cell
    from mofography.compute import get_repetitions_to_meet_cutoff
    from mofography.rust_fun import _unit_cell_grid_geometry

    frac, cell, nx_unit, ny_unit, nz_unit, ex, ey, ez = _unit_cell_grid_geometry(
        ase_atoms, mat_atoms_repeat, spacing, cutoff
    )
    nx_c, ny_c, nz_c = get_repetitions_to_meet_cutoff(ase_atoms, cutoff=cutoff)
    n_repeat = nx_c * ny_c * nz_c
    if len(charges_repeat) != len(mat_atoms_repeat):
        raise ValueError(
            f"charges_repeat length {len(charges_repeat)} != "
            f"replicated atoms {len(mat_atoms_repeat)}"
        )
    charges = np.asarray(charges_repeat, dtype=np.float64)
    return frac, cell, nx_unit, ny_unit, nz_unit, charges, ex, ey, ez


def coulomb_tricubic_from_atoms_rust(
    ase_atoms,
    spacing: float,
    cutoff: float = 12.5,
    alpha: Optional[float] = None,
    return_ncells: bool = False,
) -> Union[np.ndarray, Tuple[Optional[np.ndarray], Tuple[int, int, int]], None]:
    """Build an 8-channel Coulomb potential grid φ_K using Wolf DSF.

    Channel 0 stores φ_K (Kelvin per unit probe charge +1 e). At query time:
    ``E_coul_K = q_guest * interpolate(grid)``.

    Returns ``None`` (or ``(None, ncells)`` if ``return_ncells``) when the
    structure has no ``_atom_site_charge``. Coulomb is opt-in — charged CIFs
    do not auto-build this grid.
    """
    import mgr_rust
    from mofography.compute import get_repetitions_to_meet_cutoff
    from mofography.constants import COULOMB_CONVERSION_K
    from mofography.charges import CHARGE_INFO_KEY, normalize_charge_info

    try:
        normalize_charge_info(ase_atoms)
    except KeyError:
        logger.info("No %s — Coulomb grid skipped.", CHARGE_INFO_KEY)
        if return_ncells:
            return None, get_repetitions_to_meet_cutoff(ase_atoms, cutoff=cutoff)
        return None

    charges_unit = np.asarray(ase_atoms.info[CHARGE_INFO_KEY], dtype=np.float64)
    total_charge = float(charges_unit.sum())
    if abs(total_charge) > 1e-6:
        logger.warning(
            "Framework charge non-neutral: Sum q = %.2e. Wolf method may be inaccurate.",
            total_charge,
        )

    nx_c, ny_c, nz_c = get_repetitions_to_meet_cutoff(ase_atoms, cutoff=cutoff)
    n_repeat = nx_c * ny_c * nz_c
    mat_atoms_repeat = ase_atoms.repeat((nx_c, ny_c, nz_c))
    charges_repeat = np.tile(charges_unit, n_repeat)

    frac, cell, nx_unit, ny_unit, nz_unit, charges, ex, ey, ez = _coulomb_unit_cell_grid_inputs(
        ase_atoms, mat_atoms_repeat, charges_repeat, spacing, cutoff
    )

    if alpha is None:
        alpha = 4.0 / cutoff

    logger.info(
        "Building Coulomb grid: %dx%dx%d, alpha=%.4f, R_c=%.2f, alpha*R_c=%.2f, "
        "COULOMB_CONVERSION_K=%.6g",
        nx_unit,
        ny_unit,
        nz_unit,
        alpha,
        cutoff,
        alpha * cutoff,
        COULOMB_CONVERSION_K,
    )

    egrid8_flat = mgr_rust.coulomb_tricubic_from_atoms(
        np.asarray(frac, dtype=np.float64),
        np.asarray(cell, dtype=np.float64),
        charges,
        int(nx_unit),
        int(ny_unit),
        int(nz_unit),
        float(cutoff),
        alpha=float(alpha),
        extent_x=float(ex),
        extent_y=float(ey),
        extent_z=float(ez),
    )
    egrid8 = np.asarray(egrid8_flat, dtype=np.float64).reshape(8, nx_unit, ny_unit, nz_unit)

    if return_ncells:
        return egrid8, (nx_c, ny_c, nz_c)
    return egrid8

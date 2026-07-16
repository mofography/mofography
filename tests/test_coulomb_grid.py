"""Tests for Wolf DSF Coulomb tricubic grid builder."""

import numpy as np
import pytest
from ase import Atoms
from ase.geometry import cellpar_to_cell
from scipy.special import erfc

from mofography.constants import COULOMB_CONVERSION_K
from mofography.coulomb_grid import coulomb_tricubic_from_atoms_rust
from mofography.compute import guest_energy_kelvin, get_energy_interpolator_tricubic


def _charged_diatomic():
    cell = 20.0 * np.eye(3)
    positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
    charges = [0.5, -0.5]
    atoms = Atoms(symbols="CC", positions=positions, cell=cell, pbc=True)
    atoms.info["_atom_site_charge"] = charges
    return atoms


def test_no_charges_returns_none():
    atoms = Atoms("He", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
    assert coulomb_tricubic_from_atoms_rust(atoms, spacing=1.0) is None


def test_coulomb_grid_shape_and_finite():
    atoms = _charged_diatomic()
    grid8 = coulomb_tricubic_from_atoms_rust(atoms, spacing=2.0, cutoff=12.5)
    assert grid8 is not None
    assert grid8.shape[0] == 8
    # Interior channel 0 should be finite (corners may blow up near 1/r singularity)
    assert np.all(np.isfinite(grid8[0, 1:-1, 1:-1, 1:-1]))


def test_direct_sum_matches_grid_center():
    """Channel 0 at a grid node vs brute-force Wolf sum."""
    atoms = _charged_diatomic()
    spacing = 2.0
    cutoff = 12.5
    alpha = 4.0 / cutoff
    grid8 = coulomb_tricubic_from_atoms_rust(atoms, spacing=spacing, cutoff=cutoff, alpha=alpha)
    nx, ny, nz = grid8.shape[1:]

    frac = atoms.get_scaled_positions() - 0.5
    cell = cellpar_to_cell(atoms.cell.cellpar()).T
    charges = np.asarray(atoms.info["_atom_site_charge"], dtype=np.float64)

    i, j, k = nx // 2, ny // 2, nz // 2
    sx = -0.5 + i / (nx - 1)
    sy = -0.5 + j / (ny - 1)
    sz = -0.5 + k / (nz - 1)

    erfc_aR = erfc(alpha * cutoff)
    S = erfc_aR / cutoff
    M = S / cutoff + 2 * alpha * np.exp(-alpha**2 * cutoff**2) / (np.sqrt(np.pi) * cutoff)

    e_direct = 0.0
    for af, q in zip(frac, charges):
        dfx = sx - af[0]
        dfy = sy - af[1]
        dfz = sz - af[2]
        dfx -= np.round(dfx)
        dfy -= np.round(dfy)
        dfz -= np.round(dfz)
        d = cell @ np.array([dfx, dfy, dfz])
        r = np.linalg.norm(d)
        if r >= cutoff or r < 1e-12:
            continue
        e_direct += COULOMB_CONVERSION_K * q * (erfc(alpha * r) / r - S + M * (r - cutoff))

    grid_val = grid8[0, i, j, k]
    assert abs(grid_val - e_direct) / max(abs(e_direct), 1e-6) < 1e-5


def test_guest_energy_kelvin_combines():
    atoms = _charged_diatomic()
    coulomb8 = coulomb_tricubic_from_atoms_rust(atoms, spacing=2.0, cutoff=12.5)
    ch0 = coulomb8[0]
    coulomb_ip = get_energy_interpolator_tricubic(ch0, coulomb8)

    def lj_ip(frac):
        xi = np.atleast_2d(np.asarray(frac, dtype=np.float64))
        return np.zeros(xi.shape[0])

    r = np.array([[0.6, 0.6, 0.6]])  # ASE fractional coords
    e = guest_energy_kelvin(r, q_guest=0.3, lj_interpolator=lj_ip, coulomb_interpolator=coulomb_ip)
    phi = coulomb_ip(r)
    assert abs(float(np.asarray(e).ravel()[0]) - 0.3 * float(np.asarray(phi).ravel()[0])) < 1e-6


def test_constants_value():
    assert abs(COULOMB_CONVERSION_K - 138_935.4834964017) < 1e-6

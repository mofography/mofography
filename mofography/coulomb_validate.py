"""Validation helpers for Wolf DSF Coulomb tricubic grids."""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
from ase.geometry import cellpar_to_cell
from scipy.special import erfc

from mofography.constants import COULOMB_CONVERSION_K
from mofography.compute import get_repetitions_to_meet_cutoff


def grid_frac_to_ase_frac(grid_frac: Sequence[float]) -> np.ndarray:
    """Convert grid-build fractional coords (−0.5 origin) to ASE [0, 1) for tricubic query."""
    g = np.asarray(grid_frac, dtype=np.float64)
    return g + 0.5


def ase_frac_to_grid_frac(ase_frac: Sequence[float]) -> np.ndarray:
    """Inverse of :func:`grid_frac_to_ase_frac` (extent = 1 unit cell)."""
    a = np.asarray(ase_frac, dtype=np.float64)
    return a - 0.5


def wolf_phi_k_direct(
    ase_atoms,
    frac_point: Sequence[float],
    cutoff: float,
    alpha: Optional[float] = None,
) -> float:
    """Brute-force Wolf DSF φ_K at one fractional probe point (probe q = +1 e)."""
    if alpha is None:
        alpha = 4.0 / cutoff

    nx_c, ny_c, nz_c = get_repetitions_to_meet_cutoff(ase_atoms, cutoff=cutoff)
    if nx_c == ny_c == nz_c == 1:
        frac_atoms = ase_atoms.get_scaled_positions() - 0.5
        cell = cellpar_to_cell(ase_atoms.cell.cellpar()).T
    else:
        rep = ase_atoms.repeat((nx_c, ny_c, nz_c))
        frac_atoms = rep.get_scaled_positions() - 0.5
        cell = cellpar_to_cell(rep.cell.cellpar()).T

    charges = np.asarray(ase_atoms.info["_atom_site_charge"], dtype=np.float64)
    charges = np.tile(charges, nx_c * ny_c * nz_c)

    sx, sy, sz = float(frac_point[0]), float(frac_point[1]), float(frac_point[2])
    erfc_aR = erfc(alpha * cutoff)
    shift = erfc_aR / cutoff
    slope = shift / cutoff + 2 * alpha * np.exp(-(alpha * cutoff) ** 2) / (np.sqrt(np.pi) * cutoff)

    energy = 0.0
    for af, q in zip(frac_atoms, charges):
        dfx, dfy, dfz = sx - af[0], sy - af[1], sz - af[2]
        dfx -= np.round(dfx)
        dfy -= np.round(dfy)
        dfz -= np.round(dfz)
        d = cell @ np.array([dfx, dfy, dfz])
        r = float(np.linalg.norm(d))
        if r >= cutoff or r < 1e-12:
            continue
        energy += COULOMB_CONVERSION_K * q * (
            erfc(alpha * r) / r - shift + slope * (r - cutoff)
        )
    return float(energy)


def validate_interpolation_vs_direct(
    ase_atoms,
    coulomb_interpolator,
    cutoff: float = 12.5,
    alpha: Optional[float] = None,
    n_samples: int = 12,
    rng: Optional[np.random.Generator] = None,
    frac_min: float = 0.08,
    frac_max: float = 0.92,
) -> Tuple[float, float, int]:
    """Max and RMS relative error: tricubic interpolate vs direct Wolf sum.

    Returns (max_rel_err, rms_rel_err, n_ok).
    """
    rng = rng or np.random.default_rng(0)
    errs = []
    for _ in range(n_samples):
        ase_frac = rng.uniform(frac_min, frac_max, size=3)
        grid_frac = ase_frac_to_grid_frac(ase_frac)
        direct = wolf_phi_k_direct(ase_atoms, grid_frac, cutoff=cutoff, alpha=alpha)
        interp = float(np.asarray(coulomb_interpolator(np.array([ase_frac]))).ravel()[0])
        if not np.isfinite(direct) or not np.isfinite(interp):
            continue
        denom = max(abs(direct), 1.0)
        errs.append(abs(interp - direct) / denom)
    if not errs:
        return float("nan"), float("nan"), 0
    errs = np.asarray(errs)
    return float(errs.max()), float(np.sqrt(np.mean(errs**2))), len(errs)


def validate_grid_channel0_vs_direct(
    ase_atoms,
    coulomb8: np.ndarray,
    extent_x: float = 1.0,
    extent_y: float = 1.0,
    extent_z: float = 1.0,
    cutoff: float = 12.5,
    alpha: Optional[float] = None,
    margin: int = 2,
) -> Tuple[float, float, int]:
    """Compare grid channel 0 at interior nodes to direct Wolf sum."""
    _, nx, ny, nz = coulomb8.shape
    errs = []
    for i in range(margin, nx - margin):
        for j in range(margin, ny - margin):
            for k in range(margin, nz - margin):
                sx = -0.5 + (i / (nx - 1)) * extent_x
                sy = -0.5 + (j / (ny - 1)) * extent_y
                sz = -0.5 + (k / (nz - 1)) * extent_z
                direct = wolf_phi_k_direct(ase_atoms, (sx, sy, sz), cutoff=cutoff, alpha=alpha)
                grid = float(coulomb8[0, i, j, k])
                if not np.isfinite(direct) or not np.isfinite(grid):
                    continue
                denom = max(abs(direct), 1.0)
                errs.append(abs(grid - direct) / denom)
    if not errs:
        return float("nan"), float("nan"), 0
    errs = np.asarray(errs)
    return float(errs.max()), float(np.sqrt(np.mean(errs**2))), len(errs)


def validate_derivative_channels_fd(
    coulomb8: np.ndarray,
    extent_x: float = 1.0,
    extent_y: float = 1.0,
    extent_z: float = 1.0,
    margin: int = 2,
    min_magnitude: float = 50.0,
) -> Tuple[float, float]:
    """FD check: ch1–3 vs central differences on ch0 at interior nodes.

    Channels store dφ/d(grid_index). Skips nodes where |ch1| and |FD| are both
    below ``min_magnitude`` (near Coulomb singularities). Returns (max_rel, rms_rel).
    """
    _, nx, ny, nz = coulomb8.shape
    errs = []
    ch0 = coulomb8[0]
    for i in range(margin, nx - margin):
        for j in range(margin, ny - margin):
            for k in range(margin, nz - margin):
                checks = [
                    ((ch0[i + 1, j, k] - ch0[i - 1, j, k]) / 2.0, coulomb8[1, i, j, k]),
                    ((ch0[i, j + 1, k] - ch0[i, j - 1, k]) / 2.0, coulomb8[2, i, j, k]),
                    ((ch0[i, j, k + 1] - ch0[i, j, k - 1]) / 2.0, coulomb8[3, i, j, k]),
                ]
                for fd, stored in checks:
                    if not np.isfinite(fd) or not np.isfinite(stored):
                        continue
                    if abs(stored) < min_magnitude and abs(fd) < min_magnitude:
                        continue
                    denom = max(abs(stored), abs(fd), 1.0)
                    errs.append(abs(fd - stored) / denom)
    if not errs:
        return float("nan"), float("nan")
    errs = np.asarray(errs)
    return float(errs.max()), float(np.sqrt(np.mean(errs**2)))

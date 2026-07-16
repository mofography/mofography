"""PACMOF + Coulomb grid integration (requires pacmof + Model_RF_DDEC.pkl)."""

from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CUBTC = ROOT / "mofography/examples/TestCIFs/P1-Cu-BTC.cif"


def _pacmof_model_available() -> bool:
    try:
        from mofography.charges import resolve_pacmof_model_path

        resolve_pacmof_model_path()
        return True
    except (ImportError, FileNotFoundError):
        return False


pytestmark = pytest.mark.skipif(
    not _pacmof_model_available(),
    reason="pacmof or Model_RF_DDEC.pkl not available",
)


@pytest.mark.slow
def test_pacmof_charges_neutral_cubtc():
    from mofography.charges import get_pacmof_charges_single

    atoms = get_pacmof_charges_single(CUBTC)
    ch = np.asarray(atoms.info["_atom_site_charge"], dtype=np.float64)
    assert len(ch) == len(atoms)
    assert abs(ch.sum()) < 1e-3


@pytest.mark.slow
def test_pacmof_coulomb_grid_cubtc():
    from mofography.coulomb_grid import coulomb_tricubic_from_cif_rust

    grid8 = coulomb_tricubic_from_cif_rust(CUBTC, spacing=1.5, cutoff=12.5)
    assert grid8 is not None
    assert grid8.shape[0] == 8
    assert np.all(np.isfinite(grid8[0, 1:-1, 1:-1, 1:-1]))

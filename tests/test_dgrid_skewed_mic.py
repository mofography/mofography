"""dgrid_from_atoms must match an exact minimum image on skewed cells."""
import itertools

import numpy as np
import pytest

mgr_rust = pytest.importorskip("mgr_rust")
cellpar_to_cell = pytest.importorskip("ase.geometry").cellpar_to_cell

SHIFTS = np.array(list(itertools.product(range(-2, 3), repeat=3)), dtype=float)


def _reference(frac_atoms, cell_rows, radii, shape):
    axes = [-0.5 + np.arange(n) / (n - 1) for n in shape]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    best = np.full(len(grid), np.inf)
    for atom, radius in zip(frac_atoms, radii):
        diff = grid - atom
        diff -= np.round(diff)
        cand = diff[:, None, :] + SHIFTS[None]
        dist = np.linalg.norm(cand @ cell_rows, axis=-1).min(axis=1) - radius
        best = np.minimum(best, dist)
    return best.reshape(shape)


@pytest.mark.parametrize(
    "cellpar",
    [
        (10, 12, 8, 70, 80, 60),
        (10, 11, 9, 40, 50, 45),
        (10, 12, 8, 90, 90, 30),
        (10, 12, 8, 90, 90, 150),
        (12, 12, 12, 90, 90, 90),
    ],
)
def test_dgrid_matches_exact_minimum_image(cellpar):
    rng = np.random.default_rng(3)
    cell_rows = cellpar_to_cell(cellpar)
    atoms = rng.random((2, 3)) - 0.5
    radii = np.full(len(atoms), 1.0)
    shape = (24, 24, 24)
    grid = mgr_rust.dgrid_from_atoms(atoms, np.ascontiguousarray(cell_rows.T), radii, *shape)
    expected = _reference(atoms, cell_rows, radii, shape)
    np.testing.assert_allclose(grid, expected, atol=1e-4)


def test_refinement_reports_changed_voxels_on_sparse_triclinic_cell():
    cell_rows = cellpar_to_cell((10, 12, 8, 70, 80, 60))
    atoms = np.zeros((1, 3))
    radii = np.array([1.0])
    cell = np.ascontiguousarray(cell_rows.T)
    raw, _ = mgr_rust.dgrid_from_atoms_stats(atoms, cell, radii, 30, 30, 30, False)
    fixed, stats = mgr_rust.dgrid_from_atoms_stats(atoms, cell, radii, 30, 30, 30, True)
    assert stats["changed_voxels"] > 0
    assert np.all(fixed <= raw + 1e-6)

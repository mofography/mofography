import numpy as np

from mofography.compute import check_percolation_periodic, find_pore_limiting_radius


def _channel_grid(shift):
    """Diagonal channel (winds along x and y once) with distance-like profile."""
    n = 24
    x, y, z = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    dy = np.minimum((y - x) % n, (x - y) % n)
    grid = np.clip(1.0 - 0.2 * np.hypot(dy, z - n // 2), 0, None)
    return np.roll(grid, shift, axis=(0, 1, 2))


def test_percolation_independent_of_cell_origin():
    plrs = {find_pore_limiting_radius(_channel_grid(s)) for s in [(0, 0, 0), (12, 0, 0), (7, 5, 3)]}
    assert None not in plrs
    assert len(plrs) == 1


def test_disconnected_void_does_not_percolate():
    grid = np.zeros((10, 10, 10))
    grid[3:6, 3:6, 3:6] = 1.0
    assert not check_percolation_periodic(grid, 0.5)
    assert find_pore_limiting_radius(grid) is None


def test_segmentation_raises_clear_error_when_threshold_exceeds_grid():
    import pytest
    from ase.build import bulk

    from mofography.rust_fun import dgrid_from_atoms_rust, do_periodic_segmentation_rust

    atoms = bulk("Cu", cubic=True, a=6.0)
    dgrid = dgrid_from_atoms_rust(atoms, spacing=0.5)
    with pytest.raises(ValueError, match="threshold_abs"):
        do_periodic_segmentation_rust(
            dgrid, atoms, spacing=0.5, threshold_abs=float(dgrid.max()) + 1.0, disable_tqdm=True
        )

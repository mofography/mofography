"""
Unit and regression test for the mofography package.
"""

# Import package, test suite, and other packages as needed
import sys

import numpy as np
import pytest

import mofography


def test_mofography_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "mofography" in sys.modules


def test_accessible_psd_functions_exported():
    """Verify new accessible PSD wrapper functions are exported."""
    assert hasattr(mofography, "compute_psd_global_accessible")
    assert hasattr(mofography, "compute_psds_pocketwise_accessible")
    assert callable(mofography.compute_psd_global_accessible)
    assert callable(mofography.compute_psds_pocketwise_accessible)


@pytest.mark.slow
def test_compute_psd_global_accessible_return_format():
    """Test compute_psd_global_accessible returns same format as compute_psd_uniform_rust."""
    from ase.io import read

    mat = read("mofography/examples/TestCIFs/P1-IRMOF-1.cif")
    dgrid = mofography.dgrid_from_atoms_rust(mat, spacing=0.5)
    results = mofography.do_periodic_segmentation_rust(
        dgrid, mat, spacing=0.5, return_dgrid=True, compute_connections=True
    )

    # Without fig: returns psd_hist (array)
    psd_hist = mofography.compute_psd_global_accessible(
        results, mat, probe_radius=1.0, samples_per_region=500
    )
    assert isinstance(psd_hist, np.ndarray)
    assert psd_hist.ndim == 1

    # With fig: returns (psd_hist, fig)
    psd_hist2, fig = mofography.compute_psd_global_accessible(
        results, mat, probe_radius=1.0, samples_per_region=500, return_fig=True
    )
    assert isinstance(psd_hist2, np.ndarray)
    assert fig is not None
    assert psd_hist2.shape == psd_hist.shape


@pytest.mark.slow
def test_compute_psds_pocketwise_accessible_return_format():
    """Test compute_psds_pocketwise_accessible returns same format as compute_psds_pocketwise_rust."""
    from ase.io import read

    mat = read("mofography/examples/TestCIFs/P1-IRMOF-1.cif")
    dgrid = mofography.dgrid_from_atoms_rust(mat, spacing=0.5)
    results = mofography.do_periodic_segmentation_rust(
        dgrid, mat, spacing=0.5, return_dgrid=True, compute_connections=True
    )

    # Without fig: returns pocket_psds dict
    pocket_psds = mofography.compute_psds_pocketwise_accessible(
        results, mat, probe_radius=1.0, num_samples_per_region=200
    )
    assert isinstance(pocket_psds, dict)
    for rid, data in pocket_psds.items():
        assert "hist" in data
        assert "bin_centers" in data
        assert "num_samples" in data
        assert "mean_pore_size" in data
        assert "std_pore_size" in data

    # With fig: returns (pocket_psds, fig)
    pocket_psds2, fig = mofography.compute_psds_pocketwise_accessible(
        results, mat, probe_radius=1.0, num_samples_per_region=200, return_fig=True
    )
    assert isinstance(pocket_psds2, dict)
    assert fig is not None

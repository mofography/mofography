import numpy as np
import pytest

from skimage.segmentation import watershed as sk_watershed

import mgr_rust


def _make_simple_volume(shape=(5, 5, 5)):
    z, y, x = np.indices(shape)
    center1 = np.array(shape) // 2
    center2 = center1 + np.array([0, 0, 1])
    dist1 = np.sqrt(((z - center1[0]) ** 2 + (y - center1[1]) ** 2 + (x - center1[2]) ** 2))
    dist2 = np.sqrt(((z - center2[0]) ** 2 + (y - center2[1]) ** 2 + (x - center2[2]) ** 2))
    # Take minimum distance to create two basins that meet
    dgrid = np.minimum(dist1, dist2).astype(np.float32)
    return dgrid


@pytest.mark.parametrize("backend", ["rust_det"])
def test_seeded_watershed_det_matches_skimage_small(backend):
    dgrid = _make_simple_volume()
    seeds = np.zeros_like(dgrid, dtype=np.uint32)
    seeds[2, 2, 2] = 1
    seeds[2, 2, 3] = 2

    mask = np.ones_like(dgrid, dtype=np.uint8)

    # skimage uses lower values first, but seeded_watershed_det expects +dgrid.
    sk_labels = sk_watershed(
        -dgrid,
        seeds.astype(np.int32),
        mask=mask.astype(bool),
        connectivity=2,
        compactness=0,
    ).astype(np.uint32)

    rust_labels = mgr_rust.seeded_watershed_det(
        dgrid.astype(np.float32),
        seeds.astype(np.uint32),
        mask.astype(np.uint8),
    )

    assert rust_labels.shape == sk_labels.shape
    # Direct equality since seed labels are identical
    assert np.array_equal(rust_labels, sk_labels)


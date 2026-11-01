import numpy as np

import mofography as mgr


def _poreblazer_reference(d, max_pore=20.0, binsize=0.25):
    """Direct port of PoreBlazer's pore_distribution (psd_pb_style.poreblazer_style_psd)."""
    nbins = int(max_pore / binsize)
    cum = np.zeros(nbins)
    for x in d:
        idx = min(int(x / binsize) + 1, nbins)
        cum[:idx] += 1.0
    cum /= cum[0]
    psd = np.zeros(nbins)
    for i in range(1, nbins - 1):
        psd[i] = -(cum[i + 1] - cum[i - 1]) / (2.0 * binsize)
    return (binsize * np.arange(nbins) - binsize / 2.0)[1:-1], psd[1:-1]


def test_matches_poreblazer_post_processing():
    rng = np.random.default_rng(0)
    d = np.concatenate([rng.normal(6, 0.4, 3000), rng.normal(12, 0.6, 5000), rng.uniform(0.1, 25, 200)])
    binsize, max_pore = 0.25, 20.0
    # fold the overflow into the last bin (PoreBlazer caps the bin index)
    counts, _ = np.histogram(np.minimum(d, max_pore - 1e-9), bins=np.arange(0, max_pore + binsize / 2, binsize))
    x, y = mgr.psd_gelb_gubbins(counts, binsize)
    x_ref, y_ref = _poreblazer_reference(d, max_pore, binsize)
    assert np.allclose(x, x_ref) and np.allclose(y, y_ref)

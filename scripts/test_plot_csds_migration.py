#!/usr/bin/env python3
"""Smoke tests for plot_csds / plot_psds deprecation / pocket PSD helpers."""

import warnings
import networkx as nx


def test_plot_csds_and_deprecated_plot_psds():
    import mofography as mgr

    rag = nx.Graph()
    rag.add_node(
        1,
        labels=[1],
        cluster_geom_pbc="0",
        csd_x_pbc=[1, 2, 3],
        csd_y_pbc=[0.1, 0.2, 0.1],
    )
    traces = mgr.plot_csds(rag, [{1}], cluster_key="cluster_geom_pbc")
    assert len(traces) == 1

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        mgr.plot_psds(rag, [{1}], cluster_key="cluster_geom_pbc")
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)


def test_plot_pocket_psds_on_rag():
    import mofography as mgr

    rag = nx.Graph()
    rag.add_node(
        1,
        cluster_geom_pbc="0",
        pocket_psd_x_pbc=[1.0, 2.0, 3.0],
        pocket_psd_y_pbc=[0.2, 0.5, 0.3],
    )
    traces = mgr.plot_pocket_psds_on_rag(rag, [{1}], cluster_key="cluster_geom_pbc")
    assert len(traces) == 1


if __name__ == "__main__":
    test_plot_csds_and_deprecated_plot_psds()
    test_plot_pocket_psds_on_rag()
    print("All plot_csds migration tests passed.")

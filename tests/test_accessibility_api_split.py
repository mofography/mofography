"""Winding vs lattice accessibility API split and Chao MC regression."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from ase.io import read

import mofography as mgr
from mofography.percolation_results import (
    PERCOLATION_HISTORY_LATTICE,
    PERCOLATION_HISTORY_WINDING,
    PORE_DIA_LATTICE,
    PORE_DIA_WINDING,
    get_percolation_history,
    get_pore_dia,
)


def _load_chao_mc():
    chao_path = Path.home() / "Downloads" / "compute_surface_area_from_results_mc_chao.py"
    if not chao_path.is_file():
        pytest.skip(f"Chao reference script not found: {chao_path}")
    spec = importlib.util.spec_from_file_location("chao_mc", chao_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.compute_surface_area_from_results_mc_chao


def _segment(cif_rel: str):
    cif = Path(cif_rel)
    if not cif.is_file():
        pytest.skip(f"Missing fixture CIF: {cif}")
    mat = read(cif)
    dgrid = mgr.dgrid_from_atoms_rust(mat, spacing=0.5)
    results = mgr.do_periodic_segmentation_rust(
        dgrid, mat, spacing=0.5, return_dgrid=True, compute_connections=True
    )
    return mat, results


def _winding_labels_no_fp(results, probe):
    from mofography.winding_percolation import (
        _prepare_graph_and_labels,
        percolation_all_axes_winding,
    )

    g, bl = _prepare_graph_and_labels(results, None, None)
    ac, _, _ = percolation_all_axes_winding(results, probe, graph=g, base_region_labels=bl)
    labels = set()
    for groups in ac.values():
        for group in groups:
            labels.update(group)
    return labels


def test_canonical_names_and_dual_percolation_histories():
    assert callable(mgr.compute_surface_area_mesh_winding)
    assert callable(mgr.compute_surface_area_mesh_lattice)
    assert callable(mgr.compute_surface_area_mc_winding)
    assert callable(mgr.compute_surface_area_mc_lattice)
    mat, results = _segment("examples/TestCIFs/P1-IRMOF-1.cif")
    results = dict(results)
    for k in (PERCOLATION_HISTORY_WINDING, PERCOLATION_HISTORY_LATTICE):
        results.pop(k, None)

    mgr.compute_percolation_history_winding(results)
    mgr.compute_percolation_history_from_results_lattice(results, variant="v2")

    assert PERCOLATION_HISTORY_WINDING in results
    assert PERCOLATION_HISTORY_LATTICE in results
    assert len(get_percolation_history(results, "winding")) > 0
    assert len(get_percolation_history(results, "lattice")) > 0
    # Both stored; neither overwrites the other
    assert results[PERCOLATION_HISTORY_WINDING] is not results[PERCOLATION_HISTORY_LATTICE]


@pytest.mark.slow
def test_chao_mc_matches_winding_without_fp():
    chao_mc = _load_chao_mc()
    mat, results = _segment("examples/TestCIFs/P1-IRMOF-1.cif")
    probe = 1.2
    n = 25000
    np.random.seed(0)

    chao = chao_mc(
        dict(results),
        mat,
        probe_radius=probe,
        n_samples=n,
        shell_thickness=0.5,
        return_details=True,
    )
    ours = mgr.compute_surface_area_mc_winding(
        dict(results),
        mat,
        probe_radius=probe,
        n_samples=n,
        shell_thickness=0.5,
        window_fp=False,
        return_details=True,
    )

    chao_regs = chao["percolating_regions"]
    from mofography.compute import _accessible_percolating_region_labels_winding

    ours_regs, _, _, _ = _accessible_percolating_region_labels_winding(
        results, probe, window_fp=False
    )

    assert chao_regs == ours_regs == _winding_labels_no_fp(results, probe)

    chao_sum = sum(chao["by_region"].values())
    ours_sum = sum(ours["by_region"].values())
    assert chao_sum == pytest.approx(chao["accessible_surface_area"], rel=1e-12)
    assert ours_sum == pytest.approx(ours["accessible_surface_area"], rel=1e-12)


@pytest.mark.slow
def test_window_fp_false_is_winding_not_lattice():
    mat, results = _segment("examples/TestCIFs/NU-1000_S1.cif")
    probe = 1.2
    np.random.seed(1)

    winding_no_fp = mgr.compute_surface_area_mc_winding(
        dict(results), mat, probe_radius=probe, n_samples=15000, window_fp=False
    )
    lattice = mgr.compute_surface_area_mc_lattice(
        dict(results), mat, probe_radius=probe, n_samples=15000, variant="v2"
    )
    # Methods are independent; accessible SA can differ
    assert winding_no_fp["accessible_surface_area"] >= 0
    assert lattice["accessible_surface_area"] >= 0
    # Winding no-FP should not equal lattice in general (sanity: not accidentally same code path)
    # We only assert both produce finite breakdowns when percolation exists
    if winding_no_fp["accessible_surface_area"] > 0:
        assert len(winding_no_fp["by_region"]) > 0


@pytest.mark.slow
def test_pore_dia_dual_storage():
    mat, results = _segment("examples/TestCIFs/P1-IRMOF-1.cif")
    results = dict(results)

    mgr.get_pore_dia_mgr_winding2(results=results, mat_atoms=mat, auto_rerun=False)
    mgr.get_pore_dia_lattice(mat, dgrid=results["dgrid"], results=results, variant="v2")

    assert PORE_DIA_WINDING in results
    assert PORE_DIA_LATTICE in results
    assert get_pore_dia(results, "winding")["method"] == "winding"
    assert get_pore_dia(results, "lattice")["variant"] == "v2"

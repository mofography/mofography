"""Regression: regions_pbc[rep] must match union of regions in each pbc_group."""

from __future__ import annotations

import numpy as np
import pytest
from ase.io import read

import mofography as mgr
from mofography.pocket_typing import DEFAULT_MASK_THICKNESS_A

CUBTC_CIF = "mofography/examples/TestCIFs/P1-Cu-BTC.cif"


def _pbc_group_union_matches_regions_pbc_rep(results) -> list[bool]:
    regions = np.asarray(results["regions"])
    regions_pbc = np.asarray(results["regions_pbc"])
    matches = []
    for group in results["pbc_groups"]:
        rep = min(group)
        union = np.isin(regions, sorted(int(x) for x in group))
        rep_domain = regions_pbc == rep
        matches.append(np.array_equal(union, rep_domain))
    return matches


@pytest.mark.slow
def test_regions_pbc_labels_match_pbc_group_unions_cubtc():
    atoms = read(CUBTC_CIF)
    dgrid = mgr.dgrid_from_atoms_rust(atoms, 0.5).astype(np.float32)
    results = mgr.do_periodic_segmentation_rust(
        dgrid,
        atoms,
        spacing=0.5,
        mask_thickness=DEFAULT_MASK_THICKNESS_A,
        h=0.2,
        min_distance=1.5,
        threshold_abs=1.5,
        extra_padding=5.0,
        return_dgrid=True,
        compute_connections=False,
    )
    matches = _pbc_group_union_matches_regions_pbc_rep(results)
    n_ok = sum(matches)
    assert n_ok == len(matches), (
        f"regions_pbc[rep] != union(regions[pbc_group]) for "
        f"{len(matches) - n_ok}/{len(matches)} groups on Cu-BTC"
    )

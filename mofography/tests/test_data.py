"""
Tests for mofography.data representative sampling, plotting, and CIF I/O.
"""

import numpy as np
import pandas as pd
import pytest
from pathlib import Path

import mofography as mgr


def _make_mofdb_style_df(n_rows=200):
    """Synthetic MOFdb-style DataFrame with strat columns."""
    rng = np.random.default_rng(42)
    return pd.DataFrame({
        "surface_area_m2g": rng.uniform(100, 3000, n_rows),
        "void_fraction": rng.uniform(0.2, 0.95, n_rows),
        "lcd": rng.uniform(5, 25, n_rows),
        "pld": rng.uniform(3, 15, n_rows),
        "name": [f"MOF-{i}" for i in range(n_rows)],
    })


def test_sample_representative_mofdb_returns_correct_size():
    """sample_representative_mofdb returns exactly n materials (when n < N)."""
    df = _make_mofdb_style_df(200)
    sample = mgr.sample_representative_mofdb(df, n=50, random_state=0)
    assert len(sample) == 50
    assert set(sample.columns) >= set(mgr.MOFDB_STRAT_COLS)
    assert sample.index.is_unique


def test_sample_representative_mofdb_preserves_columns():
    """Sampled rows retain all original columns."""
    df = _make_mofdb_style_df(100)
    sample = mgr.sample_representative_mofdb(df, n=20, random_state=1)
    for c in df.columns:
        assert c in sample.columns


def test_sample_representative_mofdb_n_larger_than_population():
    """When n >= N, returns all rows (no duplicate rows)."""
    df = _make_mofdb_style_df(30)
    sample = mgr.sample_representative_mofdb(df, n=50, random_state=0)
    assert len(sample) == 30
    assert sample.index.is_unique


def test_sample_representative_mofdb_missing_columns_raises():
    """Raises ValueError if strat columns are missing."""
    df = pd.DataFrame({"void_fraction": [0.5], "lcd": [10]})  # missing surface_area_m2g, pld
    with pytest.raises(ValueError, match="missing columns"):
        mgr.sample_representative_mofdb(df, n=1)


def test_plot_representative_sample_subplots_returns_figure():
    """plot_representative_sample_subplots returns a Plotly Figure with expected structure."""
    df = _make_mofdb_style_df(80)
    sample = mgr.sample_representative_mofdb(df, n=20, random_state=2)
    fig = mgr.plot_representative_sample_subplots(df, sample, n_sample=20)
    assert fig is not None
    # 4 subplots, 2 traces each (full + sample) = 8 traces
    assert len(fig.data) == 8
    assert fig.layout.title.text == "Representative subset (N=20)"


def test_plot_representative_sample_subplots_custom_title():
    """Custom title is used when provided."""
    df = _make_mofdb_style_df(50)
    sample = mgr.sample_representative_mofdb(df, n=10, random_state=3)
    fig = mgr.plot_representative_sample_subplots(
        df, sample, title="Custom title (N=10)"
    )
    assert fig.layout.title.text == "Custom title (N=10)"


# Minimal valid CIF for testing write
_MINIMAL_CIF = """
data_test
_cell_length_a 10
_cell_length_b 10
_cell_length_c 10
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Zn 0 0 0
"""


def test_write_cifs_from_mofdb_df_basic(tmp_path):
    """write_cifs_from_mofdb_df writes one .cif per row and returns paths."""
    df = pd.DataFrame({
        "name": ["A", "B", "C"],
        "cif": [_MINIMAL_CIF.strip()] * 3,
    })
    paths = mgr.write_cifs_from_mofdb_df(df, tmp_path)
    assert len(paths) == 3
    for p in paths:
        assert p.exists()
        assert p.suffix == ".cif"
        assert p.read_text(encoding="utf-8").strip().startswith("data_")
    names = {p.stem for p in paths}
    assert names == {"A", "B", "C"}


def test_write_cifs_from_mofdb_df_duplicate_names(tmp_path):
    """Duplicate names get _1, _2 suffix."""
    df = pd.DataFrame({
        "name": ["same", "same", "same"],
        "cif": [_MINIMAL_CIF] * 3,
    })
    paths = mgr.write_cifs_from_mofdb_df(df, tmp_path)
    assert len(paths) == 3
    stems = sorted(p.stem for p in paths)
    assert stems == ["same", "same_1", "same_2"]


def test_write_cifs_from_mofdb_df_skip_failed(tmp_path):
    """Rows with missing/empty CIF are skipped when skip_failed=True."""
    df = pd.DataFrame({
        "name": ["ok", "bad", "also_ok"],
        "cif": [_MINIMAL_CIF, "", _MINIMAL_CIF],
    })
    paths = mgr.write_cifs_from_mofdb_df(df, tmp_path, skip_failed=True)
    assert len(paths) == 2
    assert {p.stem for p in paths} == {"ok", "also_ok"}


def test_write_cifs_from_mofdb_df_missing_column_raises():
    """Raises if name or cif column is missing."""
    df = pd.DataFrame({"name": ["x"], "other": ["y"]})
    with pytest.raises(ValueError, match="cif"):
        mgr.write_cifs_from_mofdb_df(df, "/tmp/out", cif_column="cif")


@pytest.mark.slow
def test_write_cifs_from_mofdb_df_representative_sample_2000():
    """Write CIFs from the 2000 representative sample (requires ~/Datasets/coredf_sample_n2000.feather)."""
    sample_path = Path.home() / "Datasets" / "coredf_sample_n2000.feather"
    if not sample_path.exists():
        pytest.skip(f"Sample file not found: {sample_path}")
    df = pd.read_feather(sample_path)
    out_dir = Path.home() / "Datasets" / "coredf_sample_n2000_cifs"
    paths = mgr.write_cifs_from_mofdb_df(df, out_dir)
    assert len(paths) <= len(df)
    assert all(p.exists() for p in paths)
    # Should have one CIF per row when none are empty
    assert len(paths) >= len(df) - 10, "Most rows should yield a CIF file"

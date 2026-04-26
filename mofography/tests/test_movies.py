"""Unit tests for the mofography.movies module."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

import mofography
from mofography import movies


def test_movies_module_imported():
    """Ensure the movies module is reachable via mofography namespace."""
    assert hasattr(mofography, "MovieFigure")
    assert hasattr(mofography, "make_framework_trace")
    assert hasattr(mofography, "make_guest_trace")
    assert hasattr(mofography, "build_fixed_framework_movie")
    assert hasattr(mofography, "read_raspa_pdb_movie")


def test_moviefigure_basic_api():
    """Smoke-test the Plotly-like mutation API."""
    fig = movies.MovieFigure()
    fig.add_trace(movies._make_trace("t1", "points", "Test", {"positions": [0, 0, 0]}))
    assert len(fig.data) == 1
    assert fig.data[0]["id"] == "t1"

    fig.add_traces(
        [
            movies._make_trace("t2", "mesh", "Mesh", {"positions": [], "indices": []}),
            movies._make_trace("t3", "line_segments", "Lines", {"positions": []}),
        ]
    )
    assert len(fig.data) == 3

    fig.update_traces({"visible": False}, selector={"type": "points"})
    assert fig.data[0]["visible"] is False
    assert fig.data[1]["visible"] is True

    fig.set_camera(eye=(1, 2, 3), center=(0, 0, 0), up=(0, 0, 1))
    assert fig.layout["scene"]["camera"]["eye"] == [1, 2, 3]

    d = fig.to_dict()
    assert d["version"] == "1.0"
    assert len(d["static"]["traces"]) == 3
    assert "styles" in d
    assert "groups" in d


def test_moviefigure_from_dict_roundtrip():
    """Round-trip serialization."""
    fig = movies.MovieFigure()
    fig.add_trace(movies._make_trace("x", "points", "X", {"positions": [1, 2, 3]}))
    d = fig.to_dict()
    fig2 = movies.MovieFigure.from_dict(d)
    assert len(fig2.data) == 1
    assert fig2.data[0]["id"] == "x"


def test_default_styles_registry():
    """Style registry keys must cover all trace types we emit."""
    required = {
        "framework",
        "guests",
        "surface",
        "cell",
        "internal_maxima",
        "periodic_maxima",
        "internal_windows",
        "periodic_windows",
        "internal_connections",
        "periodic_connections",
        "internal_nodes",
        "periodic_nodes",
    }
    assert required.issubset(movies.DEFAULT_STYLES.keys())


def test_trace_schema_validation():
    """Every trace must have id, type, name, data."""
    trace = movies._make_trace("abc", "mesh", "Framework", {"positions": []})
    assert trace["id"] == "abc"
    assert trace["type"] == "mesh"
    assert "data" in trace
    assert "style" in trace


def test_serialize_trimesh_empty():
    """Empty mesh produces valid payload with zero counts."""
    import trimesh

    mesh = trimesh.Trimesh(vertices=[], faces=[])
    payload = movies._serialize_trimesh(mesh)
    assert payload["n_vertices"] == 0
    assert payload["n_triangles"] == 0
    assert payload["positions"] == []
    assert payload["indices"] == []


def test_serialize_cell_lines_shape():
    """Cell line serialization must have 12 edges -> 24 points."""
    from ase import Atoms

    atoms = Atoms("H", positions=[[0, 0, 0]], cell=[5, 5, 5], pbc=True)
    payload = movies._serialize_cell_lines(atoms)
    assert payload["n_segments"] == 12
    assert len(payload["positions"]) == 12 * 2 * 3  # 12 edges * 2 pts * 3 coords


def test_make_guest_trace_shape():
    """Guest trace data must match atoms count."""
    from ase import Atoms

    atoms = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]], cell=[10, 10, 10], pbc=True)
    trace = movies.make_guest_trace(atoms)
    assert trace["type"] == "instanced_spheres"
    assert trace["data"]["n_atoms"] == 2
    assert len(trace["data"]["positions"]) == 2 * 3


def test_read_raspa_pdb_movie_not_found():
    """Graceful error on missing file."""
    with pytest.raises(FileNotFoundError):
        movies.read_raspa_pdb_movie("/nonexistent/path.pdb")


def test_make_movie_from_raspa_pdb_missing_files():
    """High-level helper must raise on missing inputs."""
    with pytest.raises(FileNotFoundError):
        movies.make_movie_from_raspa_pdb("/nonexistent.pdb", "/nonexistent.cif")


def test_to_html_writes_file(tmp_path: Path) -> None:
    """to_html must produce a valid HTML file."""
    fig = movies.MovieFigure()
    fig.add_trace(movies._make_trace("t1", "points", "T1", {"positions": [0, 0, 0]}))
    out = tmp_path / "test_movie.html"
    result = fig.to_html(str(out))
    assert result.exists()
    text = result.read_text(encoding="utf-8")
    assert "<!doctype html>" in text
    assert "three@0.128.0" in text
    assert "MovieFigure" not in text  # ensure template substitution happened

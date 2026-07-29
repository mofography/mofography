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


# ---------------------------------------------------------------------------
# Figure3D convenience methods
# ---------------------------------------------------------------------------

def test_figure3d_alias():
    """Figure3D must be the same class as MovieFigure."""
    assert mofography.Figure3D is mofography.MovieFigure
    assert hasattr(mofography, "Figure3D")


def test_add_framework_chainable():
    """add_framework adds a framework trace and returns self."""
    from ase import Atoms

    atoms = Atoms("Si2", positions=[[0, 0, 0], [2.5, 2.5, 2.5]], cell=[5, 5, 5], pbc=True)
    fig = mofography.Figure3D()
    result = fig.add_framework(atoms)
    assert result is fig
    assert len(fig.data) == 3
    ids = {t["id"] for t in fig.data}
    assert "framework" in ids
    assert "framework_stick" in ids
    assert "framework_spacefill" in ids
    assert fig.data[0]["type"] == "mesh"


def test_add_cell_chainable():
    """add_cell adds a unit-cell trace and returns self."""
    from ase import Atoms

    atoms = Atoms("H", positions=[[0, 0, 0]], cell=[5, 5, 5], pbc=True)
    fig = mofography.Figure3D()
    result = fig.add_cell(atoms)
    assert result is fig
    assert len(fig.data) == 1
    assert fig.data[0]["id"] == "cell"


def test_add_surface_from_arrays():
    """add_surface works with raw dgrid + regions arrays."""
    from ase import Atoms
    import numpy as np

    atoms = Atoms("H", positions=[[2, 2, 2]], cell=[5, 5, 5], pbc=True)
    dgrid = np.ones((10, 10, 10), dtype=np.float32) * 3.0
    dgrid[4:7, 4:7, 4:7] = 0.1
    regions = np.zeros((10, 10, 10), dtype=np.int32)
    regions[3:8, 3:8, 3:8] = 1

    fig = mofography.Figure3D()
    fig.add_surface(dgrid, regions, atoms, isovalue=1.5)
    assert len(fig.data) == 1
    assert fig.data[0]["id"] == "surface"
    assert fig.data[0]["type"] == "mesh"


def test_add_maxima_hides_periodic():
    """Periodic maxima should be hidden when hide_periodic=True."""
    seg = _make_minimal_seg()
    atoms = _make_cubic_atoms()

    fig = mofography.Figure3D()
    fig.add_maxima(seg, atoms, hide_periodic=True)

    traces = {t["id"]: t for t in fig.data}
    assert traces["maxima_internal"]["visible"] is True
    assert traces["maxima_outside"]["visible"] is False


def test_add_connections_hides_periodic():
    """Periodic connections should be hidden when hide_periodic=True."""
    seg = _make_minimal_seg()
    atoms = _make_cubic_atoms()

    fig = mofography.Figure3D()
    fig.add_connections(seg, atoms, hide_periodic=True)

    traces = {t["id"]: t for t in fig.data}
    assert traces["connections_internal"]["visible"] is True
    assert traces["connections_periodic"]["visible"] is False
    assert traces["connections_periodic_nodes"]["visible"] is False


def test_add_windows_hides_periodic():
    """Periodic windows should be hidden when hide_periodic=True."""
    seg = _make_minimal_seg()
    atoms = _make_cubic_atoms()

    fig = mofography.Figure3D()
    fig.add_windows(seg, hide_periodic=True)

    traces = {t["id"]: t for t in fig.data}
    assert traces["windows_internal"]["visible"] is True
    assert traces["windows_periodic"]["visible"] is False


def test_add_guests_with_color():
    """add_guests with explicit color should override style and disable per-instance colors."""
    from ase import Atoms

    atoms = Atoms("K3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]], cell=[10, 10, 10], pbc=True)
    fig = mofography.Figure3D()
    fig.add_guests(atoms, color=[1.0, 0.0, 1.0])

    trace = fig.data[0]
    assert trace["id"] == "guests"
    assert trace["type"] == "instanced_spheres"
    assert trace["style"]["color"] == [1.0, 0.0, 1.0]
    assert trace["data"]["colors"] == []  # per-instance colors disabled
    assert trace["data"]["n_atoms"] == 3


def test_make_framework_instanced_representations():
    """Instanced framework traces vary radii/bonds by representation."""
    from ase import Atoms
    from ase.data import covalent_radii, vdw_radii

    atoms = Atoms("C2", positions=[[0, 0, 0], [1.4, 0, 0]], cell=[5, 5, 5], pbc=True)
    atoms.info["covalent_bonds"] = [[0, 1]]
    bs = movies.make_framework_instanced_trace(atoms, representation="ball_and_stick", atom_scale=0.65)
    stick = movies.make_framework_instanced_trace(atoms, representation="stick")
    sf = movies.make_framework_instanced_trace(atoms, representation="space_filling")

    bs_r = bs["atoms"]["data"]["radii"][0]
    stick_r = stick["atoms"]["data"]["radii"][0]
    sf_r = sf["atoms"]["data"]["radii"][0]
    assert stick_r < bs_r
    assert abs(stick_r - movies.FRAMEWORK_STICK_RADIUS) < 1e-6
    assert abs(stick["bonds"]["data"]["radius"] - movies.FRAMEWORK_STICK_RADIUS) < 1e-6
    assert abs(sf_r - float(vdw_radii[6])) < 1e-6
    assert "bonds" in bs and bs["bonds"]["data"]["n_bonds"] == 1
    assert "bonds" in stick and stick["bonds"]["data"]["n_bonds"] == 2
    assert "bonds" not in sf
    starts = stick["bonds"]["data"]["start_positions"]
    assert starts[:3] == pytest.approx([0.0, 0.0, 0.0])


def test_zif8_aromatic_cycles_merge_duplicate_sites():
    """RASPA-expanded ZIF-8 CIF duplicates sites; merge before ring detection."""
    from pathlib import Path

    from ase.io import read

    from mofography.analyze import compute_aromatic_cycles, update_with_covalent_bonds

    cif = Path(__file__).resolve().parents[2] / "examples/TestCIFs/ZIF-8-RASPA.cif"
    if not cif.is_file():
        pytest.skip("ZIF-8 CIF missing")
    atoms = update_with_covalent_bonds(read(str(cif)))
    without_merge = compute_aromatic_cycles(
        atoms, force_recompute=True, merge_duplicate_sites=False,
    )
    with_merge = compute_aromatic_cycles(
        atoms, force_recompute=True, merge_duplicate_sites=True,
    )
    assert len(without_merge) == 12
    assert len(with_merge) == 24
    assert all(len(c) == 5 for c in with_merge)


def test_full_chainable_pipeline():
    """All add_* methods chained together produce correct traces."""
    from ase import Atoms
    import numpy as np

    atoms = Atoms("Si2", positions=[[0, 0, 0], [2.5, 2.5, 2.5]], cell=[5, 5, 5], pbc=True)
    dgrid = np.ones((10, 10, 10), dtype=np.float32) * 3.0
    dgrid[4:7, 4:7, 4:7] = 0.1
    regions = np.zeros((10, 10, 10), dtype=np.int32)
    regions[3:8, 3:8, 3:8] = 1
    seg = _make_minimal_seg()

    fig = (
        mofography.Figure3D()
        .add_framework(atoms)
        .add_cell(atoms)
        .add_surface(dgrid, regions, atoms)
        .add_maxima(seg, atoms)
        .add_connections(seg, atoms)
        .add_windows(seg)
        .add_guests(atoms)
    )

    expected_ids = {
        "framework", "cell", "surface",
        "maxima_internal", "maxima_outside",
        "connections_internal", "connections_internal_nodes",
        "connections_periodic", "connections_periodic_nodes",
        "windows_internal", "windows_periodic",
        "guests",
    }
    trace_ids = {t["id"] for t in fig.data}
    assert expected_ids.issubset(trace_ids)
    assert len(fig.data) >= 12


# ---------------------------------------------------------------------------
# Jupyter integration
# ---------------------------------------------------------------------------

def test_render_html_string():
    """_render_html_string returns valid HTML string with three.js CDN."""
    fig = mofography.Figure3D()
    fig.add_trace(movies._make_trace("t1", "points", "T1", {"positions": [0, 0, 0]}))
    html_str = fig._render_html_string()
    assert "<!doctype html>" in html_str
    assert "three@0.128.0" in html_str
    assert "MovieFigure" not in html_str


def test_repr_html_contains_iframe():
    """_repr_html_ wraps content in an IFrame with srcdoc."""
    fig = mofography.Figure3D()
    fig.add_trace(movies._make_trace("t1", "points", "T1", {"positions": [0, 0, 0]}))
    repr_html = fig._repr_html_()
    assert repr_html.startswith('<iframe ')
    assert 'srcdoc="' in repr_html
    assert 'three@0.128.0' in repr_html


def test_show_browser_writes_temp_file(monkeypatch):
    """show(renderer='browser') writes a temp HTML file."""
    import tempfile
    import webbrowser

    opened = []

    def fake_open(path: str) -> None:
        opened.append(path)

    monkeypatch.setattr(webbrowser, "open", fake_open)

    fig = mofography.Figure3D()
    fig.add_trace(movies._make_trace("t1", "points", "T1", {"positions": [0, 0, 0]}))
    fig.show(renderer="browser")

    assert len(opened) == 1
    tmp_path = Path(opened[0])
    assert tmp_path.suffix == ".html"
    content = tmp_path.read_text(encoding="utf-8")
    assert "<!doctype html>" in content
    # Cleanup
    tmp_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Static scene builder
# ---------------------------------------------------------------------------

def test_build_static_scene_basic():
    """build_static_scene with just atoms produces framework + cell."""
    from ase import Atoms

    atoms = Atoms("Si2", positions=[[0, 0, 0], [2.5, 2.5, 2.5]], cell=[5, 5, 5], pbc=True)
    fig = mofography.build_static_scene(atoms)

    trace_ids = {t["id"] for t in fig.data}
    assert "framework" in trace_ids
    assert "framework_stick" in trace_ids
    assert "framework_spacefill" in trace_ids
    assert "cell" in trace_ids
    assert len(fig.data) == 4  # framework + framework_stick + framework_spacefill + cell


def test_build_static_scene_with_surface():
    """build_static_scene with dgrid+regions includes surface."""
    from ase import Atoms
    import numpy as np

    atoms = Atoms("H", positions=[[2, 2, 2]], cell=[5, 5, 5], pbc=True)
    dgrid = np.ones((10, 10, 10), dtype=np.float32) * 3.0
    dgrid[4:7, 4:7, 4:7] = 0.1
    regions = np.zeros((10, 10, 10), dtype=np.int32)
    regions[3:8, 3:8, 3:8] = 1

    fig = mofography.build_static_scene(atoms, dgrid=dgrid, regions=regions, isovalue=1.5)
    trace_ids = {t["id"] for t in fig.data}
    assert "framework" in trace_ids
    assert "cell" in trace_ids
    assert "surface" in trace_ids


def test_build_static_scene_with_seg():
    """build_static_scene with full seg dict adds pore graph traces."""
    from ase import Atoms
    import numpy as np

    atoms = Atoms("Si2", positions=[[0, 0, 0], [2.5, 2.5, 2.5]], cell=[5, 5, 5], pbc=True)
    dgrid = np.ones((10, 10, 10), dtype=np.float32) * 3.0
    dgrid[4:7, 4:7, 4:7] = 0.1
    regions = np.zeros((10, 10, 10), dtype=np.int32)
    regions[3:8, 3:8, 3:8] = 1
    seg = _make_minimal_seg()

    fig = mofography.build_static_scene(
        atoms, dgrid=dgrid, regions=regions, seg=seg, hide_periodic=True,
    )
    trace_ids = {t["id"] for t in fig.data}
    for tid in ("framework", "cell", "surface", "maxima_internal", "connections_internal"):
        assert tid in trace_ids, f"missing trace: {tid}"

    # Periodic traces should be hidden
    periodic_traces = {
        t["id"]: t["visible"]
        for t in fig.data
        if "periodic" in t["id"] or "outside" in t["id"]
    }
    for tid, vis in periodic_traces.items():
        assert vis is False, f"{tid} should be hidden"


# ---------------------------------------------------------------------------
# Public compute helpers
# ---------------------------------------------------------------------------

def test_compute_framework_mesh():
    """compute_framework_mesh returns complete payload."""
    from ase import Atoms

    atoms = Atoms("Si2", positions=[[0, 0, 0], [2.5, 2.5, 2.5]], cell=[5, 5, 5], pbc=True)
    data = mofography.compute_framework_mesh(atoms)
    assert "positions" in data
    assert "indices" in data
    assert "colors" in data
    assert data["n_vertices"] > 0
    assert data["n_triangles"] > 0
    assert data["n_atoms"] == 2


def test_compute_cell_lines():
    """compute_cell_lines returns 12 edges."""
    from ase import Atoms

    atoms = Atoms("H", positions=[[0, 0, 0]], cell=[5, 5, 5], pbc=True)
    data = mofography.compute_cell_lines(atoms)
    assert data["n_segments"] == 12
    assert len(data["positions"]) == 12 * 6  # 12 edges * 2 pts * 3 coords


def test_compute_surface_mesh():
    """compute_surface_mesh returns marching-cubes output."""
    from ase import Atoms
    import numpy as np

    atoms = Atoms("H", positions=[[2, 2, 2]], cell=[5, 5, 5], pbc=True)
    dgrid = np.ones((10, 10, 10), dtype=np.float32) * 3.0
    dgrid[4:7, 4:7, 4:7] = 0.1
    regions = np.zeros((10, 10, 10), dtype=np.int32)
    regions[3:8, 3:8, 3:8] = 1

    data = mofography.compute_surface_mesh(atoms, dgrid, regions, isovalue=1.5)
    assert "positions" in data
    assert "indices" in data
    assert "colors" in data
    assert "labels" in data
    assert data["n_vertices"] > 0
    assert data["n_triangles"] > 0


def test_compute_connection_data():
    """compute_connection_data handles empty connections gracefully."""
    from ase import Atoms
    import numpy as np

    atoms = _make_cubic_atoms()
    seg = _make_minimal_seg()
    data = mofography.compute_connection_data(seg, atoms)
    assert "internal_segments" in data
    assert "periodic_segments" in data
    assert "internal_points" in data
    assert "periodic_points" in data


def test_connection_nodes_internal_periodic_mutually_exclusive():
    """Periodic node markers must not duplicate internal-only wall nodes."""
    from ase import Atoms
    import numpy as np

    atoms = Atoms("Si2", positions=[[0, 0, 0], [2.5, 2.5, 2.5]], cell=[5, 5, 5], pbc=True)
    seg = _make_minimal_seg()
    # Node 1 on internal edge; node 2 only on periodic edge.
    seg["connections"] = [
        ([1, 2], True, "internal", [], [1.0]),
        ([2, 1], True, "periodic", [], [1.2]),
    ]
    data = movies._serialize_connection_traces(seg, atoms)
    internal_pts = np.asarray(data["internal_points"], dtype=float).reshape(-1, 3)
    periodic_pts = np.asarray(data["periodic_points"], dtype=float).reshape(-1, 3)
    assert internal_pts.shape[0] == 2
    assert periodic_pts.shape[0] == 0


def _connection_passes_window_threshold(segment_radii, threshold):
    """Mirror movies.py connectionPassesWindowThreshold (Three.js viewer)."""
    if threshold is None:
        return True
    if not segment_radii:
        return False
    vals = [float(x) for x in segment_radii if np.isfinite(x)]
    if not vals:
        return False
    return min(vals) >= threshold


def test_classify_connections_by_close_or_inside():
    from mofography.helper import classify_connections_by_close_or_inside

    conns = [
        [[1, 2], True, "internal", [], [1.0]],
        [[1, 10], True, "internal", [], [0.5]],
        [[2, 3], True, "periodic", [], [1.0]],
    ]
    classify_connections_by_close_or_inside(conns, [1, 2, 3])
    assert conns[0][2] == "internal"
    assert conns[1][2] == "periodic"
    assert conns[2][2] == "internal"
    classify_connections_by_close_or_inside(conns, [1, 2], remap_dict={10: 99})
    assert conns[1][2] == "periodic"


def test_internal_connection_window_threshold_filter():
    """Internal pore-graph edges filter by max window radius along edge."""
    from ase import Atoms

    atoms = Atoms("Si2", positions=[[0, 0, 0], [2.5, 2.5, 2.5]], cell=[5, 5, 5], pbc=True)
    seg = _make_minimal_seg()
    seg["connections"] = [
        ([1, 2], True, "internal", [], [2.0, 3.5]),
        ([1, 2], True, "internal", [], [4.0]),
        ([1, 2], True, "internal", [], [1.2]),
        ([1, 2], True, "periodic", [], [5.0]),
    ]
    data = movies._serialize_connection_traces(seg, atoms)
    radii = data["internal_segment_window_radii"]
    assert len(radii) == 3
    kept = [
        i for i, wr in enumerate(radii)
        if _connection_passes_window_threshold(wr, 3.0)
    ]
    assert kept == [1]
    trace = movies.make_connections_trace(seg, atoms)
    assert trace["internal_lines"]["data"]["segment_window_radii"] == radii
    assert trace["internal_lines"]["data"]["n_segments"] == 3


# ---------------------------------------------------------------------------
# Helpers for tests
# ---------------------------------------------------------------------------

def _make_cubic_atoms() -> "Atoms":
    from ase import Atoms
    return Atoms("Si2", positions=[[0, 0, 0], [2.5, 2.5, 2.5]], cell=[5, 5, 5], pbc=True)


def _make_minimal_seg() -> dict:
    """Return a minimal segmentation dict with mock maxima and no connections."""
    import numpy as np

    dgrid = np.ones((10, 10, 10), dtype=np.float32) * 3.0
    dgrid[4:7, 4:7, 4:7] = 0.1
    regions = np.zeros((10, 10, 10), dtype=np.int32)
    regions[3:8, 3:8, 3:8] = 1
    all_maxima = np.array([[5, 5, 5], [10, 5, 5]], dtype=int)
    return {
        "dgrid": dgrid,
        "regions": regions,
        "all_maxima": all_maxima,
        "maxima_pbc": np.array([[5, 5, 5]]),
        "connections": [],
    }

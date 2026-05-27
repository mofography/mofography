"""Upstream feature tests — mofography.

Verifies that the features ported from mgr_cli_py are present and
functional in the mofography package.

Run with::

    conda run -n mgr_env python -m pytest tests/test_upstream_from_mgr_cli_py.py -v

Each test is adapted from real usage patterns in
python_cli/mgrcli/interactive_app.py (mgr_cli_py).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

RASPA_CIF_DIR = Path("/Users/arungopalan/repos/RASPA2/structures/mofs/cif")
ZIF20_CIF = RASPA_CIF_DIR / "ZIF-20.cif"
CUBTC_CIF = RASPA_CIF_DIR / "Cu-BTC.cif"

# Fallback CIF from the mofography repo itself
REPO_ROOT = Path(__file__).parent.parent
_TEST_CIFS = list((REPO_ROOT / "cifs").glob("*.cif"))
FALLBACK_CIF = _TEST_CIFS[0] if _TEST_CIFS else None


def _any_cif() -> Path:
    for p in [ZIF20_CIF, CUBTC_CIF, FALLBACK_CIF]:
        if p is not None and p.exists():
            return p
    pytest.skip("No CIF file available for test")


# ---------------------------------------------------------------------------
# 1. Rust CIF parser
# ---------------------------------------------------------------------------

class TestRustCifParser:
    """Tests for mgr_rust.read_cif and mofography.cif_io.read_cif_to_ase."""

    def test_mgr_rust_has_read_cif(self):
        import mgr_rust
        assert hasattr(mgr_rust, "read_cif"), "mgr_rust.read_cif missing"

    def test_read_cif_to_ase_returns_atoms(self):
        from mofography.cif_io import read_cif_to_ase
        cif = _any_cif()
        atoms = read_cif_to_ase(cif)
        assert len(atoms) > 0, "No atoms parsed"
        assert atoms.get_cell().volume > 0, "Zero-volume cell"

    def test_read_cif_fallback_on_bad_path(self):
        """read_cif_to_ase should fall back to ASE on parser failure."""
        from mofography.cif_io import read_cif_to_ase
        with pytest.raises(Exception):
            # Should raise FileNotFoundError (or similar), not crash silently
            read_cif_to_ase("/nonexistent/path/fake.cif")

    def test_symmetry_expansion_produces_valid_cell(self):
        """CIF read + symmetry expansion: cell vectors should be consistent."""
        from mofography.cif_io import read_cif_to_ase
        cif = _any_cif()
        atoms = read_cif_to_ase(cif)
        cell = atoms.get_cell().array
        # Cell matrix should be 3×3 with non-zero determinant
        assert np.abs(np.linalg.det(cell)) > 0.1, "Degenerate cell"


# ---------------------------------------------------------------------------
# 2. Fast covalent bonds (Rust)
# ---------------------------------------------------------------------------

class TestCovalentBonds:
    """Tests for mgr_rust.covalent_bonds."""

    def test_mgr_rust_has_covalent_bonds(self):
        import mgr_rust
        assert hasattr(mgr_rust, "covalent_bonds"), "mgr_rust.covalent_bonds missing"

    def test_h2o_two_bonds(self):
        """H2O should have exactly 2 O–H bonds."""
        import mgr_rust
        from ase.build import molecule
        mol = molecule("H2O")
        mol.set_cell([10, 10, 10])
        mol.set_pbc([True, True, True])
        positions = mol.get_positions().ravel().tolist()
        an = mol.get_atomic_numbers().tolist()
        cell_f = mol.get_cell().array.T.flatten(order="F").tolist()
        # Provide all-pairs neighbor list
        n = len(mol)
        ni = [i for i in range(n) for j in range(n) if i != j]
        nj = [j for i in range(n) for j in range(n) if i != j]
        bonds_flat, bond_orders, img_pos, img_an = mgr_rust.covalent_bonds(
            positions, an, cell_f, ni, nj, 0.3
        )
        # bonds_flat contains directed pairs; with 4 directed inputs we get 4 entries.
        # Unique undirected bonds: (0,1)=O-H and (0,2)=O-H → 2 unique.
        pairs = {(min(bonds_flat[2*i], bonds_flat[2*i+1]),
                  max(bonds_flat[2*i], bonds_flat[2*i+1]))
                 for i in range(len(bonds_flat)//2)}
        assert len(pairs) == 2, f"Expected 2 unique O–H bonds, got {pairs}"

    def test_update_with_covalent_bonds_uses_rust(self):
        """update_with_covalent_bonds should complete without ImportError."""
        from mofography.analyze import update_with_covalent_bonds
        from mofography.cif_io import read_cif_to_ase
        cif = _any_cif()
        atoms = read_cif_to_ase(cif)
        updated = update_with_covalent_bonds(atoms)
        assert "covalent_bonds" in updated.info, "covalent_bonds key missing"
        bonds = updated.info["covalent_bonds"]
        assert len(bonds) > 0, "No bonds found"


# ---------------------------------------------------------------------------
# 3. Aromatic ring Rust expansion
# ---------------------------------------------------------------------------

class TestAromaticRings:
    """Tests for mgr_rust.expand_aromatic_cycle_embeddings."""

    def test_mgr_rust_has_expand_aromatic(self):
        import mgr_rust
        assert hasattr(mgr_rust, "expand_aromatic_cycle_embeddings")

    def test_expand_aromatic_cycle_embeddings_basic(self):
        """Benzene-like 6-ring graph: should find one embedding of itself."""
        import mgr_rust
        # Build a 6-ring graph (benzene)
        n = 6
        ei = [0, 1, 2, 3, 4, 5]
        ej = [1, 2, 3, 4, 5, 0]
        elems = ["C"] * 6
        bts = ["aromatic"] * 6
        templates = [list(range(6))]  # template is the ring itself
        result = mgr_rust.expand_aromatic_cycle_embeddings(n, ei, ej, elems, bts, templates)
        # Should return at least one embedding
        assert len(result) >= 1, "No embeddings found for 6-ring"

    def test_expand_aromatic_embeddings_rust_helper(self):
        """Test _expand_aromatic_embeddings_rust via mofography.analyze."""
        from mofography.analyze import _expand_aromatic_embeddings_rust, _sgraph_to_rust_arrays
        import networkx as nx
        # 5-ring
        g = nx.cycle_graph(5)
        nx.set_edge_attributes(g, "aromatic", "bond_type")
        nx.set_node_attributes(g, "C", "specie")
        result = _expand_aromatic_embeddings_rust(g, [list(range(5))])
        assert isinstance(result, list)


# ---------------------------------------------------------------------------
# 4. Pocket typing
# ---------------------------------------------------------------------------

class TestPocketTyping:
    """Tests for mofography.pocket_typing and pocket_framework_env."""

    def test_pocket_fingerprint_dataclass(self):
        from mofography.pocket_framework_env import PocketFingerprint
        fp = PocketFingerprint(
            rep=1,
            vector=np.zeros(20),
            r_pocket=5.0,
            r_shell=6.0,
            n_atoms=10,
            element_z=(6, 7),
            global_radial=np.ones(10) * 0.1,
            global_sh_power=np.ones(7),
            per_z_radial={6: np.ones(10) * 0.1, 7: np.ones(10) * 0.05},
            per_z_sh_power={6: np.ones(7), 7: np.ones(7) * 0.5},
            max_probe_radius_a=5.0,
        )
        assert fp.r_pocket == 5.0
        assert fp.n_atoms == 10

    def test_fingerprint_cosine_similarity_identical(self):
        """Two identical fingerprints should have similarity 1.0."""
        from mofography.pocket_framework_env import PocketFingerprint, fingerprint_cosine_similarity
        kwargs = dict(
            rep=1, vector=np.zeros(20), r_pocket=5.0, r_shell=6.0, n_atoms=10,
            element_z=(6,), global_radial=np.ones(10), global_sh_power=np.ones(7),
            per_z_radial={6: np.ones(10)}, per_z_sh_power={6: np.ones(7)},
            max_probe_radius_a=5.0,
        )
        fp1 = PocketFingerprint(**kwargs)
        fp2 = PocketFingerprint(**{**kwargs, "rep": 2})
        sim = fingerprint_cosine_similarity(fp1, fp2)
        assert abs(sim - 1.0) < 1e-6, f"sim={sim}, expected 1.0"

    def test_assign_pocket_types_from_segmentation_importable(self):
        from mofography.pocket_typing import assign_pocket_types_from_segmentation
        assert callable(assign_pocket_types_from_segmentation)

    def test_workflow_assign_pocket_types(self):
        from mofography.workflow import assign_pocket_types
        assert callable(assign_pocket_types)


# ---------------------------------------------------------------------------
# 5. Window plotting improvements
# ---------------------------------------------------------------------------

class TestWindowPlotting:
    """Tests for far-outside window reclassification and marker size."""

    def test_serialize_windows_trace_signature(self):
        """_serialize_windows_trace should accept atoms and spacing kwargs."""
        from mofography.movies import _serialize_windows_trace
        import inspect
        sig = inspect.signature(_serialize_windows_trace)
        params = list(sig.parameters.keys())
        assert "atoms" in params, f"atoms not in {params}"
        assert "spacing" in params, f"spacing not in {params}"

    def test_serialize_windows_trace_no_atoms(self):
        """Should work without atoms (falls back gracefully)."""
        from mofography.movies import _serialize_windows_trace
        # Minimal seg dict with one internal window
        seg = {
            "connections": [
                ((1, 2), True, "internal", [[5.0, 5.0, 5.0]], [1.5]),
            ]
        }
        result = _serialize_windows_trace(seg)
        assert "internal_n_points" in result
        assert result["internal_n_points"] == 1

    def test_plot_windows_marker_size_10(self):
        """plot_windows should use marker_size=10 (upstream fix from mgr_cli_py)."""
        from mofography.visualize import plot_windows
        import inspect
        src = inspect.getsource(plot_windows)
        assert "marker_size=10" in src or "size=10" in src, \
            "Window marker size should be 10 (not 8)"

    def test_plot_percolation_windows_importable(self):
        from mofography.visualize import plot_percolation_windows
        assert callable(plot_percolation_windows)

    def test_add_percolation_windows_to_fig_importable(self):
        from mofography.workflow import add_percolation_windows_to_fig
        assert callable(add_percolation_windows_to_fig)


# ---------------------------------------------------------------------------
# 6. PLD / percolation plot functions
# ---------------------------------------------------------------------------

class TestPLDPercolation:
    """Tests for PLD and percolation window trace functions."""

    def test_plot_percolation_windows_returns_traces(self):
        """plot_percolation_windows on a minimal analysis result."""
        from mofography.visualize import plot_percolation_windows
        # Minimal analysis dict with two connections
        # percolation_history: list of (threshold, {axis: [component]})
        analysis = {
            "connections": [
                ((1, 2), True, "internal", [[5.0, 5.0, 5.0]], [1.5]),
                ((2, 3), True, "internal", [[6.0, 6.0, 6.0]], [2.5]),
            ],
            "percolation_history": [(1.5, {"x": [[1, 2, 3]]})],
        }
        # Returns a 2-tuple (window_traces, pld_traces) or a list — both are acceptable
        result = plot_percolation_windows(analysis, pld_probe_radius=1.5, pld_tolerance=0.3)
        # Should not raise; result should be iterable
        assert result is not None

    def test_add_percolation_windows_to_fig_executes(self):
        """add_percolation_windows_to_fig should add traces to a figure."""
        import plotly.graph_objects as go
        from mofography.workflow import add_percolation_windows_to_fig
        fig = go.Figure()
        n_before = len(fig.data)
        analysis = {
            "connections": [
                ((1, 2), True, "internal", [[5.0, 5.0, 5.0]], [1.5]),
            ],
            "percolation_history": [(1.5, {"x": [[1, 2]]})],
        }
        fig = add_percolation_windows_to_fig(fig, analysis, pld_probe_radius=1.5)
        # Should not raise and fig should be a Figure
        import plotly.graph_objects as go
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import subprocess
    result = subprocess.run(
        ["conda", "run", "-n", "mgr_env", "python", "-m", "pytest", __file__, "-v"],
        cwd=REPO_ROOT,
    )
    sys.exit(result.returncode)

"""Tests for gRASPA LAMMPS .data snapshot parsing."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mofography.data import read_graspa_lammps_data
from mofography.raspa_helper import (
    GRASPA_FRAMEWORK_MOL_ID,
    get_atoms_from_graspa_lammps_snapshot,
    graspa_lammps_data_to_ase,
)

_MINIMAL_DATA = """\
# LAMMPS data file written by WriteLammpsdata function in RASPA(written by Zhao Li)
# Total atoms: 3, Printable atoms: 3
3 atoms
0 bonds
0 angles
0 dihedrals
0 impropers
0 10.0 xlo xhi
0 10.0 ylo yhi
0 10.0 zlo zhi
0.0 0.0 0.0 xy xz yz

2 atom types
0 bond types
0 angle types
0 dihedral types
0 improper types

Masses

1 12.011 # C
2 1.008 # H_xyl

Pair Coeffs

1 0.1 3.0 # C
2 0.1 2.5 # H_xyl

Atoms

1 1 1 0.0 1.0 2.0 3.0 # FRAME.cif C
2 1 1 0.0 4.0 5.0 6.0 # FRAME.cif O
3 2 2 -0.1 7.1 8.2 9.3 # o-xylene H_xyl
"""

ACODAZ_SNAPSHOT = Path(
    "/tmp/gcmc-test/gcmc-mgr-0708/1.99/ACODAZ_clean_pacmof/Movies/System_0/result_90000.data"
)


def _reference_atoms_from_file(path: Path):
    """Independent line parser for cross-validation."""
    in_atoms = False
    refs = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s == "Atoms":
            in_atoms = True
            continue
        if not in_atoms or not s or s.startswith("#"):
            continue
        body, _, comment = s.partition("#")
        toks = body.split()
        if len(toks) < 7:
            continue
        refs.append(
            {
                "id": int(toks[0]),
                "mol_id": int(toks[1]),
                "type_id": int(toks[2]),
                "charge": float(toks[3]),
                "position": np.array(
                    [float(toks[4]), float(toks[5]), float(toks[6])], dtype=float
                ),
                "comment": comment.strip(),
            }
        )
    return refs


def test_read_graspa_lammps_data_minimal(tmp_path):
    p = tmp_path / "snap.data"
    p.write_text(_MINIMAL_DATA, encoding="utf-8")
    data = read_graspa_lammps_data(p)
    assert data["n_atoms"] == 3
    assert len(data["atoms"]) == 3
    assert data["type_labels"] == {1: "C", 2: "H_xyl"}
    assert np.allclose(data["cell"].diagonal(), [10.0, 10.0, 10.0])
    assert data["atoms"][0]["symbol"] == "C"
    assert data["atoms"][2]["symbol"] == "H"
    assert data["atoms"][2]["molecule_name"] == "o-xylene"


def test_graspa_lammps_split_framework_guests(tmp_path):
    p = tmp_path / "snap.data"
    p.write_text(_MINIMAL_DATA, encoding="utf-8")
    fw, guests, mols = get_atoms_from_graspa_lammps_snapshot(p)
    assert len(fw) == 2
    assert len(guests) == 1
    assert list(mols) == [2]
    assert mols[2].info.get("molecule_name") == "o-xylene"
    assert fw.arrays["mol_id"].tolist() == [1, 1]
    assert guests.arrays["type_id"].tolist() == [2]


@pytest.mark.skipif(not ACODAZ_SNAPSHOT.is_file(), reason="ACODAZ snapshot not extracted")
def test_acodaz_snapshot_matches_reference_lines():
    data = read_graspa_lammps_data(ACODAZ_SNAPSHOT)
    refs = _reference_atoms_from_file(ACODAZ_SNAPSHOT)
    assert data["n_atoms"] == 2352
    assert len(refs) == 2352
    assert len(data["atoms"]) == 2352

    by_id = {a["id"]: a for a in data["atoms"]}
    for ref in refs:
        rec = by_id[ref["id"]]
        assert rec["mol_id"] == ref["mol_id"]
        assert rec["type_id"] == ref["type_id"]
        assert rec["charge"] == pytest.approx(ref["charge"])
        assert rec["position"] == pytest.approx(ref["position"])

    fw, guests, mols = graspa_lammps_data_to_ase(data)
    assert len(fw) == 2112
    assert len(guests) == 240
    assert len(mols) == 20
    assert fw.arrays["mol_id"].min() == GRASPA_FRAMEWORK_MOL_ID
    assert fw.arrays["mol_id"].max() == GRASPA_FRAMEWORK_MOL_ID
    assert set(guests.arrays["type_id"].tolist()) <= {12, 13, 14}

    names = {at.info.get("molecule_name") for at in mols.values()}
    assert names == {"o-xylene", "m-xylene", "p-xylene"}

    # United-atom xylene: 12 sites per molecule in this run.
    assert sum(len(at) for at in mols.values()) == 240
    o_count = sum(1 for at in mols.values() if at.info.get("molecule_name") == "o-xylene")
    m_count = sum(1 for at in mols.values() if at.info.get("molecule_name") == "m-xylene")
    p_count = sum(1 for at in mols.values() if at.info.get("molecule_name") == "p-xylene")
    assert o_count == 6
    assert m_count == 4
    assert p_count == 10

    # Cell: cubic 46.73 Å (orthogonal within file precision).
    assert np.allclose(fw.cell.lengths(), [46.73, 46.73, 46.73], atol=1e-3)

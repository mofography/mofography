"""Tests for PACMOF model path resolution."""

from pathlib import Path

import pytest

from mofography.charges import DEFAULT_PACMOF_REPO_MODEL, resolve_pacmof_model_path


def test_resolve_pacmof_model_explicit(tmp_path):
    p = tmp_path / "model.pkl"
    p.write_bytes(b"x")
    assert resolve_pacmof_model_path(p) == p.resolve()


def test_resolve_pacmof_model_missing_explicit(tmp_path):
    with pytest.raises(FileNotFoundError):
        resolve_pacmof_model_path(tmp_path / "nope.pkl")


@pytest.mark.skipif(
    not DEFAULT_PACMOF_REPO_MODEL.is_file(),
    reason="repo pacmof model not present",
)
def test_resolve_pacmof_repo_default():
    p = resolve_pacmof_model_path()
    assert p.is_file()
    assert p.name == "Model_RF_DDEC.pkl"

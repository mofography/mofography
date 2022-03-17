"""
Unit and regression test for the mofography package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import mofography


def test_mofography_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "mofography" in sys.modules

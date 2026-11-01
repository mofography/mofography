from mofography.pocket_typing import DEFAULT_SPACING_A, _segmentation_spacing


def test_spacing_read_from_input_params():
    assert _segmentation_spacing({"input_params": {"spacing": 0.4}}) == 0.4


def test_spacing_top_level_wins_and_default():
    assert _segmentation_spacing({"spacing": 0.3, "input_params": {"spacing": 0.4}}) == 0.3
    assert _segmentation_spacing({}) == DEFAULT_SPACING_A

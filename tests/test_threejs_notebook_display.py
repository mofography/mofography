import html
from pathlib import Path

import pytest

ase_io = pytest.importorskip("ase.io")
from mofography import Figure3js  # noqa: E402

CIF = Path(__file__).resolve().parents[1] / "examples" / "TestCIFs" / "P1-IRMOF-1.cif"


def _fig():
    atoms = ase_io.read(CIF)
    return Figure3js().plot_framework(atoms).plot_cell(atoms), atoms


def test_repr_html_iframe_srcdoc():
    fig, _ = _fig()
    out = fig._repr_html_()
    assert out.startswith("<iframe") and 'srcdoc="' in out
    inner = html.unescape(out.split('srcdoc="', 1)[1].rsplit('"', 1)[0])
    assert "<canvas" in inner or "THREE" in inner
    assert 'id="resetCam"' in inner and "cell" in inner
    assert fig._repr_html_() != out  # unique id per figure
    assert "text/html" in fig._repr_mimebundle_()


def test_save_html_and_height(tmp_path):
    fig, _ = _fig()
    assert fig.save_html(tmp_path / "f.html").exists()
    assert "height:420px" in fig.figure._iframe_html(height=420)


def test_nbclient_display_data():
    nbformat = pytest.importorskip("nbformat")
    nbclient = pytest.importorskip("nbclient")
    nb = nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell(
        "import ase.io\nfrom mofography import Figure3js\n"
        f"a = ase.io.read({str(CIF)!r})\n"
        "fig = Figure3js().plot_framework(a).plot_cell(a)\nfig.show(height=500)\nfig")])
    nbclient.NotebookClient(nb, timeout=300).execute()
    outs = [o for o in nb.cells[0].outputs if "data" in o]
    assert len(outs) == 2, nb.cells[0].outputs
    assert all("<iframe" in o["data"]["text/html"] for o in outs)

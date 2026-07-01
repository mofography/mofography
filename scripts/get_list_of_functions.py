#!/usr/bin/env python3
"""Regenerate the mofography API reference.

What it does
------------
1. Introspects the *live* ``mofography`` package to collect public callables
   from every first-level submodule that exposes at least one.
2. Writes a fresh ``docs/_quarto.yml`` ``quartodoc.sections`` block with one
   section per module (skipping empty modules).
3. Rewrites the corresponding ``docs/_sidebar.yml`` section list to match.
4. Invokes ``python -m quartodoc build`` to (re)create per-function ``.qmd``
   files under ``docs/reference/``.

Run from the repo root::

    python internal_docs/get_list_of_functions.py
"""

from __future__ import annotations

import importlib
import inspect
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Sequence

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
DOCS_DIR = REPO_ROOT / "docs"
QUARTO_YML = DOCS_DIR / "_quarto.yml"
SIDEBAR_YML = DOCS_DIR / "_sidebar.yml"
REFERENCE_DIR = DOCS_DIR / "reference"

# Modules to include in the API reference, in display order. Add new top-level
# mofography submodules here as they grow.
MODULES: Sequence[str] = (
    "compute",
    "visualize",
    "analyze",
    "data",
    "helper",
    "rust_fun",
    "workflow",
    "threejs_viz",
    "threejs_plotting",
    "viz_backend",
    "movies",
    "winding_percolation",
    "pocket_typing",
    "pocket_framework_env",
    "widest_path_resegment",
    "window_refinement",
    "export_mp4",
    "percolation_results",
    "vdw_radii",
    "cif_io",
    "raspa_helper",
    "poreblazer",
)

# Friendly one-line descriptions used for the section `desc:` field in
# _quarto.yml and the page-level blurb in _sidebar.yml.
MODULE_DESCRIPTIONS: dict[str, str] = {
    "compute": "Core algorithms for grids, regions, RAGs, paths, percolation, and surface area.",
    "visualize": "Plotly-based 3D figures and traces for frameworks, regions, and RAGs.",
    "analyze": "Cheminformatics, ring finding, functional groups, and graph conversions.",
    "data": "MOFdb fetchers, SQLite graph store, and CIF/PDB readers.",
    "helper": "Grid utilities, PBC helpers, PSD peak finding, and small math routines.",
    "rust_fun": "Python bindings to the Rust backend (distance/energy grids, segmentation, watershed).",
    "workflow": "High-level pipelines that wire compute + visualize for common analyses.",
    "threejs_viz": "Three.js trace helpers (Figure3JS / Figure3js) that mirror the Plotly API.",
    "threejs_plotting": "Standalone Three.js HTML builders for percolation scenes.",
    "viz_backend": "Backend dispatch utilities (Plotly vs Three.js) for figures and traces.",
    "movies": "Frame builders and movie-HTML helpers for time-varying scenes.",
    "winding_percolation": "Graph-winding-based accessibility and percolation routines.",
    "pocket_typing": "Pocket classification (chemistry/shape) on segmentation results.",
    "pocket_framework_env": "Framework fingerprints and probe-radius helpers for pocket typing.",
    "widest_path_resegment": "Pore-graph widest-path resegmenting for PBC wrap repairs.",
    "window_refinement": "Morphological window refinement on connections.",
    "export_mp4": "MP4 export helpers for Plotly/Three.js figures.",
    "percolation_results": "Result-pack accessors for percolation history and pore diameter.",
    "vdw_radii": "Van der Waals radii tables (Bondi, Mendeleev, typing) and metadata.",
    "cif_io": "CIF parsing with large-file support.",
    "raspa_helper": "RASPA input/output helpers, blocking spheres, and accessibility masks.",
    "poreblazer": "External Poreblazer runner (no public functions, see class `PoreblazerRunner`).",
}

# Short titles used in the sidebar / quartodoc sections.
MODULE_TITLES: dict[str, str] = {
    "compute": "Compute",
    "visualize": "Visualize",
    "analyze": "Analyze",
    "data": "Data",
    "helper": "Helper",
    "rust_fun": "Rust backend",
    "workflow": "Workflows",
    "threejs_viz": "Three.js viz",
    "threejs_plotting": "Three.js plotting",
    "viz_backend": "Viz backend",
    "movies": "Movies",
    "winding_percolation": "Winding percolation",
    "pocket_typing": "Pocket typing",
    "pocket_framework_env": "Pocket framework env",
    "widest_path_resegment": "Widest-path resegment",
    "window_refinement": "Window refinement",
    "export_mp4": "Export MP4",
    "percolation_results": "Percolation results",
    "vdw_radii": "VDW radii",
    "cif_io": "CIF I/O",
    "raspa_helper": "RASPA helper",
    "poreblazer": "Poreblazer",
}


def _public_names(module) -> List[str]:
    """Return the module's public function names, preferring ``__all__``."""
    all_attr = getattr(module, "__all__", None)
    if all_attr:
        names = list(all_attr)
    else:
        names = [
            n
            for n, obj in inspect.getmembers(module, inspect.isfunction)
            if not n.startswith("_") and obj.__module__ == module.__name__
        ]
    # Stable order, drop leading underscores if any slipped through.
    return sorted({n for n in names if not n.startswith("_")})


def collect_module_functions(
    modules: Iterable[str], package_root: Path
) -> dict[str, List[str]]:
    """Return ``{module_name: [function_name, ...]}`` for every importable module."""
    sys.path.insert(0, str(package_root))
    collected: dict[str, List[str]] = {}
    for name in modules:
        full = f"mofography.{name}"
        try:
            mod = importlib.import_module(full)
        except Exception as exc:  # noqa: BLE001
            print(f"[skip] {full}: import failed ({exc})")
            continue
        funcs = _public_names(mod)
        if not funcs:
            print(f"[skip] {full}: no public functions")
            continue
        collected[name] = funcs
        print(f"[ok]   {full}: {len(funcs)} functions")
    return collected


def build_quartodoc_sections(
    funcs_by_module: dict[str, List[str]]
) -> List[dict]:
    """Build the ``quartodoc.sections`` list for ``_quarto.yml``."""
    sections: List[dict] = []
    for name, funcs in funcs_by_module.items():
        sections.append(
            {
                "title": MODULE_TITLES.get(name, name.replace("_", " ").title()),
                "package": f"mofography.{name}",
                "desc": MODULE_DESCRIPTIONS.get(name, "Description of this module."),
                "contents": funcs,
            }
        )
    return sections


def update_quarto_yml(
    funcs_by_module: dict[str, List[str]], path: Path = QUARTO_YML
) -> None:
    """Replace the ``quartodoc.sections`` block in ``_quarto.yml``."""
    with path.open("r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    cfg.setdefault("quartodoc", {})["sections"] = build_quartodoc_sections(
        funcs_by_module
    )
    with path.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False, default_flow_style=False)
    print(f"[write] {path}")


def build_sidebar_contents(
    funcs_by_module: dict[str, List[str]]
) -> List[dict]:
    """Build the inner ``- contents:`` entries of the reference sidebar block."""
    out: List[dict] = []
    for name, funcs in funcs_by_module.items():
        out.append(
            {
                "contents": [f"reference/{f}.qmd" for f in funcs],
                "section": MODULE_TITLES.get(name, name.replace("_", " ").title()),
            }
        )
    return out


def update_sidebar_yml(
    funcs_by_module: dict[str, List[str]], path: Path = SIDEBAR_YML
) -> None:
    raw = path.read_text(encoding="utf-8")
    # The reference block starts with the first top-level "- contents:" that
    # contains the reference/index.qmd entry, and ends at the dummy-sidebar
    # placeholder that follows it.  Because YAML round-tripping is lossy on
    # comments, do a surgical string replacement on the reference block.
    start_marker = "  - contents:\n    - reference/index.qmd\n"
    end_marker = "  - id: dummy-sidebar"
    if start_marker not in raw or end_marker not in raw:
        raise RuntimeError(
            f"Could not locate sidebar reference block markers in {path}"
        )
    start = raw.index(start_marker)
    end = raw.index(end_marker)
    new_block = (
        "  - contents:\n"
        "    - reference/index.qmd\n"
    )
    for entry in build_sidebar_contents(funcs_by_module):
        new_block += "    - contents:\n"
        for item in entry["contents"]:
            new_block += f"      - {item}\n"
        new_block += f"      section: {entry['section']}\n"
    new_block += "    id: reference\n"
    new_raw = raw[:start] + new_block + raw[end:]
    path.write_text(new_raw, encoding="utf-8")
    print(f"[write] {path}")


def regenerate_reference_qmds(docs_dir: Path = DOCS_DIR) -> None:
    """Run ``python -m quartodoc build`` to refresh per-function qmd files."""
    cmd = [sys.executable, "-m", "quartodoc", "build"]
    print(f"[exec]  {' '.join(cmd)}  (cwd={docs_dir})")
    result = subprocess.run(
        cmd, cwd=str(docs_dir), capture_output=True, text=True
    )
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise SystemExit(
            f"quartodoc build failed with exit code {result.returncode}"
        )
    # Print last 40 lines of output so the user sees progress.
    tail = "\n".join(result.stdout.splitlines()[-40:])
    print(tail)


def main() -> None:
    funcs_by_module = collect_module_functions(MODULES, REPO_ROOT)
    if not funcs_by_module:
        raise SystemExit("No modules produced any public functions; aborting.")
    update_quarto_yml(funcs_by_module)
    update_sidebar_yml(funcs_by_module)
    regenerate_reference_qmds()
    total = sum(len(v) for v in funcs_by_module.values())
    print(
        f"\n[done] {len(funcs_by_module)} modules, {total} functions written."
    )


if __name__ == "__main__":
    main()

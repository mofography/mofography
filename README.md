<div align="center">

<img src="./docs/mofography_logo_dark_background.png#gh-dark-mode-only" alt="mofography" width="420" />
<img src="./docs/mofography_logo_white_background.png#gh-light-mode-only" alt="mofography" width="420" />

<img src="docs/feature_strip.svg" alt="Segmentation · Textural Properties · Pore Windows · Percolated Clusters · Pore Graphs · 3D Visualization · High-Throughput · ML Featurization" width="880" />

*Characterization toolkit for MOFs, zeolites, polymers, and other crystalline nanoporous materials*

[![License](https://img.shields.io/badge/License-BSD--3--Clause-blue?logo=opensourceinitiative&logoColor=white)](LICENSE)
[![Release](https://img.shields.io/badge/Release-2026.10.09.1-success?logo=github&logoColor=white)](docs/releases-and-citing.md)
[![Python](https://img.shields.io/badge/Python-3.10%E2%80%933.12-3776ab?logo=python&logoColor=white)](#install)
[![Rust](https://img.shields.io/badge/Backend-Rust-dea584?logo=rust&logoColor=black)](docs/rust_backend_documentation.md)
[![Cite this code](https://img.shields.io/badge/Cite-this_code-2563eb?logo=zenodo&logoColor=white)](docs/releases-and-citing.md)
[![Web demo](https://img.shields.io/badge/Try-Web_demo-ff6b35?logo=huggingface&logoColor=white)](https://arung-northwestern-mofography-web-demo.hf.space/)
[![Docs](https://img.shields.io/badge/Docs-Quarto_(preview)-0ea5e9?logo=quarto&logoColor=white)](https://mofography.github.io/mofography/)

[Features](#features) ·
[Who is it for?](#who) ·
[Install](#install) ·
[Quick test](#quick-test) ·
[Tutorials](#tutorials) ·
[Cite](#cite)

<img src="docs/mofography_overview.gif" alt="mofography workflow: framework, pore segmentation, pore graph, unfolded views" width="100%" />

</div>

**mofography** segments the void space of crystalline nanoporous materials into pockets and builds **pore graphs** from them. One segmentation gives you global and pocket-level textural properties, pore types, windows, and percolation, plus interactive 3D views and graph descriptors for screening and machine learning (see [Features](#features)).

Core numerics run in a **Rust** backend (`mgr_rust`) behind a **Python** API.

> **Canonical code:** [mofography/mofography](https://github.com/mofography/mofography) · **Latest release:** [`2026.10.09.1`](https://github.com/mofography/mofography/releases/tag/2026.10.09.1) · branch `master`

---

<a id="features"></a>

## Features

Everything below is available from the Python module. Global properties and pocket-level properties come from the same segmentation, so you can go from a material-level number down to the individual pore that produced it.

| Capability | What you get | Learn more |
|---|---|---|
| **Pore diameters** | Largest cavity diameter (LCD), pore limiting diameter (PLD), and largest included sphere (LIS) along percolating pathways | [Compute pore diameters](docs/compute_pore_diameters.ipynb), [PLD and percolation](docs/pld_and_percolation.ipynb) |
| **Surface area and void volume** | Total (probe-center occupiable) and accessible (probe-center accessible) values, globally and per pocket | [Surface area](docs/surface_area.qmd), [Accessibility](docs/accessibility.qmd) |
| **Pore size distributions** | Global distributions, and pocket-wise distributions that do not undersample small pockets | [Pore size distributions](docs/pore_size_distributions.qmd) |
| **Pocket segmentation** | The void space split into convex pockets with periodic boundary conditions handled | [Region segmentation](docs/region_segmentation.qmd), [PBC considerations](docs/pbc_considerations.qmd) |
| **Pocket types** | Fingerprinting of pockets by their local environment, to group chemically and geometrically similar pores | [Pocket type detection](docs/pockets_detection.qmd) |
| **Pore windows** | Detection of the windows that connect pockets, including across periodic boundaries | [Connections and windows](docs/connections_and_windows.qmd), [Window fingerprinting](docs/window_fingerprinting.qmd) |
| **Pore graphs** | Pockets as nodes and windows as edges, annotated with local properties, as input for analysis and graph learning | [Region adjacency graph](docs/region_adjacency_graph.ipynb), [Graph neural network model on pore graphs (PoroNet)](docs/PoroNet.qmd) |
| **Percolation** | Percolation thresholds for a given probe, and 3D plots of the accessible, percolating clusters of pockets | [PLD and percolation](docs/pld_and_percolation.ipynb) |
| **Visualization** | Interactive 3D views of the framework, pore surfaces, pockets, and pore graphs | [Visualization](docs/framework_visualization_example.ipynb) |
| **High throughput** | Parallel characterization of whole databases with a Rust core | [High-throughput characterization](docs/hts_using_mofography.ipynb) |

Not a Python user? The same core is available as a **CLI** ([`mgr_cli_py`](https://github.com/arung-northwestern/mgr_cli_py)) for terminal users and AI agents, and as a **web app** for anyone who prefers a GUI and does not code or use a terminal. Both expose a subset of the features above; see [Who should use what](#who).

Using an AI coding agent? Copy [`skills/mofography/SKILL.md`](skills/mofography/SKILL.md) into your agent's skills folder (for Claude Code, `.claude/skills/mofography/`) so it knows the Python API, the textural-property and percolation functions, and common pitfalls.

---

<a id="who"></a>

## Who should use what?

Pick your path by **role**. Three tools, one shared science core:

- **① Library** — this repo (`mofography` org)
- **② CLI + Streamlit** — [`mgr_cli_py`](https://github.com/arung-northwestern/mgr_cli_py) (personal repo)
- **③ Web demo** — [Hugging Face Space](https://arung-northwestern-mofography-web-demo.hf.space/) (hosted)

```mermaid
%%{init: {'themeVariables': {'fontSize': '20px', 'fontFamily': 'Helvetica, Arial, sans-serif'}, 'flowchart': {'htmlLabels': true, 'curve': 'basis', 'nodeSpacing': 28, 'rankSpacing': 40, 'padding': 10}}}%%
flowchart TB
  V("<b>Just curious?</b>")
  E("<b>A few MOFs?</b>")
  G("<b>Scripts / agents?</b>")
  R("<b>Research / ML?</b>")

  HF("<b>Web demo</b> ③<br/>No install · hosted")
  CLI("<b>CLI + app</b> ②<br/>Local · mgr_cli_py")
  LIB("<b>Library</b> ①<br/>Fastest · this repo")

  V --> HF
  E --> CLI
  G --> CLI
  R --> LIB

  classDef ask fill:transparent,stroke:#94a3b8,stroke-width:1.5px,stroke-dasharray:4 3
  classDef hf fill:#f9731614,stroke:#f97316,stroke-width:2.5px
  classDef cli fill:#14b8a614,stroke:#14b8a6,stroke-width:2.5px
  classDef lib fill:#8b5cf614,stroke:#8b5cf6,stroke-width:2.5px
  class V,E,G,R ask
  class HF hf
  class CLI cli
  class LIB lib
  linkStyle 0 stroke:#f97316,stroke-width:2px
  linkStyle 1,2 stroke:#14b8a6,stroke-width:2px
  linkStyle 3 stroke:#8b5cf6,stroke-width:2px
```

| | **① Library** (this repo) | **② CLI + web app** (`mgr_cli_py`) | **③ Web demo** (HF) |
|---|---|---|---|
| **Best for** | Papers, notebooks, HTS, custom analysis | A few structures, interactive UI, automation without a full env | Try before you install |
| **Install** | Clone + build `mgr_rust` ([below](#install)) | Separate repo; `pip install -e .` + Rust ([README](https://github.com/arung-northwestern/mgr_cli_py)) | None — [open demo](https://arung-northwestern-mofography-web-demo.hf.space/) |
| **Performance** | **Fastest** on your hardware (`--release`, tuned to your CPU) | Fast locally; extra caching vs library | Shared **queue**, **~2 vCPUs** — often slower than local; **16 GB RAM** helps large cells |
| **Maintainer** | [`mofography`](https://github.com/mofography/mofography) org | Personal repo [`arung-northwestern/mgr_cli_py`](https://github.com/arung-northwestern/mgr_cli_py) (will go public; PRs via maintainer) | Hosted demo only |
| **AI / agents** | `import mofography` in code | **`mgrcli`** subprocess — no notebook context | Not for batch/agent pipelines |

> **Tip:** For serious work on medium/large cells, prefer **①** or **②** on your machine. Use **③** to explore the UI and export HTML; expect queueing at busy times.

**Self-hosted web app:** clone [`mgr_cli_py`](https://github.com/arung-northwestern/mgr_cli_py), install, then run Streamlit interactive mode (see that repo's README). It is **not** the same package as `pip install mofography` from this tree.

---

<a id="cite"></a>

## Cite this code

- **Software:** release tag **`2026.10.09.1`** + commit hash — [docs/releases-and-citing.md](docs/releases-and-citing.md), [CITATION.cff](CITATION.cff).
- **Paper:**
  - (Software) A. Gopalan*, C. Zheng, B. Azdapour, K. Shi*, L. Sarkisov, "Advanced High-Throughput Characterization of Nanoporous Materials using Pore Graphs". arXiv *link when posted* (`PLACEHOLDER_ARXIV_ID` — [docs/PAPER_METADATA_TODO.md](docs/PAPER_METADATA_TODO.md)).
  - (Pore graph) C. Zheng, A. Gopalan, K. Shi*, "PoroNet: An Intrinsically Interpretable Pore Graph Neural Network for Resolving Pore-Level Adsorption in Metal–Organic Frameworks". *J. Chem. Theory Comput.* 2026, 22 (11), 5716–5730. [link](https://doi.org/10.1021/acs.jctc.6c00100)

```python
import mofography as m
print(m.__version__, m.__git_revision__)
```

---

<a id="install"></a>

## Install the library

This repository is the **full Python package**. There is no PyPI wheel yet: you install dependencies, compile **`mgr_rust`** with [maturin](https://www.maturin.rs/) in **release** mode (important for speed on your CPU), then `pip install -e .`.

**Python:** 3.10–3.12 (3.11 recommended) · **Rust:** [rustup](https://rustup.rs/) · **OS:** Linux, macOS, Windows

**Pin a paper/release checkout (optional):**

```bash
git clone https://github.com/mofography/mofography.git && cd mofography
git checkout 2026.10.09.1
```

### One command (Linux / macOS)

**A. Fresh sandbox with [uv](https://github.com/astral-sh/uv)** (creates `.venv` in the clone):

```bash
curl -fsSL https://raw.githubusercontent.com/mofography/mofography/master/scripts/install_mofography.sh | bash -s -- --mode uv-sandbox
```

**B. Active conda / venv** (installs into the Python you already activated):

```bash
curl -fsSL https://raw.githubusercontent.com/mofography/mofography/master/scripts/install_mofography.sh | bash -s -- --mode pip-existing
```

Pin a release in the one-liner:

```bash
curl -fsSL https://raw.githubusercontent.com/mofography/mofography/master/scripts/install_mofography.sh | bash -s -- --mode uv-sandbox --ref 2026.10.09.1
```

From an existing clone: `./scripts/install_mofography.sh --mode uv-sandbox` or `--mode pip-existing --repo-dir .`

### Windows (PowerShell)

```powershell
git clone https://github.com/mofography/mofography.git
cd mofography
git checkout 2026.10.09.1 # optional
.\scripts\install_mofography.ps1
```

### Manual steps (any OS)

```bash
conda create -n mgr_env python=3.11 -y && conda activate mgr_env # or your venv
conda install -c conda-forge pymatgen ase dask networkx scikit-image scikit-learn plotly freud rdkit -y

git clone https://github.com/mofography/mofography.git && cd mofography
python -m pip install maturin
cd mgr_rust && maturin develop --release && cd ..
python -m pip install -e .
```

**Verify:** `python -c "import mgr_rust, mofography as m; print(m.__version__)"`

After `git pull`, rebuild Rust: `cd mgr_rust && maturin develop --release`

More detail: [docs/rust_backend_documentation.md](docs/rust_backend_documentation.md), [docs/install.qmd](docs/install.qmd).

We do **not** ship prebuilt `mgr_rust` wheels per OS — local `--release` builds match your machine best.

---

<a id="quick-test"></a>

## Quick test

```python
import mofography as mgr
from ase.io import read

atoms = read("examples/TestCIFs/P1-IRMOF-1.cif")
figures = mgr.quick_cluster(atoms)
```

Interactive **Three.js HTML** is not stored in this repo; export from the **web demo** or [`mgr_cli_py`](https://github.com/arung-northwestern/mgr_cli_py), or regenerate with `examples/gen_*_threejs.py`.

---

<a id="tutorials"></a>

## Tutorials & documentation

Tutorials are first-class: Jupyter notebooks and Quarto pages live under [`docs/`](docs/). The site is built with **Quarto** (+ **quartodoc** API reference) into `docs/_site/`:

```bash
cd docs && quarto render
open _site/index.html # macOS; or use a local static server
```

**Planned public URL (GitHub Pages when `mofography/mofography` is public):** [https://mofography.github.io/mofography/](https://mofography.github.io/mofography/) — the Docs badge above will point there after deploy. Until then, render locally or browse notebook sources on GitHub under `docs/`.

### Basic tutorials (recommended order)

1. [Quick start](docs/quick_start.qmd) — end-to-end workflow **(new)**
2. [Basic structure analysis](docs/basic_structure_io.ipynb)
3. [Basic structure visualization](docs/quick_structure_visualization.ipynb)
4. [Distance grid](docs/distance_grid.qmd)
5. [Region segmentation](docs/region_segmentation.qmd)
6. [PBC considerations](docs/pbc_considerations.qmd)
7. [Region adjacency graph (RAG)](docs/region_adjacency_graph.ipynb)
8. [Connections and windows](docs/connections_and_windows.qmd)
9. [Pocket type detection](docs/pockets_detection.qmd) (environment-cluster method) **(new)**

### Advanced tutorials

- [Three.js plotting](docs/threejs_plotting.qmd) **(new)**
- [Watershed backends](docs/watershed_backends.qmd) **(new)**
- [Fast Rust CIF reader](docs/fast_cif_reader.qmd) **(new)**
- [Advanced structure visualization](docs/framework_visualization_example.ipynb)
- [Framework subset visualization](docs/plot_subset_as_mesh_tutorial.ipynb)
- [PLD and percolation history](docs/pld_and_percolation.ipynb)
- [Compute pore diameters](docs/compute_pore_diameters.ipynb)
- [MOFDB example](docs/mofdb_example.ipynb)
- [High-throughput characterization](docs/hts_using_mofography.ipynb)
- [Periodic segmentation Rust vs Python (HTS)](docs/hts_periodic_seg_benchmark.ipynb)
- [Functional group detection](docs/detect_functional_groups_in_dataset.ipynb)
- [Graph neural network model on pore graphs (PoroNet)](docs/PoroNet.qmd) **(new)**

Quarto guides (not only `.ipynb`): [install](docs/install.qmd), [CLI quick start](docs/cli_quickstart.qmd), [learn hub](docs/learn.qmd), [API reference](docs/reference/index.qmd).

---

## Releases & metadata

Calendar-version tags and collaborator citations: [docs/releases-and-citing.md](docs/releases-and-citing.md). Maintainer TODOs (arXiv, Zenodo, `mgr_cli_py` links): [docs/PAPER_METADATA_TODO.md](docs/PAPER_METADATA_TODO.md).

---

## Acknowledgements

This material is based upon work supported by the ACS Petroleum Research Fund under Doctoral New Investigator Grant 68604-DNI6 and NSF CAREER Award No. 2543449.

---

## Copyright

Copyright © 2022–2026, Arun Gopalan and contributors. [BSD-3-Clause](LICENSE).

Scaffolding from the [MolSSI Python Cookiecutter](https://github.com/molssi/cookiecutter-cms).

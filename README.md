mofography
==========

[![License](https://img.shields.io/github/license/mofography/mofography)](LICENSE)
[![Release](https://img.shields.io/github/v/release/mofography/mofography?label=release)](https://github.com/mofography/mofography/releases)
[![Cite this code](https://img.shields.io/badge/Cite-this_code-2563eb)](docs/releases-and-citing.md)
[![Web demo](https://img.shields.io/badge/Try-Web_demo-ff6b35)](https://arung-northwestern-mofography-web-demo.hf.space/)

![](./docs/mofography_logo_dark_background.png#gh-dark-mode-only)
![](./docs/mofography_logo_white_background.png#gh-light-mode-only)

**mofography** segments void space into pockets, builds **pore graphs**, and supports visualization and featurization for crystalline nanoporous materials (MOFs, zeolites, polymers, and related frameworks). Core numerics use a **Rust** backend (`mgr_rust`) with a **Python** API.

<p align="center">
  <img src="docs/mofography_nu1000_exploded.png" alt="NU-1000 pore segmentation (exploded view)" width="100%" />
</p>

**Canonical code:** [github.com/mofography/mofography](https://github.com/mofography/mofography) · **Latest release:** [2026.10.09.1](https://github.com/mofography/mofography/releases/tag/2026.10.09.1) · branch `master`

---

## Three ways to use mofography

The **library** (this repo), the **CLI + self-hosted web app** ([`mgr_cli_py`](https://github.com/arung-northwestern/mgr_cli_py)), and the **public Hugging Face demo** share similar science but differ in depth, performance, and who maintains them.

```mermaid
%%{init: {'theme': 'base', 'themeVariables': {'fontSize': '14px', 'fontFamily': 'ui-sans-serif, system-ui', 'primaryBorderColor': '#6366f1', 'lineColor': '#64748b'}}}%%
flowchart LR
  START(["🎯 Start here"])

  START --> Q1{"📓 Notebooks,\nML, HTS?"}
  START --> Q2{"🖥️ Browser or\nterminal locally?"}
  START --> Q3{"☁️ Try first,\nno install?"}

  Q1 -->|research path| LIB["🧪 ① Python library\nmofography/mofography"]
  Q2 -->|app path| CLI["⚡ ② mgr_cli_py\nCLI + Streamlit"]
  Q3 -->|demo path| HF["🌐 ③ HF Space\ndemo queue"]

  LIB --> LIBd["🔬 Full API · build mgr_rust\n--release on your CPU"]
  CLI --> CLId["🚀 mgrcli · caches · self-host\narung-northwestern/mgr_cli_py"]
  HF --> HFd["⏳ ~2 vCPU · 16 GB RAM\nexport HTML · queue waits"]

  LIBd -.->|agents OK| AG1["🤖 import mofography"]
  CLId -.->|agents OK| AG2["🤖 mgrcli subprocess"]
  HFd -.->|not for batch| AG3["❌ agents / HTS"]

  classDef start fill:#eef2ff,stroke:#4f46e5,color:#312e81,stroke-width:2px
  classDef question fill:#faf5ff,stroke:#a855f7,color:#581c87,stroke-width:1px
  classDef lib fill:#4f46e5,stroke:#312e81,color:#ffffff,stroke-width:2px
  classDef cli fill:#059669,stroke:#047857,color:#ffffff,stroke-width:2px
  classDef hf fill:#ea580c,stroke:#c2410c,color:#ffffff,stroke-width:2px
  classDef detail fill:#f1f5f9,stroke:#94a3b8,color:#0f172a,stroke-width:1px
  classDef agent fill:#fef3c7,stroke:#d97706,color:#78350f,stroke-width:1px
  classDef agentbad fill:#fee2e2,stroke:#dc2626,color:#7f1d1d,stroke-width:1px

  class START start
  class Q1,Q2,Q3 question
  class LIB lib
  class CLI cli
  class HF hf
  class LIBd,CLId,HFd detail
  class AG1,AG2 agent
  class AG3 agentbad
```

| | **① Library (this repo)** | **② CLI + web app (`mgr_cli_py`)** | **③ Public web demo (HF)** |
|---|---------------------------|-------------------------------------|----------------------------|
| **Best for** | Papers, notebooks, HTS, custom analysis | A few structures, interactive UI, automation without full env | Try before you install |
| **Install** | Clone + build `mgr_rust` ([below](#install-the-library)) | Separate repo; `pip install -e .` + Rust ([mgr_cli_py README](https://github.com/arung-northwestern/mgr_cli_py)) | None — [open demo](https://arung-northwestern-mofography-web-demo.hf.space/) |
| **Performance** | **Fastest** on your hardware (`maturin develop --release` tuned to your CPU) | Fast locally; extra caching vs library | Shared **queue**, **~2 vCPUs** — often slower than local; **16 GB RAM** helps large cells |
| **Maintainer** | [mofography](https://github.com/mofography/mofography) org | Personal repo — [**arung-northwestern/mgr_cli_py**](https://github.com/arung-northwestern/mgr_cli_py) (will go public; PRs via maintainer) | Hosted demo only |
| **AI / agents** | Import `mofography` in code | **`mgrcli`** subprocess — no notebook context | Not for batch/agent pipelines |

> **Tip:** For serious work on medium/large cells, prefer **①** or **②** on your machine. Use **③** to explore the UI and export HTML; expect waiting in queue at busy times.

**Self-hosted web app:** clone [`mgr_cli_py`](https://github.com/arung-northwestern/mgr_cli_py), install, then run Streamlit interactive mode (see that repo’s README). It is **not** the same package as `pip install mofography` from this tree.

---

## Cite this code

- **Software:** release tag **`2026.10.09.1`** + commit hash — [docs/releases-and-citing.md](docs/releases-and-citing.md), [CITATION.cff](CITATION.cff).
- **Paper:** *Advanced High-Throughput Characterization of Nanoporous Materials using Pore Graphs* — Arun Gopalan, Chao Zheng, Behnam Azdapour, Kaihang Shi, Lev Sarkisov.  
  **arXiv:** *link when posted* (`PLACEHOLDER_ARXIV_ID` — [docs/PAPER_METADATA_TODO.md](docs/PAPER_METADATA_TODO.md)).

```python
import mofography as m
print(m.__version__, m.__git_revision__)
```

---

## Install the library

This repository is the **full Python package**. There is no PyPI wheel yet: you install dependencies, compile **`mgr_rust`** with [maturin](https://www.maturin.rs/) in **release** mode (important for speed on your CPU), then `pip install -e .`.

**Python:** 3.10–3.12 (3.11 recommended) · **Rust:** [rustup](https://rustup.rs/) · **OS:** Linux, macOS, Windows (see PowerShell script)

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
git checkout 2026.10.09.1   # optional
.\scripts\install_mofography.ps1
```

### Manual steps (any OS)

```bash
conda create -n mgr_env python=3.11 -y && conda activate mgr_env   # or your venv
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

## Quick test

```python
import mofography as mgr
from ase.io import read

atoms = read("examples/TestCIFs/P1-IRMOF-1.cif")
figures = mgr.quick_cluster(atoms)
```

Interactive **Three.js HTML** is not stored in this repo; export from the **web demo** or [`mgr_cli_py`](https://github.com/arung-northwestern/mgr_cli_py), or regenerate with `examples/gen_*_threejs.py`.

---

## Tutorials

Notebooks under [`docs/`](docs/): start with [quick cluster](docs/quick_cluster.ipynb) → [segmentation](docs/region_segmentation.ipynb) → [RAG](docs/region_adjacency_graph.ipynb) → [pocket types](docs/pockets_detection.ipynb). See `docs/` for HTS, percolation, and GCN examples.

---

## Releases & metadata

Calendar-version tags and collaborator citations: [docs/releases-and-citing.md](docs/releases-and-citing.md). Maintainer TODOs (arXiv, Zenodo, `mgr_cli_py` links): [docs/PAPER_METADATA_TODO.md](docs/PAPER_METADATA_TODO.md).

---

## Copyright

Copyright © 2022–2026, Arun Gopalan and contributors. [BSD-3-Clause](LICENSE).

Scaffolding from the [MolSSI Python Cookiecutter](https://github.com/molssi/cookiecutter-cms).

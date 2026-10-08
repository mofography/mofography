# Install mofography + mgr_rust (release) on Windows (PowerShell).
# Prerequisites: git, Rust (rustup), Python 3.10-3.12 on PATH.
#
#   git clone https://github.com/mofography/mofography.git
#   cd mofography
#   .\scripts\install_mofography.ps1

param(
    [ValidateSet("pip-existing", "venv")]
    [string]$Mode = "venv",
    [string]$RepoDir = (Get-Location).Path
)

$ErrorActionPreference = "Stop"
Set-Location $RepoDir

$deps = @(
    "pymatgen", "ase", "dask", "dask-jobqueue", "networkx", "scikit-image", "tqdm", "aiohttp",
    "sparse", "scikit-learn", "rich", "trimesh", "pyvista", "plotly", "furl", "requests",
    "freud-analysis", "rdkit", "optuna", "kaleido", "nbformat"
)

if ($Mode -eq "venv") {
    python -m venv .venv
    & .\.venv\Scripts\Activate.ps1
}

python -m pip install --upgrade pip maturin
python -m pip install $deps

Push-Location mgr_rust
maturin develop --release
Pop-Location

python -m pip install -e .

Write-Host "Verify: python -c `"import mgr_rust, mofography as m; print(m.__version__)`""

# Suggested Development Commands

## Environment Setup
```bash
# Create conda environment
conda create -n mgr_env python
conda activate mgr_env

# Install dependencies with conda
conda install -c conda-forge pymatgen ase dask dask-jobqueue networkx scikit-image tqdm aiohttp icecream sparse scikit-learn rich trimesh pyvista plotly furl requests freud nbformat cupy python-kaleido

# Install dependencies with pip (exclude cupy if no GPU)
python -m pip install pymatgen ase dask dask-jobqueue networkx scikit-image tqdm aiohttp icecream sparse scikit-learn rich trimesh pyvista plotly furl requests freud-analysis cupy kaleido
```

## Development Commands
```bash
# Install package in development mode
python setup.py install

# Add to Python path for development
export PYTHONPATH="${PYTHONPATH}:/path/to/mofography"

# Run tests
pytest
pytest --cov=mofography
pytest --cov=mofography --cov-report=html

# Code formatting
yapf -i *.py

# Linting
flake8 --max-line-length=119 *.py

# Build documentation (in docs directory)
cd docs
make html
make clean
```

## Testing Commands
```bash
# Create test environment
conda env create -f devtools/conda-envs/test_env.yaml
conda activate test

# Run specific test file
pytest mofography/tests/test_compute.py

# Run with coverage
pytest --cov=mofography --cov-report=term-missing

# Coverage report upload (if configured)
codecov
```

## Usage Examples
```bash
# Quick test of installation
python -c "import mofography as mgr; print('MOFography imported successfully')"

# Run basic clustering workflow
python -c "
import mofography as mgr
from ase.io import read
mat_atoms = read('examples/TestCIFs/P1-IRMOF-1.cif')
figures = mgr.quick_cluster(mat_atoms)
print('Analysis complete!')
"
```

## Git Commands
```bash
# Standard git workflow
git add .
git commit -m "Descriptive commit message"
git push origin master

# Pull latest changes
git pull origin master

# Check status
git status
git diff
```

## Dask/Parallel Computing (if using)
```bash
# Start dask scheduler for distributed computing
dask-scheduler

# Start workers
dask-worker scheduler-address:8786
```

## Rust Components (if modifying)
```bash
# Navigate to rust directory
cd mgr_rust

# Build rust components
cargo build
cargo build --release

# Run rust tests
cargo test
```
# Memory Block: Project Commands

## Build & Development Commands

### Rust Backend (mgr_rust/)
```bash
cd mgr_rust
cargo build --release  # Build optimized Rust library
cargo test               # Run Rust tests
cargo clippy             # Lint Rust code
```

### Python Package
```bash
# Development installation
pip install -e .

# Run tests (minimal currently)
pytest mofography/tests/test_mofography.py

# Build documentation
cd docs
quarto render
```

### Environment Setup
```bash
# From README - recommended dependencies
conda install -c conda-forge pymatgen ase dask dask-jobqueue networkx scikit-image tqdm aiohttp icecream sparse scikit-learn rich trimesh pyvista plotly furl requests freud nbformat cupy python-kaleido

# Alternative: pip
python -m pip install pymatgen ase dask dask-jobqueue networkx scikit-image tqdm aiohttp icecream sparse scikit-learn rich trimesh pyvista plotly furl requests freud-analysis cupy kaleido
```

### Documentation Generation
```bash
cd docs
quarto render  # Generate HTML site
quarto preview  # Local preview server
```

## Testing Commands
```bash
# Current basic test
pytest mofography/tests/ -v

# With coverage (if expanded)
pytest --cov=mofography mofography/tests/
```

## Key File Patterns
- **Distance Grid**: Use `dgrid_from_atoms_aabb_rust()` for new development (AABB acceleration)
- **Periodic Analysis**: Use `do_periodic_segmentation()` for crystalline materials
- **Old Workflow**: `quick_cluster()` still available but deprecated for periodic systems

## Performance Profiling
```python
import time
import mofography as mgr

# Compare old vs new methods
start = time.time()
dgrid_old = mgr.dgrid_from_atoms_cpu(atoms, spacing)
old_time = time.time() - start

start = time.time() 
dgrid_new = mgr.dgrid_from_atoms_aabb_rust(atoms, spacing)
new_time = time.time() - start

speedup = old_time / new_time
print(f"Speedup: {speedup:.2f}x")
```

## Last Updated: 2025-12-19
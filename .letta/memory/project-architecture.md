# Memory Block: Project Architecture

## Module Structure & Relationships

### Core Python Modules (`mofography/`)
1. **workflow.py** - Main entry point with `quick_cluster()` function
2. **compute.py** - All computational algorithms (segmentation, connections, RAGs)
3. **rust_fun.py** - Python-Rust bridge functions with comprehensive docstrings
4. **winding_percolation.py** - Advanced topologically correct percolation algorithms
5. **visualize.py** - 3D plotting and interactive visualization
6. **analyze.py** - Chemical analysis (functional groups, rings, bonds)
7. **data.py** - Data management and forcefield calculations

### Rust Performance Backend (`mgr_rust/`)
- **lib.rs** contains:
  - `dgrid_from_atoms()` - Basic distance field computation
  - `dgrid_from_atoms_aabb()` - **AABB accelerated version** (major performance win)
  - `dgrid_from_atoms_aabb_interp()` - Interpolated smoothing
  - `egrid_from_atoms()` - Lennard-Jones energy fields
  - Uses Rayon for parallel processing, R*-tree for neighbor search

### Documentation (`docs/`)
- **Quarto-based** site with `_quarto.yml` configuration
- **Auto-generated API reference** using quartodoc + sphinx parser
- **Tutorial notebooks** for learning path from basic to advanced

### Key Data Flow
```
ASE Atoms → Distance Grid (Rust AABB) → Periodic Segmentation → 
Regions → Connections → RAG + Percolation Analysis → Visualization
```

## Critical Design Patterns
1. **Performance First**: Computational bottlenecks moved to Rust
2. **PBC Awareness**: All algorithms handle periodic boundary conditions
3. **Modular Design**: Clear separation between compute, analyze, visualize layers
4. **Documentation Driven**: Extensive tutorials and examples

## Recent Development Focus
- AABB neighbor search implementation and optimization
- Periodic segmentation for crystalline accuracy
- Winding percolation for topological correctness

## Testing Strategy
- Minimal pytest setup currently (just import tests)
- Need expansion for comprehensive coverage

## Dependencies
- **Core**: ASE, NetworkX, NumPy, SciPy, scikit-image
- **Performance**: Dask (parallel), Rust/PyO3 (compute)
- **Visualization**: Plotly, PyVista, Trimesh
- **Optional**: CuPy (GPU acceleration)

## Last Updated: 2025-12-19
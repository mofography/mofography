# MOFography Project Overview

**Project Purpose**: MOFography is a powerful Python package for the analysis, visualization, featurization, and machine learning on molecular pores in nanoporous materials, particularly Metal-Organic Frameworks (MOFs).

**Tech Stack**:
- **Primary Language**: Python (3.10+ recommended)
- **Key Dependencies**: 
  - pymatgen, ase (materials science)
  - dask, dask-jobqueue (parallel computing)
  - networkx, scikit-image, scikit-learn (analysis/ML)
  - pyvista, trimesh, plotly (visualization)
  - freud-analysis (spatial analysis)
  - rich (terminal output)
  - cupy (GPU acceleration, optional)
- **Rust Integration**: mgr_rust directory contains Rust components for performance-critical operations
- **Documentation**: Quarto-based documentation system

**Code Style & Conventions**:
- **Line Length**: 119 characters (configured in setup.cfg)
- **Indentation**: 4 spaces, no tabs
- **Formatting**: YAPF (Yet Another Python Formatter)
- **Linting**: Flake8 with max-line-length of 119
- **Testing**: pytest with pytest-cov for coverage
- **Version Control**: Git with versioneer for automatic version numbering
- **Code Coverage**: Target 50% threshold via codecov

**Package Structure**:
```
mofography/
├── __init__.py          # Package initialization and imports
├── compute.py          # Core computational functions
├── analyze.py          # Analysis and featurization
├── visualize.py        # Visualization utilities
├── workflow.py         # Complete analysis workflows
├── helper.py           # Helper functions
├── data.py             # Data handling
├── rust_fun.py         # Rust function interfaces
├── raspa_helper.py     # RASPA simulation helper
├── tests/              # Test suite
├── examples/           # Example files
└── data/              # Package data
```

**Development Status**: Active development, not yet stable API
**Installation Methods**: 
1. Add to Python path (recommended for development)
2. Hard install via `python setup.py install`

**Key Features**:
- Pore size analysis and visualization
- Region segmentation and pore network analysis
- Periodic boundary condition handling
- High-throughput characterization with Dask
- Graph-based analysis (Region Adjacency Graphs)
- Machine learning featurization
- Integration with MOFDB and other materials databases
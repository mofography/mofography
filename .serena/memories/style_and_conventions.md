# MOFography Style and Conventions

## Code Style
- **Line Length**: 119 characters maximum
- **Indentation**: 4 spaces, no tabs
- **Import Style**: 
  - Standard library imports first
  - Third-party imports second
  - Local imports third
  - Use `from .module import *` for package modules in __init__.py
- **Docstrings**: Use standard Python docstring format
- **Type Hints**: Recommended but not strictly required

## Code Organization
- **Module Structure**: Each module focuses on a specific aspect (compute, analyze, visualize, etc.)
- **Function Naming**: Use descriptive names, preferably with snake_case
- **Class Naming**: Use CamelCase for classes
- **Constants**: Use UPPER_CASE for constants

## Development Practices
- **Testing**: Write tests for new functionality in the tests/ directory
- **Documentation**: Update documentation when adding new features
- **Version Control**: Use descriptive commit messages
- **Code Review**: Ensure code passes flake8 and yapf formatting before commits

## Performance Considerations
- **Rust Integration**: Performance-critical operations implemented in Rust (mgr_rust directory)
- **Parallel Computing**: Use Dask for high-throughput computations
- **GPU Support**: Optional CuPy support for GPU acceleration
- **Memory Efficiency**: Consider memory usage for large datasets

## File Naming Conventions
- **Python Files**: Use snake_case for module names
- **Test Files**: Prefix with test_ (e.g., test_compute.py)
- **Documentation**: Use .ipynb for tutorials, .md for documentation
- **Data Files**: Store in mofography/data/ directory

## Import Conventions
```python
# Standard library imports
import os
import sys
from typing import List, Optional

# Third-party imports
import numpy as np
import pandas as pd
from ase import Atoms

# Local imports
from .compute import compute_pore_sizes
from .analyze import analyze_structure
```

## Error Handling
- Use appropriate exception types
- Provide informative error messages
- Handle edge cases gracefully
- Use try-except blocks for operations that might fail

## Logging and Output
- Use `rich` for terminal output when available
- Use standard `logging` module for debugging
- Provide progress indicators for long-running operations
- Use `tqdm` for progress bars

## Configuration
- Use setup.cfg for package configuration
- Use .style.yapf for YAPF formatting rules
- Use .flake8 for linting configuration
- Environment-specific configs in devtools/ directory
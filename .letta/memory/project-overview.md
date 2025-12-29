# Memory Block: Project Overview

## High-Level Description
**MOFOGRAPHY** is a powerful Python package for analysis, visualization, featurization, and machine learning on molecular pores in nanoporous materials, with a focus on Metal-Organic Frameworks (MOFs) and zeolites.

## Key Innovation: Hybrid Python/Rust Architecture
- **Python Frontend**: User-friendly API, scientific workflows, visualization
- **Rust Backend**: High-performance computational kernels using PyO3 bindings
- **Performance Focus**: AABB neighbor search for distance fields, parallel processing

## Core Scientific Capabilities
1. **Distance Grid Computation**: CPU and GPU-accelerated with AABB neighbor search optimization
2. **Periodic Segmentation**: PBC-aware watershed segmentation for crystalline materials
3. **Region Adjacency Graphs (RAGs)**: Graph-based pore connectivity analysis
4. **Winding Percolation**: Topologically correct percolation using winding numbers
5. **Visualization**: 3D interactive plotting with Plotly/PyVista/Trimesh
6. **Machine Learning Integration**: Feature extraction for ML on pore structures

## Recent Algorithmic Breakthroughs
- **AABB Neighbor Search**: Major performance improvement for distance field computation
- **Periodic Segmentation**: Solves fundamental accuracy issues in crystalline pore analysis
- **Winding Number Percolation**: Avoids dead-end detection issues in amorphous pores

## Key Applications
- High-throughput materials screening
- Pore size distribution analysis
- Percolation pathway identification
- Functional group detection
- Molecular simulations integration

## Target Users
- Materials science researchers
- Computational chemists
- MOF/zeolite researchers
- Porous materials engineers

## Documentation Structure
- **Quarto-based** site with auto-generated API reference
- **Tutorial-driven** learning path from basic to advanced
- **Example-rich** with real MOF structures

## Last Updated: 2025-12-19
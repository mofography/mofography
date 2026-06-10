#!/usr/bin/env python3
"""
Generate publication-quality subfigures for mofography workflow flowchart
Following the standard data preparation workflow and plotting parameters
"""

from pathlib import Path
import logging
import numpy as np
from ase.io import read
from ase.build import make_supercell
from skimage.measure import regionprops
import plotly.express as px
import plotly.graph_objects as go
import mofography as mgr

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_subfigure(traces, filename, width=800, height=600, scale=2, orthographic=False):
    """Create publication-quality subfigure"""
    fig = go.Figure(data=traces, layout=mgr.white_layout)

    # Apply orthographic projection only if specifically requested
    if orthographic:
        fig = mgr.make_fig_orthographic(fig)

    # Export as high-resolution SVG
    fig.write_image(filename, width=width, height=height, scale=scale)
    logger.info(f"Created subfigure: {filename}")
    return fig

def main():
    """Generate three subfigures for the mofography workflow flowchart"""

    # Load and prepare material
    logger.info("Loading NU-1000 material structure...")
    mat_atoms = read('/Users/arungopalan/repos/mofography/mofography/examples/NU-1000_S1.cif')

    # Standard parameters from documentation
    spacing = 0.5
    mask_thickness = 1.5
    h = 0.2
    threshold_abs = 1.0
    min_distance = 1.5

    # Create appropriate supercell
    logger.info("Creating supercell...")
    nx_cells, ny_cells, nz_cells = mgr.get_repetitions_to_meet_cutoff(mat_atoms, 6)
    mat_atoms = make_supercell(mat_atoms, [[nx_cells, 0, 0], [0, ny_cells, 0], [0, 0, nz_cells]])

    # Generate distance grid and perform segmentation
    logger.info("Generating distance grid and performing segmentation...")
    dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=spacing)
    results = mgr.do_periodic_segmentation_rust(dgrid, mat_atoms,
                                               mask_thickness=mask_thickness, h=h,
                                               spacing=spacing, return_dgrid=True,
                                               threshold_abs=threshold_abs, min_distance=min_distance)

    # Extract key results
    dgrid = results['dgrid']
    regions = results['regions']
    regions_pbc = results['regions_pbc']
    maxima_info = results['maxima_info']
    pbc_groups = results['pbc_groups']
    connections = results['connections']
    all_maxima = results['all_maxima']
    unique_regions = results['unique_regions']

    # Prepare consistent maxima data
    maxima = all_maxima[:regions.max()]
    all_maxima_radii = dgrid[tuple(maxima_info[:,1:4].T)]
    maxima_radii = all_maxima_radii[:regions.max()]

    # CRITICAL: Use trimmed connections for consistency
    connections_trimmed = [c for c in connections if np.all(np.isin(c[0], unique_regions))]

    # Build RAG with proper sequence
    logger.info("Building region adjacency graph...")
    rag = mgr.rag_from_connections(connections_trimmed, maxima)
    rag = mgr.add_maxima_to_rag(rag, maxima, maxima_radii, dgrid.shape, mat_atoms)
    rag = mgr.add_pbc_to_rag(rag, pbc_groups)

    # Add analytical data
    props = regionprops(regions)
    bins = np.linspace(0, 50, 30)
    rag = mgr.add_psd_to_rag_with_pbc(rag, props, pbc_groups, bins, mat_atoms, dgrid.shape)
    rag = mgr.add_clustering_to_rag_geometry_js(rag, props, pbc_groups, mat_atoms, dgrid.shape)
    regions_cluster = mgr.clustering_grid_from_regions(regions, rag, cluster_key='cluster_geom_pbc')

    # Set consistent parameters
    isoval = 2.0
    cscale = px.colors.qualitative.Pastel  # Default color scheme

    logger.info("Generating traces for subfigures...")

    # Core traces (use connections_trimmed consistently)
    cell_trace = mgr.plot_cell(mat_atoms)
    frame_trace = mgr.plot_framework_as_mesh(mat_atoms, stylize=True)
    surface_trace_regions = mgr.plot_surface_dgrid(dgrid, isoval, mat_atoms, regions,
                                                 select='all', cscale=cscale, stylize=True)
    surface_trace_regions_pbc = mgr.plot_surface_dgrid(dgrid, isoval, mat_atoms, regions_pbc,
                                                     select='all', cscale=cscale, stylize=True)
    maxima_trace = mgr.plot_maxima(dgrid, maxima, mat_atoms)
    windows_trace = mgr.plot_windows(connections_trimmed, marker_color='purple')
    connections_traces = mgr.plot_connections_as_it_is(connections_trimmed, all_maxima, dgrid.shape, mat_atoms)

    # Create Subfigure A: Unit cell and framework
    logger.info("Creating Subfigure A: Unit cell and framework...")
    subfig_a_traces = [*frame_trace, cell_trace]
    create_subfigure(subfig_a_traces, '/Users/arungopalan/repos/mofography/subfigure_A_framework.svg',
                    width=600, height=450, scale=2)

    # Create Subfigure B: Unit cell maxima surface colored by regions
    logger.info("Creating Subfigure B: Unit cell maxima surface colored by regions...")
    subfig_b_traces = [surface_trace_regions, maxima_trace]
    create_subfigure(subfig_b_traces, '/Users/arungopalan/repos/mofography/subfigure_B_regions.svg',
                    width=600, height=450, scale=2)

    # Create Subfigure C: Unit cell maxima surface with regions, connections, and windows
    logger.info("Creating Subfigure C: Unit cell maxima surface with regions, connections, and windows...")
    subfig_c_traces = [surface_trace_regions_pbc, *connections_traces, maxima_trace, windows_trace]
    create_subfigure(subfig_c_traces, '/Users/arungopalan/repos/mofography/subfigure_C_network.svg',
                    width=600, height=450, scale=2)

    logger.info("All subfigures generated successfully!")
    logger.info("Files created:")
    logger.info("  - subfigure_A_framework.svg")
    logger.info("  - subfigure_B_regions.svg")
    logger.info("  - subfigure_C_network.svg")

if __name__ == "__main__":
    main()
"""A powerful python package for the analysis, visualization, featurization and machine learning on molecular pores in nanoporous materials."""

# Add imports here
from .compute import *
from .analyze import *
from .visualize import *
from .helper import *
from .data import *
from .workflow import *
from .raspa_helper import *
from .rust_fun import *
from .winding_percolation import *
from .movies import *
from .threejs_viz import (
    Figure3D,
    Figure3js,
    Figure3JS,
    MovieFigure,
    add_trace as add_trace_3js,
    add_traces as add_traces_3js,
    plot_surface_dgrid_3js,
    plot_cell_3js,
    draw_framework_3js,
    plot_framework_as_mesh_3js,
    plot_maxima_3js,
    plot_maxima_from_grid_3js,
    plot_connections_3js,
    plot_connections_as_it_is_3js,
    plot_windows_3js,
    plot_windows_from_results_3js,
    plot_percolated_cluster_meshes_3js,
    plot_percolation_windows_3js,
)

from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

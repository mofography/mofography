"""
 Functions for creating 
 ██╗    ██╗ ██████╗ ██████╗ ██╗  ██╗███████╗██╗      ██████╗ ██╗    ██╗
 ██║    ██║██╔═══██╗██╔══██╗██║ ██╔╝██╔════╝██║     ██╔═══██╗██║    ██║
 ██║ █╗ ██║██║   ██║██████╔╝█████╔╝ █████╗  ██║     ██║   ██║██║ █╗ ██║
 ██║███╗██║██║   ██║██╔══██╗██╔═██╗ ██╔══╝  ██║     ██║   ██║██║███╗██║
 ╚███╔███╔╝╚██████╔╝██║  ██║██║  ██╗██║     ███████╗╚██████╔╝╚███╔███╔╝
  ╚══╝╚══╝  ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚══════╝ ╚═════╝  ╚══╝╚══╝ 
                                                                       
 
"""


def build_padded_threejs_scene(cif_path, output_path=None, **kwargs):
    """Build a Three.js HTML scene with padded transparent segmentation surfaces.

    Thin wrapper that delegates to the stand-alone script's builder so it
    stays in one place.  Any keyword argument is forwarded to the builder
    (``spacing``, ``isovalue``, ``mask_thickness``, ``h``, ``threshold_abs``,
    ``min_distance``, ``extra_padding``, ``preserve_threshold``, …).
    """
    import sys
    from pathlib import Path
    repo_root = Path(__file__).resolve().parent.parent
    scripts_dir = repo_root / "scripts"
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))
    from make_padded_threejs_scene import build_padded_threejs_scene as _builder
    return _builder(cif_path=cif_path, output_path=output_path, **kwargs)


def quick_cluster(ase_atoms,

                  spacing=0.5,
                  mask_thickness=1.5,
                  h=0.2,
                  compactness=0,
                  zoom_ratio=2,
                  wall_thickness=1.0,
                  isoval=1.75,
                  bins=None,
                  block_size=50000, a=1, b=1, c=1, 
                  density=True,
                  js_threshold = 0.15,
                  distance_threshold=0.4, minimum_overlap_fraction=0.5, make_figures=True, renderer='browser', sample_fraction=0.5, cscale=None,return_these_vars=[]):  # distance threshold for agglomerative clustering
    # # # * Pick a random material from the database
    # if chosen is None:
    #     chosen = random.choice(np.arange(len(iza_atoms)))
    # # chosen = 108
    # # chosen = 107
    # print(chosen, df['name'][chosen])
    # mat_atoms = iza_atoms[chosen]
    # from ase.io import read
    # mat_atoms = read('../examples/TestCIFs/NU-1000_S1.cif')
    """Run the legacy Python dgrid segmentation and geometry-JS clustering pipeline with optional figures.

    :param ase_atoms: ASE atoms object
    :param spacing: the spacing of the distance grid
    :param mask_thickness: The thickness of the mask used to create the distance grid
    :param h: float, optional
    :param compactness: float, optional
    :param zoom_ratio: The distance grid is made finer by zooming in, here we are making it twice as
    fine, defaults to 2 (optional)
    :param wall_thickness: The thickness of the walls of the material
    :param bins: the bins for the cavity size distribution histogram
    :param block_size: the size of the blocks in the dask array, this is the size of the chunks that are
    computed in parallel, defaults to 50000 (optional)
    :param a: number of unit cells in the a direction, defaults to 1 (optional)
    :param b: int, defaults to 1 (optional)
    :param c: number of unit cells in the c direction, defaults to 1 (optional)
    :param distance_threshold: distance threshold for agglomerative clustering
    :param minimum_overlap_fraction: The minimum fraction of overlap between the projections of two regions on opposite walls of the unit cell for them to be
    considered as part of the same pbc group
    :param show_figures: whether to show the figures or not, defaults to True (optional)
    :param renderer: 'browser' or 'png', defaults to browser (optional)
    :param list_return_strings: list of strings, the names of the variables you want to return
    """


    

    import numpy as np 

    #* to find and return variables 
    import inspect

    #* for atomic radii
    from ase.data import covalent_radii

    # * cif manipulation
    import ase
    from ase.build import make_supercell

    # * For manipulating the distance grid efficiently
    import dask.array as da

    #* make grid finer
    from scipy.ndimage import zoom 

    #* measure some properties ofthe regions 
    from skimage.measure import regionprops

    #* Plotting and visualization
    import plotly.graph_objects as go
    # Import necessary libraries for plotting and analysis
    import plotly.express as px 
    from .analyze import update_with_covalent_bonds
    from .compute import (dgrid_from_atoms_cpu, regions_from_dgrid, apply_pbc, 
                           connections_from_regions_and_dgrid, rag_from_connections, 
                           add_pbc_to_rag, add_centroids_to_rag_with_pbc, 
                           add_psd_to_rag_with_pbc, clustering_grid_from_regions, 
                           add_clustering_to_rag_geometry_js, add_maxima_to_rag)
    from .data import mix_lorentz_berthelot
    from .visualize import (plot_surface_dgrid, plot_framework_as_mesh, 
                            plot_cell, plot_maxima, plot_csds, 
                            white_layout, plot_vdw_hists)
    from .helper import zoom_grid

    # Set the colorscale for plotting if not already defined
    if cscale is None:
        # Use Plotly's qualitative colorscale for consistency in figures
        cscale = px.colors.qualitative.Plotly

    # Define bins for the cavity feature histogram if not provided
    if bins is None:
        # Create 50 evenly spaced bins between 0 and 30
        bins = np.linspace(0, 30, 50)

    # Update the atomic structure with covalent bonds information
    mat_updated = update_with_covalent_bonds(ase_atoms)

    # Generate a distance grid from the atomic structure
    dgrid = dgrid_from_atoms_cpu(ase_atoms, spacing=spacing, precompute=True)

    # Extract regions and maxima from the distance grid
    regions, maxima = regions_from_dgrid(dgrid, mask_thickness=mask_thickness, 
                                         compactness=compactness, h=h)

    # Get the radii of the maxima from the distance grid
    maxima_radii = dgrid[tuple(maxima.T)]

    # Apply periodic boundary conditions (PBC) to the regions and maxima
    pbc_groups = apply_pbc(regions, maxima, maxima_radii, ase_atoms, 
                            wall_thickness=wall_thickness, return_conn=True, 
                            minimum_overlap_fraction=minimum_overlap_fraction)  # Returns a list of PBC groups

    # Get the regions with PBC applied, returning a 3D array
    regions_pbc = apply_pbc(regions, maxima, maxima_radii, ase_atoms, 
                             wall_thickness=wall_thickness, return_conn=False, 
                             minimum_overlap_fraction=minimum_overlap_fraction)

    # Create connections from regions and the distance grid
    connections = connections_from_regions_and_dgrid(regions, dgrid, maxima,ase_atoms)

    # Generate a region adjacency graph (RAG) from the connections
    rag = rag_from_connections(connections, maxima)

    # Add maxima information to the RAG
    rag = add_maxima_to_rag(rag, maxima, maxima_radii, regions.shape, ase_atoms)

    # Measure properties of the regions
    props = regionprops(regions)

    # Incorporate PBC information into the RAG
    rag = add_pbc_to_rag(rag, pbc_groups)   

    # Add centroids to the RAG with PBC consideration
    rag = add_centroids_to_rag_with_pbc(rag, props, pbc_groups,ase_atoms, regions.shape)

    # Add cavity size distribution (PSD) information to the RAG with PBC
    rag = add_psd_to_rag_with_pbc(rag, props, pbc_groups, bins=bins, 
                                   ase_atoms=ase_atoms, shape=regions.shape)

    # Add clustering information to the RAG geometry
    rag =  add_clustering_to_rag_geometry_js(rag, props, pbc_groups, ase_atoms, dgrid.shape)

    # Generate clusters from the regions based on the RAG
    regions_cluster =  clustering_grid_from_regions(regions, rag, cluster_key='cluster_geom_pbc')

    if make_figures:
        # vdw_traces = plot_vdw_hists(rag, pbc_groups, cscale=cscale, cluster_key='cluster_geom_pbc')
        csd_traces = plot_csds(rag, pbc_groups, cscale=cscale, cluster_key='cluster_geom_pbc')
        fig_csd = go.Figure(data=csd_traces)
        fig_csd.update_layout(xaxis_title='Distance to centroid (A)', yaxis_title='Normalized histogram')
        # fig_vdw.show(renderer=renderer)

    if make_figures:


        # * Plot the distance isosurface at 1.75 A colored in different ways
        surface_trace =  plot_surface_dgrid(dgrid,isoval, select='all', ase_atoms=mat_updated,
                                               region_labels=regions, cscale=cscale, stylize=True, name='Raw segmentation', hover_label='Region :')
        surface_trace_pbc =  plot_surface_dgrid(dgrid,isoval, select='all', ase_atoms=mat_updated,
                                                   region_labels=regions_pbc-1, cscale=cscale, stylize=True, name='PBC groups', hover_label='PBC Group :')
        


        surface_trace_cluster =  plot_surface_dgrid(dgrid,isoval, select='all', ase_atoms=mat_updated,
                                                       region_labels=regions_cluster-1, cscale=cscale, stylize=True, name='Cavity types', hover_label='Cavity type :')

        # * plot the framework atoms and bonds as a mesh
        frame_trace =  plot_framework_as_mesh(mat_updated, show_rings=True, stylize=True)
        # surface_trace= plot_surface_dgrid(grid.compute(), 1.75, select=list(pbc_groups[0]),ase_atoms=mat_updated,region_labels=regions, stylize=True)

        # * plot the unit cell as a box
        cell_trace =  plot_cell(mat_updated)

        # * plot the local maxima as points colored by their size (radius)
        maxima_trace =  plot_maxima(maxima=maxima, dgrid=dgrid, ase_atoms=mat_updated)

        #* plot the windows of the material
        # connections =  get_connections_for_rag(rag, region_labels=regions, ase_atoms = mat_atoms, maxima=maxima, dgrid=grid.compute())
        # windows_trace =  plot_windows(connections)
        # * collect the traces into a figure object
        fig_full_visual = go.Figure(data=[surface_trace, surface_trace_pbc,
                                    surface_trace_cluster, *frame_trace, cell_trace, maxima_trace])
        # fig_full_visual = go.Figure(data=[surface_trace, surface_trace_pbc, surface_trace_cluster, *frame_trace, cell_trace, maxima_trace, windows_trace])
        fig_full_visual.update(layout= white_layout)
        fig_full_visual.update_layout(scene=dict(aspectmode='data'))
        # fig_full_visual.show(renderer=renderer)

        
    # # # * compute and plot the energy histograms for each pbc group
    # rag =  add_vdw_hist_to_rag(rag, props, mat_atoms, forcefield_mixed, regions.shape,use_pbc=True,  pbc_groups=pbc_groups, density=True)
    # ehist_traces =  plot_vdw_hists(rag, pbc_groups, cscale=cscale)
    # fig_ehist = go.Figure(data=ehist_traces)
    # fig_ehist.update_layout(xaxis_title='Energy (K)', yaxis_title='Normalized histogram')
    # fig_ehist.show(renderer=renderer)
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    # callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    if len(return_these_vars) > 0:
        print('Some return objects specified, possible options to pass in list_return_strings are:')
        # list_return_objects = [val for name, val in callers_local_vars if name in return_these_vars]
        # list_return_objects = [val for i,val in enumerate(values) if args[i] in return_these_vars]
        list_return_objects = [frame.f_locals[var] for var in return_these_vars]

        return list_return_objects + [[fig_full_visual, fig_csd]]
    else:
        print('No return objects specified, possible options to pass in list_return_strings are:')
        # print(callers_local_vars)
        print(frame.f_locals.keys())
        if make_figures:
            print('Returning the figures if make_figures=True')
            return [fig_full_visual, fig_csd]

# def quick_cluster_copy(ase_atoms,

#                   spacing=0.5,
#                   mask_thickness=1.5,
#                   h=0.2,
#                   compactness=0,
#                   zoom_ratio=2,
#                   wall_thickness=1.0,
#                   isoval=1.75,
#                   bins=None,
#                   block_size=50000, a=1, b=1, c=1, 
#                   density=True,
#                   distance_threshold=0.4, minimum_overlap_fraction=0.5, make_figures=True, renderer='browser', use_energy_for_clustering=False, sample_fraction=0.5, cscale=None,return_these_vars=[]):  # distance threshold for agglomerative clustering
#     # # # * Pick a random material from the database
#     # if chosen is None:
#     #     chosen = random.choice(np.arange(len(iza_atoms)))
#     # # chosen = 108
#     # # chosen = 107
#     # print(chosen, df['name'][chosen])
#     # mat_atoms = iza_atoms[chosen]
#     # from ase.io import read
#     # mat_atoms = read('../examples/TestCIFs/NU-1000_S1.cif')
#     """ 
#     It takes in a material (in the form of an ASE atoms object), and returns a list of objects that you
#     can specify in the list_return_strings argument.
    
#     :param ase_atoms: ASE atoms object
#     :param spacing: the spacing of the distance grid
#     :param mask_thickness: The thickness of the mask used to create the distance grid
#     :param h: float, optional
#     :param compactness: float, optional
#     :param zoom_ratio: The distance grid is made finer by zooming in, here we are making it twice as
#     fine, defaults to 2 (optional)
#     :param wall_thickness: The thickness of the walls of the material
#     :param bins: the bins for the cavity size distribution histogram
#     :param block_size: the size of the blocks in the dask array, this is the size of the chunks that are
#     computed in parallel, defaults to 50000 (optional)
#     :param a: number of unit cells in the a direction, defaults to 1 (optional)
#     :param b: int, defaults to 1 (optional)
#     :param c: number of unit cells in the c direction, defaults to 1 (optional)
#     :param distance_threshold: distance threshold for agglomerative clustering
#     :param minimum_overlap_fraction: The minimum fraction of overlap between the projections of two regions on opposite walls of the unit cell for them to be
#     considered as part of the same pbc group
#     :param show_figures: whether to show the figures or not, defaults to True (optional)
#     :param renderer: 'browser' or 'png', defaults to browser (optional)
#     :param list_return_strings: list of strings, the names of the variables you want to return
#     """


    

#     import numpy as np 

#     #* to find and return variables 
#     import inspect

#     #* for atomic radii
#     from ase.data import covalent_radii

#     # * cif manipulation
#     import ase

#     # * For manipulating the distance grid efficiently
#     import dask.array as da

#     #* make grid finer
#     from scipy.ndimage import zoom 

#     #* measure some properties ofthe regions 
#     from skimage.measure import regionprops

#     #* Plotting and visualization
#     import plotly.graph_objects as go

#     #* colorscale for plotting 
#     import plotly.express as px 
#     from .analyze import update_with_covalent_bonds
#     from .compute import dgrid_from_atoms,regions_from_dgrid,apply_pbc, get_rag_from_regions_and_grid,add_pbc_to_rag, add_psd_to_rag, add_psd_to_rag_with_pbc, add_clustering_to_rag_pca, clustering_grid_from_regions,add_clustering_to_rag_energy, add_vdw_hist_to_rag
#     from .data import mix_lorentz_berthelot
#     from .visualize import plot_surface_dgrid, plot_framework_as_mesh, plot_cell, plot_maxima, plot_csds, white_layout, plot_vdw_hists, draw_framework
#     if cscale is None:
#         #* set colorscale for plotting, make things consistent across the two figures
#         cscale = px.colors.qualitative.Alphabet

#     #* bins for the cavity feature histogram
#     if bins is None:
#         bins = np.linspace(0,30,50)

#     mat_updated =  update_with_covalent_bonds(ase_atoms)

#     # * compute the distance grid, locate the local maxima and make a list of their size (radii)
#     radii = covalent_radii[ase_atoms.get_atomic_numbers()]  # assign covalent radii
#     # compute the distance grid (dask array) of given spacing and partition size
#     grid =  dgrid_from_atoms(ase_atoms, radii=radii, spacing=spacing, block_size=block_size)
#     # (optional) make a supercell of the material
#     ase_atoms = ase.build.make_supercell(ase_atoms, [[a, 0, 0], [0, b, 0], [0, 0, c]])
#     grid = da.from_array(np.tile(grid.compute(), [a, b, c]))
#     # (optional) make the distance grid finer by zooming in, here we are aamking it twice as fine
#     grid = da.from_array(zoom(grid.compute(), zoom_ratio, mode='grid-wrap',
#                          order=1, grid_mode=False), chunks=(40, 40, 40))
#     # compute the region segementaion and array indices of the maxima corresponding to each region
#     regions, maxima =  regions_from_dgrid(
#         grid.compute(), mask_thickness=mask_thickness, compactness=compactness, h=h)
#     maxima_radii = grid.compute()[tuple(maxima.T)]  # get teh distance values (size) of the maxima

#     #* get the connectivity of the pockets, each pbc group is a list of regions that are split up by the PBC, that come together to form a complete cavity.
#     # * return_conn = True returns the list of pbc groups, return_conn = False returns the a 3D array of regions where region index is replaced by the pbc group index
#     pbc_groups =  apply_pbc(regions, maxima, maxima_radii, ase_atoms, wall_thickness=wall_thickness,
#                                return_conn=True, minimum_overlap_fraction=minimum_overlap_fraction)  # this returns a list of pbc groups
#     regions_pbc =  apply_pbc(regions, maxima, maxima_radii, ase_atoms, wall_thickness=wall_thickness,
#                                 return_conn=False, minimum_overlap_fraction=minimum_overlap_fraction)  # this returns a 3D array
#     # * construct the region adjacency graph, a.k.a the pore graph form the regions, grid, maxima indices and the ase atoms object
#     rag =  get_rag_from_regions_and_grid(regions, grid.compute(), ase_atoms, maxima)

#     # * in this part we will be adding properties to the region adjacency graph
#     rag =  add_pbc_to_rag(rag, pbc_groups)
#     props = regionprops(regions)  # * get the region properties from skimage.measure.regionprops
#     rag =  add_psd_to_rag(rag, region_props=props,  bins=bins, ase_atoms=ase_atoms,
#                               shape=regions.shape)  # *  add the cavity histogram without pbc to the RAG
#     rag =  add_psd_to_rag_with_pbc(rag=rag, props=props, pbc_groups=pbc_groups, bins=bins,
#                                        ase_atoms=ase_atoms, shape=regions.shape)  # * add the cavity histogram with pbc to the RAG
#     # rag=  add_clustering_to_rag(rag) # * add the clustering labels to the RAG without pbc

#     # * add energy histogram to the RAG
#     # forcefield_mixed =  mix_lorentz_berthelot(mat_atoms)
#     # rag =  add_vdw_hist_to_rag_with_pbc(rag, props, pbc_groups, mat_atoms, forcefield_mixed, regions.shape)
#     if use_energy_for_clustering:
#         # * use the energy criterion for clustering
#         props = regionprops(regions) # * get the  region properties from skimage.measure.regionprops
#         forcefield_mixed = mix_lorentz_berthelot(ase_atoms) # * mix the lorentz-berthelot forcefield using UFF and helium probe
#         # * do the clustering using the sampled point energies with pbc, 
#         # * using the kolmogorov-smirnov two-sided test with 95% confidence interval to tell if the two distributions are same or not
#         rag = add_vdw_hist_to_rag(rag, props, ase_atoms, forcefield_mixed, regions.shape, use_pbc=True, pbc_groups=pbc_groups, sample_fraction=sample_fraction, density=density)
#         rag = add_clustering_to_rag_energy(rag=rag, props=props, pbc_groups=pbc_groups, ase_atoms=ase_atoms, forcefield_mixed =forcefield_mixed, shape=regions.shape, use_pbc=True, sample_fraction=sample_fraction)
#         regions_cluster =  clustering_grid_from_regions(regions, rag, cluster_key='cluster_energy_pbc')
#         if make_figures:
#             vdw_traces = plot_vdw_hists(rag, pbc_groups, cscale=cscale, cluster_key='cluster_energy_pbc')
#             fig_csd = go.Figure(data=csd_traces)
#             fig_vdw.update_layout(xaxis_title='Distance to centroid (A)', yaxis_title='Normalized histogram')
#             # fig_vdw.show(renderer=renderer)
#     else:
#         pass
#         # rag, cluster_obj =  add_clustering_to_rag_pca(rag, no_pbc=False, kwargs_cluster=dict(distance_threshold=distance_threshold, kernel='linear'), return_cluster_obj=True, kwargs_pca=dict(
#         #     kernel='cosine'))  # * do clustering of the regions with PCA and default settings (clustering kernel='linear', linkage='single', distance='euclidean', pca kernel ='cosine')

#         # #---------------------------------------------------#
#         # # #* add energy histogram to the RAG
#         # # props = regionprops(regions)
#         # # forcefield_mixed =  mix_lorentz_berthelot(mat_atoms)
#         # # rag =  add_vdw_hist_to_rag(rag, props, mat_atoms, forcefield_mixed, regions.shape, use_pbc=True,pbc_groups=pbc_groups,density=False)

#         # # #* add clustering to rag using energy with pbc
#         # # rag, cluster_obj  = add_clustering_to_rag_pca2(rag, no_pbc=False, kwargs_cluster=dict(distance_threshold=distance_threshold, kernel = 'linear'), return_cluster_obj=True, kwargs_pca=dict(kernel='cosine')) # * do clustering of the regions with PCA and default settings (clustering kernel='linear', linkage='single', distance='euclidean', pca kernel ='cosine')
#         # #-------------------------------------------------#
#         # # * get the 3D array of cluster labels
#         # regions_cluster =  clustering_grid_from_regions(regions, rag)
#         if make_figures:
#             # * plot the cavity histogram for each pbc group
#             csd_traces =  plot_csds(rag, pbc_groups, cscale=cscale)
#             fig_vdw = go.Figure(data=csd_traces)
#             fig_vdw.update_layout(xaxis_title='Distance to centroid (A)', yaxis_title='Normalized histogram')
#             # fig_vdw.show(renderer=renderer)

#     if make_figures:


#         # * Plot the distance isosurface at 1.75 A colored in different ways
#         surface_trace =  plot_surface_dgrid(grid.compute(),isoval, select='all', ase_atoms=mat_updated,
#                                                region_labels=regions, cscale=cscale, stylize=True, name='Raw segmentation', hover_label='Region :')
#         surface_trace_pbc =  plot_surface_dgrid(grid.compute(),isoval, select='all', ase_atoms=mat_updated,
#                                                    region_labels=regions_pbc-1, cscale=cscale, stylize=True, name='PBC groups', hover_label='PBC Group :')
        

#         # ## * mesh 
#         # # * Compute the isosurface using Marching Cubes Lewiner as usual
#         # from skimage.measure import marching_cubes as mc
#         # verts, faces, normals, values = mc(grid.compute(), isoval, step_size=1)

#         # # * create a dictionary of colors for the cluster types
#         # import networkx as nx
#         # unique_clusters = np.unique(list(nx.get_node_attributes(rag, 'cluster').values())).astype(int)
#         # cmap = dict(zip(unique_clusters, cscale[:len(unique_clusters)]))
#         # cmap[-1]=cscale[-1]


#         # # * create an interpolator for the cluster grid
#         # from scipy.interpolate import RegularGridInterpolator as RGI
#         # xx = np.linspace(0,regions_cluster.shape[0] - 1,regions_cluster.shape[0])
#         # yy = np.linspace(0,regions_cluster.shape[1] - 1,regions_cluster.shape[1])
#         # zz = np.linspace(0,regions_cluster.shape[2] - 1,regions_cluster.shape[2])
#         # rl = RGI((xx, yy, zz),regions_cluster-1, method='nearest') # starts from -1 for backgroudn and then 0, 1, 2 etc  for pocket types
        

#         # colours = rl(verts)
#         # #  Find the actual coordinates of the vertices for plotting
#         # from ase.geometry import complete_cell
#         # hovertext = ['Pocket type:' + str(int(c)) for c in colours]
#         # cell = complete_cell(mat_updated.get_cell())
#         # true_verts = np.dot(cell, (verts / regions_cluster.shape).T).T
#         # vertexcolors = [cmap[c] for c in colours]
#         # surface_trace_cluster = go.Mesh3d(x=true_verts[:, 0], y=true_verts[:, 1], z=true_verts[:, 2], vertexcolor=vertexcolors,
#         #                       hoverinfo='all', hovertext=hovertext, flatshading=False, colorscale=cscale, i=faces[:, 0],
#         #                       j=faces[:, 1], k=faces[:, 2], opacity=1, name='Pocket types', showlegend=False)


#         surface_trace_cluster =  plot_surface_dgrid(grid.compute(),isoval, select='all', ase_atoms=mat_updated,
#                                                        region_labels=regions_cluster-1, cscale=cscale, stylize=True, name='Cavity types', hover_label='Cavity type :')

#         # * plot the framework atoms and bonds as a mesh
#         # frame_trace =  plot_framework_as_mesh(mat_updated, show_rings=True, stylize=True)
#         frame_trace =draw_framework(mat_updated, show_aromatics=True, bond_color='gray')
#         # surface_trace= plot_surface_dgrid(grid.compute(), 1.75, select=list(pbc_groups[0]),ase_atoms=mat_updated,region_labels=regions, stylize=True)

#         # * plot the unit cell as a box
#         cell_trace =  plot_cell(mat_updated)

#         # * plot the local maxima as points colored by their size (radius)
#         maxima_trace =  plot_maxima(maxima=maxima, dgrid=grid.compute(), ase_atoms=mat_updated)

#         #* plot the windows of the material
#         # connections =  get_connections_for_rag(rag, region_labels=regions, ase_atoms = mat_atoms, maxima=maxima, dgrid=grid.compute())
#         # windows_trace =  plot_windows(connections)
#         # * collect the traces into a figure object
#         fig_full_visual = go.Figure(data=[surface_trace, surface_trace_pbc,
#                                     surface_trace_cluster, *frame_trace, cell_trace, maxima_trace])
#         # fig_full_visual = go.Figure(data=[surface_trace, surface_trace_pbc, surface_trace_cluster, *frame_trace, cell_trace, maxima_trace, windows_trace])
#         fig_full_visual.update(layout= white_layout)
#         fig_full_visual.update_layout(scene=dict(aspectmode='data'))
#         # fig_full_visual.show(renderer=renderer)

        
#     # # # * compute and plot the energy histograms for each pbc group
#     # rag =  add_vdw_hist_to_rag(rag, props, mat_atoms, forcefield_mixed, regions.shape,use_pbc=True,  pbc_groups=pbc_groups, density=True)
#     # ehist_traces =  plot_vdw_hists(rag, pbc_groups, cscale=cscale)
#     # fig_ehist = go.Figure(data=ehist_traces)
#     # fig_ehist.update_layout(xaxis_title='Energy (K)', yaxis_title='Normalized histogram')
#     # fig_ehist.show(renderer=renderer)
#     frame = inspect.currentframe()
#     args, _, _, values = inspect.getargvalues(frame)
#     # callers_local_vars = inspect.currentframe().f_back.f_locals.items()
#     if len(return_these_vars) > 0:
#         print('Some return objects specified, possible options to pass in list_return_strings are:')
#         # list_return_objects = [val for name, val in callers_local_vars if name in return_these_vars]
#         # list_return_objects = [val for i,val in enumerate(values) if args[i] in return_these_vars]
#         list_return_objects = [frame.f_locals[var] for var in return_these_vars]

#         return list_return_objects # + [[fig_full_visual, fig_vdw]]
#     else:
#         print('No return objects specified, possible options to pass in list_return_strings are:')
#         # print(callers_local_vars)
#         print(frame.f_locals.keys())
#         if make_figures:
#             print('Returning the figures if make_figures=True')
#             return [fig_full_visual, fig_csd]
def enrich_segmentation_with_pocket_psds(
    segmentation,
    atoms,
    *,
    inplace=True,
    raise_on_error=False,
    probe_radius_psd=None,
    **kwargs,
):
    """Run pocket-wise Rust PSDs once and store on the segmentation dict.

    Call before :func:`assign_pocket_types` and/or
    :func:`~mofography.compute.add_pocket_psds_to_rag_from_segmentation`.
    """
    from mofography.pocket_typing import DEFAULT_PROBE_RADIUS_PSD_A
    from mofography.rust_fun import compute_pocket_psds_for_segmentation

    pr = (
        float(probe_radius_psd)
        if probe_radius_psd is not None
        else DEFAULT_PROBE_RADIUS_PSD_A
    )
    try:
        pocket_psds = compute_pocket_psds_for_segmentation(
            segmentation,
            atoms,
            probe_center_min_radius=pr,
            **kwargs,
        )
        pocket_psd_error = None
    except Exception as exc:
        pocket_psd_error = str(exc)
        if raise_on_error:
            raise
        pocket_psds = None

    patch = {
        "pocket_psd_probe_radius_a": pr,
        "pocket_psd_error": pocket_psd_error,
    }
    if pocket_psds is not None:
        patch["pocket_psds"] = pocket_psds

    if inplace:
        segmentation.update(patch)
        return segmentation
    return {**dict(segmentation), **patch}


def assign_pocket_types(segmentation_results, atoms, method='env_cluster', **kwargs):
    """Assign pocket types to segmentation results.

    Wraps pocket_typing.assign_pocket_types_from_segmentation.

    Reuses ``segmentation_results['pocket_psds']`` when set (see
    :func:`enrich_segmentation_with_pocket_psds`).

    By default ``typing_radii='bondi_crc_fallback'``: ``max_r`` bands use a
    second Bondi/CRC distance grid while segmentation ``dgrid`` (CRC/Mendeleev/
    user radii) still drives PSD, windows, and framework fingerprints.  Pass
    ``typing_radii='segmentation'`` for the legacy behaviour where typing reads
    the same ``dgrid`` as the rest of the analysis.

    :param segmentation_results: dict from do_periodic_segmentation_rust
    :param atoms: ASE Atoms object
    :param method: 'env_cluster' (default, SH fingerprint + connectivity)
    :param typing_radii: VDW table for typing max_r only (see pocket_typing)
    :return: PocketTypingResult
    """
    from mofography.pocket_typing import assign_pocket_types_from_segmentation
    return assign_pocket_types_from_segmentation(segmentation_results, atoms, method=method, **kwargs)

def add_percolation_windows_to_fig(fig, analysis_result, pld_probe_radius=None, pld_tolerance=0.25):
    """Add window traces (regular + PLD-critical diamonds) to an existing figure.

    :param fig: plotly Figure (go.Figure)
    :param analysis_result: dict with 'connections' key from segmentation
    :param pld_probe_radius: if given, highlight PLD-critical windows as red diamonds
    :param pld_tolerance: radius tolerance for PLD-critical classification (default 0.25)
    :return: figure with window traces added
    """
    from mofography.visualize import plot_percolation_windows
    traces = plot_percolation_windows(
        analysis_result,
        probe_radius=pld_probe_radius,
        pld_tolerance=pld_tolerance,
    )
    for t in traces:
        fig.add_trace(t)
    return fig


def pocket_regions_cluster(typed_result, seg):
    """Return the pocket-type region label grid for surface colouring.

    ``assign_pocket_types`` stores the typed region array as
    ``PocketTypingResult.regions_cluster`` (a ``numpy.ndarray``).  The naive
    ``typed.regions_cluster or seg.get("regions_pbc")`` fails with
    ``ValueError`` because numpy arrays have no single boolean truth value.

    This helper does the safe ``is not None`` check and falls back to the raw
    PBC-merged region grid from the segmentation dict when pocket typing did
    not populate the array (e.g. when ``result.ok`` is False).

    Usage::

        from mofography.workflow import assign_pocket_types, pocket_regions_cluster

        result  = assign_pocket_types(seg, atoms)
        regions = pocket_regions_cluster(result, seg)
        surface = plot_surface_dgrid(dgrid, 1.75, select="all",
                                     region_labels=regions, ...)

    :param typed_result: ``PocketTypingResult`` from ``assign_pocket_types``
    :param seg: segmentation dict from ``do_periodic_segmentation_rust``
    :return: numpy int array with one label per voxel, or None
    """
    rc = getattr(typed_result, "regions_cluster", None)
    if rc is not None:
        return rc
    return seg.get("regions_pbc")


def _cache_aromatic_cycles_if_needed(fw_atoms, *, show_rings: bool) -> None:
    """Compute aromatic cycles once and cache on ``fw_atoms.info`` (shared by Plotly/Three.js)."""
    if not show_rings or fw_atoms is None:
        return
    from mofography.analyze import compute_aromatic_cycles

    compute_aromatic_cycles(fw_atoms)


def show_mof_plotly(seg, atoms, typed_result=None, stylize=True, framework_atoms=None, **options):
    """Build the full mofography visualization figure (Plotly backend).

    Reproduces the default interactive-app main visualization:

    - Framework mesh with aromatic rings
    - Unit cell
    - Raw region surface
    - PBC-merged region surface
    - Pocket-type surface (when *typed_result* is supplied and ``result.ok``)
    - Pore-graph connections trimmed to accessible (unique) regions
    - **Internal** windows only — red spheres with far-outside reclassification
    - Local maxima

    Parameters
    ----------
    seg : dict
        Output of ``do_periodic_segmentation_rust``.
    atoms : ase.Atoms
        Structure to visualize.
    typed_result : PocketTypingResult, optional
        Output of ``assign_pocket_types``.  When supplied and ``result.ok``,
        adds a pocket-type surface coloured by type id.
    stylize : bool
        Apply PBR lighting tweaks (default ``True``).
    framework_atoms : ase.Atoms, optional
        Atoms object used *only* for framework mesh and aromatic ring rendering.
        Pass the output of ``update_with_covalent_bonds(atoms)`` here to get
        proper ring highlighting, while keeping *atoms* (un-modified) for all
        grid / segmentation / pocket-typing operations.  When omitted, *atoms*
        is used for both.
    **options
        Override any ``plot_*`` flag, e.g.::

            show_mof_plotly(seg, atoms, plot_surface=False, plot_layout="white")

        Recognised keys (all optional, shown with defaults):

        .. code-block:: text

            plot_cell               True
            plot_framework          True
            plot_framework_show_rings True
            plot_maxima             True
            plot_connections        True
            plot_windows            True
            plot_surface            True
            plot_surface_pbc        True
            plot_surface_isovalue   1.0
            plot_surface_opacity    1.0
            plot_surface_pbc_opacity  1.0
            plot_surface_types_opacity 1.0
            plot_layout             "dark"   # or "white"

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import numpy as np
    import plotly.graph_objects as go
    import plotly.express as px
    from mofography.visualize import (
        plot_surface_dgrid, plot_framework_as_mesh, plot_cell,
        plot_maxima, plot_windows, plot_windows_from_results,
        dark_layout, white_layout,
        plot_connections_as_it_is,
    )

    defaults = {
        "plot_cell": True,
        "plot_framework": True,
        "plot_framework_show_rings": True,
        "plot_maxima": True,
        "plot_connections": True,
        "plot_windows": True,
        "plot_surface": True,
        "plot_surface_pbc": True,
        "plot_surface_isovalue": 1.0,
        "plot_surface_opacity": 1.0,
        "plot_surface_pbc_opacity": 1.0,
        "plot_surface_types_opacity": 1.0,
        # Taubin smoothing — match app defaults (intensity=0.05 → 3 iterations)
        "plot_surface_taubin_smooth": True,
        "plot_surface_taubin_iterations": 3,
        "plot_surface_taubin_lambda": 0.2175,
        "plot_surface_taubin_nu": -0.2175,
        # Marching-cubes step size: 1 = full resolution ("medium" / "high" quality)
        # Set plot_surface_step_size=2 for faster / lower-fidelity meshes
        "plot_surface_step_size": 1,
        "plot_windows_reclassify_far_outside": True,
        # Periodic overlays hidden by default (HUD toggles); set True to show on load
        "plot_connections_periodic": False,
        "plot_windows_periodic": False,
        "plot_maxima_pbc": True,
        "plot_layout": "dark",
    }
    defaults.update({k: v for k, v in options.items() if v is not None})

    dgrid         = np.asarray(seg["dgrid"])
    regions       = np.asarray(seg["regions"])
    regions_pbc   = seg.get("regions_pbc")
    all_maxima    = np.asarray(seg.get("all_maxima", []))
    connections   = seg.get("connections", [])
    unique_regions = np.asarray(seg.get("unique_regions", []), dtype=int)

    _taubin_kw = dict(
        taubin_smooth=bool(defaults.get("plot_surface_taubin_smooth", True)),
        taubin_iterations=int(defaults.get("plot_surface_taubin_iterations", 3)),
        taubin_lambda=float(defaults.get("plot_surface_taubin_lambda", 0.2175)),
        taubin_nu=float(defaults.get("plot_surface_taubin_nu", -0.2175)),
    )
    _stepsize = int(defaults.get("plot_surface_step_size", 1))

    # Trim: at least one endpoint in unique_regions (periodic edges with one outside image kept)
    from mofography.helper import trim_connections_to_unique_regions

    connections_trimmed = trim_connections_to_unique_regions(connections, unique_regions)

    # Cap maxima to regions.max() — excludes padding-only sites
    n_cap   = int(regions.max()) if regions.size else 0
    maxima  = all_maxima[:n_cap] if n_cap > 0 else all_maxima[:0]

    traces  = []
    cscale  = px.colors.qualitative.Pastel
    isovalue = float(defaults["plot_surface_isovalue"])

    if defaults.get("plot_framework", True):
        _fw_atoms = framework_atoms if framework_atoms is not None else atoms
        _cache_aromatic_cycles_if_needed(
            _fw_atoms,
            show_rings=bool(defaults.get("plot_framework_show_rings", True)),
        )
        fw = plot_framework_as_mesh(
            _fw_atoms, stylize=stylize,
            show_rings=bool(defaults.get("plot_framework_show_rings", True)),
        )
        traces.extend(fw)

    if defaults.get("plot_cell", True):
        traces.append(plot_cell(atoms))

    if defaults.get("plot_surface", True):
        traces.append(plot_surface_dgrid(
            dgrid=dgrid, isovalue=isovalue, ase_atoms=atoms,
            region_labels=regions,
            opacity=float(defaults["plot_surface_opacity"]),
            name="Regions Surface", cscale=cscale, stylize=stylize,
            stepsize=_stepsize, **_taubin_kw,
        ))

    if defaults.get("plot_surface_pbc", True) and regions_pbc is not None:
        traces.append(plot_surface_dgrid(
            dgrid=dgrid, isovalue=isovalue, ase_atoms=atoms,
            region_labels=np.asarray(regions_pbc),
            opacity=float(defaults["plot_surface_pbc_opacity"]),
            name="Regions Surface (PBC)", cscale=cscale, stylize=stylize,
            stepsize=_stepsize, **_taubin_kw,
        ))

    # Pocket types surface (requires typed_result.ok and .regions_cluster)
    rc = None
    if typed_result is not None and getattr(typed_result, "ok", False):
        rc = getattr(typed_result, "regions_cluster", None)
    if rc is not None:
        traces.append(plot_surface_dgrid(
            dgrid=dgrid, isovalue=isovalue, ase_atoms=atoms,
            region_labels=rc - 1,
            opacity=float(defaults["plot_surface_types_opacity"]),
            name="Pocket Types", cscale=cscale, stylize=stylize,
            hover_label="Pocket type :",
            stepsize=_stepsize, **_taubin_kw,
        ))

    if defaults.get("plot_connections", True) and connections_trimmed:
        conn_traces = plot_connections_as_it_is(
            connections_trimmed, all_maxima, dgrid.shape, atoms,
            plot_periodic=bool(defaults.get("plot_connections_periodic", True)),
        )
        traces.extend(conn_traces)

    # Windows from all trimmed connections; far-outside coords reclassified when enabled.
    if defaults.get("plot_windows", True) and connections_trimmed:
        _cell    = atoms.get_cell().array.T
        _lengths = np.linalg.norm(_cell, axis=0)
        _shape   = np.asarray(dgrid.shape, dtype=float)
        _spacing = float(np.mean(_lengths / np.maximum(_shape, 1.0)))
        _do_reclassify = defaults.get("plot_windows_reclassify_far_outside", True)
        win_traces = plot_windows_from_results(
            connections_trimmed,
            dgrid=dgrid,
            atoms=atoms,
            spacing=_spacing,
            outside_windows_to_periodic=_do_reclassify,
            name="windows",
            marker_color="red",
            marker_size=10,
            opacity=0.6,
            show_periodic=bool(defaults.get("plot_windows_periodic", False)),
        )
        traces.extend(win_traces)

    if defaults.get("plot_maxima", True) and len(maxima):
        mx = plot_maxima(dgrid, maxima, atoms, marker_size=14)
        if isinstance(mx, (list, tuple)):
            traces.extend(mx)
        else:
            traces.append(mx)

    _layout = str(defaults.get("plot_layout", "dark")).lower()
    layout  = white_layout if _layout == "white" else dark_layout
    fig = go.Figure(data=traces)
    fig.update_layout(**layout)
    fig.update_layout(
        scene=dict(aspectmode="data"),
        title=dict(text=getattr(atoms, "info", {}).get("name", "MOF")),
    )
    return fig


def show_mof_threejs(seg, atoms, typed_result=None, stylize=True, framework_atoms=None, output_html=None, show=False, renderer="browser", **options):
    """Build a standalone Three.js HTML visualization.

    Reproduces the default interactive-app main visualization using the
    Three.js / ``MovieFigure`` backend (same as the app's Three.js toggle).

    Includes framework (ball-and-stick + aromatic rings), unit cell, raw and
    PBC-merged region surfaces, optional pocket-type surface, pore-graph
    connections, internal windows (periodic hidden by default), and local
    maxima.

    Parameters
    ----------
    seg : dict
        Output of ``do_periodic_segmentation_rust``.
    atoms : ase.Atoms
        Structure to visualize.
    typed_result : PocketTypingResult, optional
        When supplied and ``result.ok``, adds the pocket-type surface.
    stylize : bool
        Apply PBR roughness/metalness tweaks (default ``True``).
    output_html : str or Path, optional
        If given, write the HTML to this path and return the path string.
        Otherwise return the HTML string directly.
    **options
        Same ``plot_*`` flags as :func:`show_mof_plotly` plus Three.js extras:

        .. code-block:: text

            plot_surface_taubin_smooth      True
            plot_surface_taubin_iterations  10
            plot_surface_taubin_lambda      0.5
            plot_surface_taubin_nu          -0.53
            plot_layout                     "dark"   # or "white"

    Returns
    -------
    str
        Self-contained HTML string (or the output path when *output_html* is set).
    """
    import numpy as np
    import tempfile, os
    from mofography.movies import MovieFigure, make_surface_trace

    defaults = {
        "plot_cell": True,
        "plot_framework": True,
        "plot_framework_show_rings": True,
        "plot_maxima": True,
        "plot_connections": True,
        "plot_windows": True,
        "plot_surface": True,
        "plot_surface_pbc": True,
        "plot_surface_isovalue": 1.0,
        "plot_surface_opacity": 1.0,
        "plot_surface_pbc_opacity": 1.0,
        "plot_surface_types_opacity": 1.0,
        # Taubin smoothing — match app defaults (intensity=0.05 → 3 iterations)
        "plot_surface_taubin_smooth": True,
        "plot_surface_taubin_iterations": 3,
        "plot_surface_taubin_lambda": 0.2175,
        "plot_surface_taubin_nu": -0.2175,
        # Marching-cubes step size: 1 = full resolution ("medium" / "high" quality)
        "plot_surface_step_size": 1,
        # Framework representation (plot_framework_style alias supported):
        #   "mesh"      — baked ball-and-stick mesh (default)
        #   "scatter"   — points + line segments (interactive-app simple scatter)
        #   "instanced" — instanced spheres + cylinders (interactive-app instanced)
        "plot_framework_representation": "mesh",
        "plot_windows_reclassify_far_outside": True,
        # Periodic overlays hidden by default (HUD toggles); set True to show on load
        "plot_connections_periodic": False,
        "plot_windows_periodic": False,
        "plot_maxima_pbc": True,
        "plot_layout": "dark",
    }
    defaults.update({k: v for k, v in options.items() if v is not None})

    dgrid          = np.asarray(seg["dgrid"])
    regions        = np.asarray(seg["regions"])
    regions_pbc    = seg.get("regions_pbc")
    all_maxima_raw = seg.get("all_maxima")
    maxima_pbc     = seg.get("maxima_pbc")
    maxima_info    = seg.get("maxima_info")
    maxima_pbc_dict = seg.get("maxima_pbc_dict", {})
    connections_raw = seg.get("connections", [])
    unique_regions  = np.asarray(seg.get("unique_regions", []), dtype=int)

    from mofography.helper import trim_connections_to_unique_regions

    connections = trim_connections_to_unique_regions(connections_raw, unique_regions)

    # Cap maxima to regions.max()
    n_cap = int(regions.max()) if regions.size else 0
    am = np.asarray(all_maxima_raw, dtype=float) if all_maxima_raw is not None else np.zeros((0, 3), dtype=float)
    all_maxima_plot = am[:n_cap] if n_cap > 0 else am[:0]
    m_plot = len(all_maxima_plot)
    if maxima_info is not None:
        mi = np.asarray(maxima_info, dtype=int)
        maxima_info_plot = mi[:m_plot] if mi.ndim == 2 and mi.shape[0] >= m_plot and m_plot > 0 else np.zeros((0, 7), dtype=int)
    else:
        maxima_info_plot = np.zeros((0, 7), dtype=int)

    # Estimate grid spacing for window reclassification (override with plot_spacing).
    _cell    = atoms.get_cell().array.T
    _lengths = np.linalg.norm(_cell, axis=0)
    _shape   = np.asarray(dgrid.shape, dtype=float)
    spacing  = float(defaults["plot_spacing"]) if defaults.get("plot_spacing") is not None else float(
        np.mean(_lengths / np.maximum(_shape, 1.0))
    )

    fig = MovieFigure()

    taubin_kw = dict(
        taubin_smooth=bool(defaults.get("plot_surface_taubin_smooth", True)),
        taubin_iterations=int(defaults.get("plot_surface_taubin_iterations", 3)),
        taubin_lambda=float(defaults.get("plot_surface_taubin_lambda", 0.2175)),
        taubin_nu=float(defaults.get("plot_surface_taubin_nu", -0.2175)),
    )
    _surf_step_size = int(defaults.get("plot_surface_step_size", 1))
    isovalue = float(defaults["plot_surface_isovalue"])

    if defaults.get("plot_framework", True):
        _fw_atoms = framework_atoms if framework_atoms is not None else atoms
        _cache_aromatic_cycles_if_needed(
            _fw_atoms,
            show_rings=bool(defaults.get("plot_framework_show_rings", True)),
        )
        _fw_repr_raw = defaults.get("plot_framework_representation")
        if _fw_repr_raw is None:
            _fw_repr_raw = defaults.get("plot_framework_style", "mesh")
        _fw_repr = str(_fw_repr_raw).lower()
        _fw_scale = float(defaults.get("plot_framework_mesh_scale", 0.65))
        _fw_bond_lw = float(defaults.get("plot_framework_bond_linewidth", 8.0))
        if _fw_repr == "scatter":
            from mofography.movies import make_framework_scatter_trace
            _scatter = make_framework_scatter_trace(
                _fw_atoms,
                bonds="auto",
                style="framework",
                atom_scale=_fw_scale,
                bond_linewidth=_fw_bond_lw,
            )
            fig.add_trace(_scatter["atoms"])
            fig.add_trace(_scatter["bonds"])
        elif _fw_repr == "instanced":
            from mofography.movies import make_framework_instanced_trace
            _bs = make_framework_instanced_trace(
                _fw_atoms,
                representation="ball_and_stick",
                trace_id="framework",
                atom_scale=_fw_scale,
            )
            fig.add_trace(_bs["atoms"])
            fig.add_trace(_bs["bonds"])
            fig.groups["framework"] = {
                "traces": ["framework_atoms", "framework_bonds"],
                "label": "Framework",
            }
            # Stick variant (hidden by default; toggled by the Representation selector)
            _st = make_framework_instanced_trace(
                _fw_atoms,
                representation="stick",
                trace_id="framework_stick",
                name="Framework (stick)",
                atom_scale=_fw_scale,
            )
            _st["atoms"]["visible"] = False
            _st["bonds"]["visible"] = False
            fig.add_trace(_st["atoms"])
            fig.add_trace(_st["bonds"])
            _sf = make_framework_instanced_trace(
                _fw_atoms,
                representation="space_filling",
                trace_id="framework_spacefill",
                name="Framework (spacefill)",
                atom_scale=_fw_scale,
            )
            _sf["atoms"]["visible"] = False
            fig.add_trace(_sf["atoms"])
        else:
            # Default: baked mesh (ball-and-stick + hidden stick/spacefill variants)
            fig.add_framework(_fw_atoms, representation="ball_and_stick")
        if defaults.get("plot_framework_show_rings", True):
            fig.add_aromatic_rings(_fw_atoms, cycles=_fw_atoms.info.get("aromatic_cycles"))

    if defaults.get("plot_cell", True):
        fig.add_cell(atoms)

    if defaults.get("plot_surface", True):
        t = make_surface_trace(
            dgrid=dgrid, regions=regions, atoms=atoms,
            isovalue=isovalue, name="Regions (raw)",
            trace_id="surface", style="surface",
            hover_prefix="Region:", step_size=_surf_step_size, **taubin_kw,
        )
        fig.add_trace(t)
        fig.update_traces(
            {"opacity": float(defaults["plot_surface_opacity"])},
            selector={"id": "surface"},
        )

    if defaults.get("plot_surface_pbc", True) and regions_pbc is not None:
        t = make_surface_trace(
            dgrid=dgrid, regions=np.asarray(regions_pbc), atoms=atoms,
            isovalue=isovalue, name="Regions (PBC merged)",
            trace_id="surface_regions_pbc", style="surface",
            hover_prefix="PBC:", step_size=_surf_step_size, **taubin_kw,
        )
        fig.add_trace(t)
        fig.update_traces(
            {"opacity": float(defaults["plot_surface_pbc_opacity"])},
            selector={"id": "surface_regions_pbc"},
        )

    # Pocket types surface
    rc = None
    if typed_result is not None and getattr(typed_result, "ok", False):
        rc = getattr(typed_result, "regions_cluster", None)
    if rc is not None:
        t = make_surface_trace(
            dgrid=dgrid, regions=rc, atoms=atoms,
            isovalue=isovalue, name="Pocket types (cluster)",
            trace_id="surface_regions_cluster", style="surface",
            hover_prefix="Pocket type:", step_size=_surf_step_size, **taubin_kw,
        )
        fig.add_trace(t)
        fig.update_traces(
            {"opacity": float(defaults["plot_surface_types_opacity"])},
            selector={"id": "surface_regions_cluster"},
        )

    seg_dict = {
        "dgrid": dgrid,
        "regions": regions,
        "all_maxima": all_maxima_plot,
        "all_maxima_conn": am,
        "maxima_pbc": maxima_pbc if maxima_pbc is not None else np.zeros((0, 3), dtype=int),
        "maxima_info": maxima_info_plot,
        "maxima_info_conn": (
            np.asarray(maxima_info, dtype=int)
            if maxima_info is not None
            else np.zeros((0, 7), dtype=int)
        ),
        "maxima_pbc_dict": maxima_pbc_dict if maxima_pbc_dict else {},
        "connections": connections,
    }

    hide_periodic_conn = not bool(defaults.get("plot_connections_periodic", False))
    hide_periodic_win = not bool(defaults.get("plot_windows_periodic", False))
    hide_periodic_max = not bool(defaults.get("plot_maxima_pbc", True))

    if defaults.get("plot_connections", True) and connections:
        fig.add_connections(seg_dict, atoms, hide_periodic=hide_periodic_conn)

    if defaults.get("plot_windows", True) and connections:
        _win_style = defaults.get("plot_windows_style") or {}
        fig.add_windows(
            seg_dict,
            hide_periodic=hide_periodic_win,
            atoms=atoms,
            spacing=spacing,
            reclassify_far_outside=defaults.get("plot_windows_reclassify_far_outside", True),
            **_win_style,
        )

    if defaults.get("plot_maxima", True) and len(all_maxima_plot):
        fig.add_maxima(seg_dict, atoms, hide_periodic=hide_periodic_max)

    # Legend metadata (optional); HUD legend is schema-driven in the HTML template.
    trace_ids = {t.get("id") for t in fig.data}
    fig.groups = {}

    if stylize:
        fig.update_traces(
            {"roughness": 0.05, "metalness": 0.3},
            selector={"type": "mesh"},
        )

    layout = defaults.get("plot_layout", "dark")
    if layout == "white":
        fig.set_appearance(background="#ffffff", scene_bg=0xffffff)
    else:
        fig.set_appearance(
            background="radial-gradient(circle at center, #7a7a7a 0%, #111111 100%)",
            scene_bg=0x111111,
        )

    # show=True — mirrors plotly's fig.show() / fig.show(renderer="browser")
    # renderer="browser" opens a temp HTML file in the system browser
    # renderer="notebook" displays via IPython.display inside Jupyter
    if show:
        fig.show(renderer=renderer)

    if output_html is not None:
        from pathlib import Path as _P
        _p = _P(output_html)
        _p.parent.mkdir(parents=True, exist_ok=True)
        fig.to_html(str(_p))
        return str(_p)

    # No output_html — return HTML as a string (or None when show was the only goal)
    if not show:
        tmp = tempfile.NamedTemporaryFile(suffix=".html", delete=False)
        try:
            fig.to_html(tmp.name)
            with open(tmp.name, "r", encoding="utf-8") as f:
                return f.read()
        finally:
            os.unlink(tmp.name)

    return None  # show=True, no output_html — browser/notebook opened above


def show_percolation_threejs(
    analysis,
    atoms,
    probe_radius: float,
    output_html=None,
    show: bool = False,
    renderer: str = "browser",
    **options,
):
    """Three.js PLD / percolation scene (app percolation tab parity).

    Thin wrapper around :func:`mofography.threejs_plotting.show_percolation_threejs`.
    """
    from mofography.threejs_plotting import show_percolation_threejs as _build

    html = _build(analysis, atoms, probe_radius=probe_radius, options=options or {})

    if show:
        import tempfile
        import webbrowser

        tmp = tempfile.NamedTemporaryFile(suffix=".html", delete=False)
        tmp.write(html.encode("utf-8"))
        tmp.close()
        webbrowser.open(f"file://{tmp.name}")

    if output_html is not None:
        from pathlib import Path as _P

        _p = _P(output_html)
        _p.parent.mkdir(parents=True, exist_ok=True)
        _p.write_text(html, encoding="utf-8")
        return str(_p)

    return html if not show else None

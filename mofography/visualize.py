"""
 
 Visualization subroutines for mofography
 ██╗   ██╗██╗███████╗██╗   ██╗ █████╗ ██╗     ██╗███████╗███████╗
 ██║   ██║██║██╔════╝██║   ██║██╔══██╗██║     ██║╚══███╔╝██╔════╝
 ██║   ██║██║███████╗██║   ██║███████║██║     ██║  ███╔╝ █████╗  
 ╚██╗ ██╔╝██║╚════██║██║   ██║██╔══██║██║     ██║ ███╔╝  ██╔══╝  
  ╚████╔╝ ██║███████║╚██████╔╝██║  ██║███████╗██║███████╗███████╗
   ╚═══╝  ╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝╚══════╝╚══════╝
                                                                 
 
"""

# TODO draw pore graph
# todo show the region wise properties well ?   



from tenacity import retry

from mofography.analyze import sgraph_from_atoms


global white_layout 
white_layout = dict(template='plotly_white',  scene=dict(xaxis=dict(zeroline=False, visible=False), yaxis_visible=False, zaxis_visible=False), margin=dict(t=20, b=20, l=20, r=20), showlegend=True, legend=dict(x=1.1, y=0.75))
global dark_layout 
dark_layout = dict(template='plotly_dark',  scene=dict(xaxis=dict(zeroline=False, visible=False), yaxis_visible=False, zaxis_visible=False), margin=dict(t=20, b=20, l=20, r=20), showlegend=True, legend=dict(x=0.1, y=1.1))
global epic_layout
epic_layout = dict(template='plotly_dark',  scene=dict(xaxis=dict(zeroline=False, visible=False), yaxis_visible=False, zaxis_visible=False), margin=dict(t=20, b=20, l=20, r=20), showlegend=True, legend=dict(x=0.1, y=0.5, font=dict(size=10),bgcolor='rgba(0,0,0,0)'), plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')


def plot_surface_dgrid(dgrid, isovalue,  ase_atoms,region_labels=None, select='all', opacity=1.0, cscale=None, color='blue',  stepsize=1, stylize=False, random_seed=0):
    """
    It takes a 3D grid, an isovalue, an ASE atoms object, and a few other optional parameters and
    returns a plotly Mesh3D trace object that can be plotted using the plotly library.
    
    :param dgrid: the density grid
    :param isovalue: the value of the isosurface to be plotted
    :param ase_atoms: ASE atoms object
    :param region_labels: a 3D array of the same shape as the dgrid. This array contains the region
    labels for each voxel
    :param select: the region labels to be plotted. If 'all' is selected, all regions are plotted,
    defaults to all (optional)
    :param opacity: opacity of the surface
    :param cscale: a colorscale for the regions. If None, a random colorscale is generated
    :param color: color of the surface, defaults to blue (optional)
    :param stepsize: The step size for the marching cubes algorithm, defaults to 1 (optional)
    :param stylize: if True, will apply some random style choices, purely personal preference, defaults
    to False (optional)
    :param random_seed: int, optional, defaults to 0 (optional)
    :return: A plotly Mesh3D trace
    """


    import numpy as np
    import copy
    import plotly.graph_objects as go
    import ase 
    import plotly.express as px

    dgrid_copy = copy.deepcopy(dgrid) # deepcopy of the grid
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()).T #get cell matrix
    #* Exclude the regions that are not selected
    if select != 'all':
        dgrid_copy[np.isin(region_labels, select, assume_unique=False, invert=True)] = -10

    #* Compute the isosurface using Marching Cubes Lewiner
    from skimage.measure import marching_cubes as mc
    verts, faces, normals, values = mc(dgrid_copy, isovalue, step_size=stepsize)  
    np.random.seed(random_seed)
    if np.any(region_labels) != None: # if region labels are provided

        #*  Make a colorscale using random colors if no colorscale is provided
        if cscale is None:
            cvals = np.linspace(0, 1, (region_labels.max()).astype(int))
            rand_ints = np.random.randint(low=0, high=256, size=(region_labels.max().astype(int),3))
            
            cscale = [[cvals[i],
                    "rgb(" + str(rand_ints[i][0]) + ',' + str(rand_ints[i][1]) + ',' + str(
rand_ints[i][2]) + ')'] for i in range((region_labels.max()).astype(int))]
        if type(cscale) == str:
            try:
                cscale = getattr(px.colors.qualitative, cscale)
            except:
                try:
                    cscale = getattr(px.colors.sequential, cscale)
                except:
                    pass


            # cscale = [[cvals[i],
            #         "rgb(" + str(np.random.randint(0, 256)) + ',' + str(np.random.randint(0, 256)) + ',' + str(
            #             np.random.randint(0, 256)) + ')'] for i in range((region_labels.max()).astype(int))]

        #  Find the colour for the vertices
        from scipy.interpolate import RegularGridInterpolator as RGI
        xx = np.linspace(0,region_labels.shape[0] - 1,region_labels.shape[0])
        yy = np.linspace(0,region_labels.shape[1] - 1,region_labels.shape[1])
        zz = np.linspace(0,region_labels.shape[2] - 1,region_labels.shape[2])
        rl = RGI((xx, yy, zz),region_labels, method='nearest')
        colours = rl(verts)

        #  Create a plotly Mesh3D trace
#        ------------------------------------------

        #  Find the actual coordinates of the vertices for plotting
        hovertext = ["Region: " + str(int(c)) for c in colours]
        true_verts = np.dot(cell, (verts /region_labels.shape).T).T
        data_plot = go.Mesh3d(x=true_verts[:, 0], y=true_verts[:, 1], z=true_verts[:, 2], intensity=colours.astype(int),
                        hoverinfo='all', hovertext=hovertext, flatshading=False, colorscale=cscale, i=faces[:, 0],
                        j=faces[:, 1], k=faces[:, 2], opacity=opacity, name='isosurface',showlegend=True)

    else: #if no region labels are provided color by a flat color

        true_verts = np.dot(cell, (verts /dgrid.shape).T).T
        data_plot = go.Mesh3d(x=true_verts[:, 0], y=true_verts[:, 1], z=true_verts[:, 2],
                        hoverinfo='all', 
                        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2], opacity=opacity, color=color, name='pore surface', showlegend=True)
    
    if stylize:    #some random style choices, purely personal preference. 
        #can be modified after ploting as well. see plotly docs.                
        data_plot.update(               lighting=dict(ambient=0.18,
                                                    diffuse=1,
                                                    fresnel=0.1,
                                                    specular=1,
                                                    roughness=0.05,
                                                    facenormalsepsilon=1e-15,
                                                    vertexnormalsepsilon=1e-15),
                                        lightposition=dict(x=10,
                                                    y=10,
                                                    z=30
                                                        ),showscale=False
                        )

    return data_plot

def plot_cell(ase_atoms, color='purple', name='cell'):
    """Returns a a plotly scatter3D line trace for the 3D unit cell

    :param ase_atoms: ASE atoms object
    :type ase_atoms: ase.atoms
    :param color: line color for the unit cell. Can be any of the plotly colors, refer to plotly docs or see here 
    https://community.plotly.com/t/plotly-colours-list/11730/3  for a list, defaults to 'purple'
    :type color: str, optional
    :param name: Name of the trace, defaults to 'cell'
    :type name: str, optional
    """

    def draw_cell(A):
        import numpy as np
        avec = A[:, 0].T
        bvec = A[:, 1].T
        cvec = A[:, 2].T
        origin = np.array([0, 0, 0])
        cell_lines = []
        cell_lines.append(origin)
        cell_lines.append(avec)
        cell_lines.append([None, None, None])
        cell_lines.append(origin)
        cell_lines.append(bvec)
        cell_lines.append([None, None, None])
        cell_lines.append(origin)
        cell_lines.append(cvec)
        cell_lines.append([None, None, None])
        cell_lines.append(avec)
        cell_lines.append(avec + bvec)
        cell_lines.append([None, None, None])
        cell_lines.append(bvec)
        cell_lines.append(avec + bvec)
        cell_lines.append([None, None, None])
        cell_lines.append(bvec)
        cell_lines.append(cvec + bvec)
        cell_lines.append([None, None, None])
        cell_lines.append(cvec)
        cell_lines.append(cvec + bvec)
        cell_lines.append([None, None, None])
        cell_lines.append(cvec)
        cell_lines.append(cvec + avec)
        cell_lines.append([None, None, None])
        cell_lines.append(avec)
        cell_lines.append(cvec + avec)
        cell_lines.append([None, None, None])
        cell_lines.append(avec + bvec)
        cell_lines.append(cvec + avec + bvec)
        cell_lines.append([None, None, None])
        cell_lines.append(avec + cvec)
        cell_lines.append(cvec + avec + bvec)
        cell_lines.append([None, None, None])
        cell_lines.append(cvec + bvec)
        cell_lines.append(cvec + avec + bvec)
        cell_lines.append([None, None, None])
        return np.array(cell_lines)

    import plotly.graph_objects as go
    from ase.io import read

    A_unit =ase_atoms.get_cell().T
    cell_lines = draw_cell(A_unit)
    data_cell = go.Scatter3d(z=cell_lines[:, 2], y=cell_lines[:, 1], x=cell_lines[:, 0], mode='lines',
                            line=dict(color=color, width=8), name=name)
    return data_cell

def get_default_colors():
    """[summary]
    Uses the ASE's jmol colorscheme as the defual for plotting the atoms.     
    :return: [NOrmalized rgb values of the jmol colorscheme as included within ASE python module. ]
    :rtype: [(N,3) numpy array of floats all within [0,1]]
    """
    from ase.data import colors
    return colors.jmol_colors
def draw_framework(ase_atoms, colorscheme=get_default_colors(), opacity=[1.0, 1.0],dia_scale=1):
    """ Take a grid3D object and plot the atoms and bonds after detecting them using PyMatGen and return
    the Plotly scatter trace to be plotted alongside isosurfaces, windows and diffusion paths.  
    :param t1: [grid3D object that needs to be plotted]
    :type t1: [grid3D object]
    """
    from ase.io import read
    from ase import neighborlist
    # data = read(path_to_cif)
    A_unit = ase_atoms.get_cell().T
    coord = ase_atoms.get_positions()
    names = ase_atoms.get_chemical_symbols()

    def find_covalent_bonds(molecule):
        import itertools
        from pymatgen.core.bonds import CovalentBond
        number_of_sites = len(molecule.sites)
        bonds = []
        for a, b in itertools.combinations(range(number_of_sites), 2):
            if CovalentBond.is_bonded(molecule.sites[a], molecule.sites[b], default_bl=2.0):
                bonds.append([a, b])
        return bonds

    import numpy as np
    from pymatgen.core.structure import Molecule
    print("1. Creating Molecule ...")
    # coord = np.dot(A_unit, (coord+[0.5,0.5,0.5]).T).T
    mol = Molecule(names, coord)
    print("2. Computing bonds ...")
    bonds = find_covalent_bonds(mol)
    # print(bonds)
    # print(np.hstack(bonds))
    edges_master = []
    # bonds_atoms_list = []
    for b in bonds:
        edges_master.append(coord[b[0]])
        edges_master.append(coord[b[1]])
        edges_master.append([None, None, None])

        # * Append the atom id's that are in each bond, None is a spacer.
        # bonds_atoms_list.append(b[0])
        # bonds_atoms_list.append(b[1])
        # bonds_atoms_list.append(None)

    edges_master = np.vstack(edges_master)
    # print(bonds_atoms_list)
    print("3. Creating Plotly figure traces")
    # * Figure out how to plot the framework atoms
    import plotly.graph_objects as go
    # atom_colors           = {'C':'teal', 'H':'yellow', 'O':'red', 'N':'blue ', 'Zn':'red', 'Zr': 'pink'} # this line needs to be updated.
    # scatter_colors        = [atom_colors[t] for t in np.array(names)]

    # * Create a color scale using ASE cpk colorlist.
    from ase.data import colors
    from sklearn.preprocessing import minmax_scale
    cols = np.round(colorscheme * 255, decimals=0).astype(int)
    colorlist = ["rgb(" + ",".join(c.astype(str)) + ")" for c in cols]
    # atomic_numbers        = np.vstack(data.get_atomic_numbers())
    # atom_colors           =   minmax_scale(np.linspace(0,103, 104), feature_range=(0,1))
    atom_colors = np.linspace(0, 1, len(cols))
    cscale = [[atom_colors[i], colorlist[i]] for i in range(len(colorlist))]

    # * Hovertext for each point
    hovertext = ['Atom No: ' + str(i) + ', Atom type: ' + names[i] for i in range(len(names))]

    # * Scale the atom sizes beautifully.
    # atom_dia              = f1.sigma_array*2-f1.sigma
    atom_dia = 2 * neighborlist.natural_cutoffs(ase_atoms)
#   print(atom_dia)
    scatter_size = np.round(minmax_scale(atom_dia, feature_range=(12, 30))*dia_scale)

    bond_atoms_diameter = np.array(atom_dia)[np.hstack(bonds)]
    bond_atoms_Z = ase_atoms.get_atomic_numbers()[np.hstack(bonds)] - 1  # The indices in the color scale start at 0, not 1.
    bond_scatter_size = np.round(minmax_scale(bond_atoms_diameter, feature_range=(12, 35)))

    data_scatter = go.Scatter3d(x=coord[:, 0], y=coord[:, 1], z=coord[:, 2], mode='markers', hovertext=hovertext,
                                marker=dict(opacity=opacity[0], size=scatter_size,
                                            color=atom_colors[ase_atoms.get_atomic_numbers() - 1], colorscale=cscale,
                                            sizemode='diameter'), showlegend=True, name='framework atoms')
    data_edges = go.Scatter3d(x=edges_master[:, 0], y=edges_master[:, 1], z=edges_master[:, 2],
                            line=dict(width=15, color='grey'), mode='lines', opacity=opacity[1],
                            name='framework bonds', showlegend=True)
    # data_framework          = go.Scatter3d(x=edges_master[:,0] , y=edges_master[:,1], z=edges_master[:,2] , line=dict(width=15, color='grey'), mode='markers+lines', hovertext=hovertext, marker=dict(opacity=opacity[0], size=bond_scatter_size, color=bond_atoms_Z, colorscale=cscale), opacity=opacity[1], name='Framework')
    return data_scatter, data_edges

def scatter3d(X, name='trace name here', mode='markers', marker_color = 'blue', marker_size=12, marker_line_color ='white', marker_symbol = 'circle', marker_line_width=1, line_color='black', line_width=2, hovertext=None, text=None, **kwargs):        
    
        """ This is a short hand for the plotly go.Scatter3d function which returs a 3D scatter trace. Just keeping the basic options for a quick plot of scatter traces.
        :type X: Arrya of 3D coordinates (N,3)
        :param X: numpy.array
    
        :type name: name for the trace
        :param name: str
    
        :type mode: One of 'markers', 'lines' or 'markers+lines'
        :param mode: str
    
        :type marker_color: marker color, check plotly docs for acceptable colors.this list can be found here:https://community.plotly.com/t/plotly-colours-list/11730/3
        :param marker_color: str
    
        :type marker_size: size of the marker, defaults to 12
        :param marker_size: int
    
        :type marker_line_color: line color of the marker, defaults to white
        :param marker_line_color: str
    
        :type marker_symbol: marker symbol, defaults to circle. check plotly docs for all the options
        :param marker_symbol: str
    
        :type marker_line_width: line width of the marker outline. defaults to 1
        :param marker_line_width: float
    
        :type line_color: color of the conncting lines in the scatter plot  defaults to black
        :param line_color: str
    
        :type line_width: width of the line  defaults to 1.0
        :param line_width: float
    
        :type hovertext: text to display on hover
        :param hovertext: string or a list of strings of the same length as the data array
    
        :raises:
    
        :rtype: a scatter trace of the type plotly scatter3d, which can then be added to a plotly go.Figure or go.FigureWidget
        """
        import plotly.graph_objects as go 
        return go.Scatter3d(x=X[:,0], y=X[:,1], z=X[:,2], mode=mode, marker =dict(color=marker_color, size=marker_size, symbol=marker_symbol, line=dict(color=marker_line_color, width=marker_line_width)), name=name, line=dict(color=line_color, width=line_width), text=None, hovertext=hovertext, **kwargs)

def plot_pore_graph(G):


    import networkx as nx 
    import plotly.graph_objects as go 
    import numpy as np
    coords = [n[1]['node_coords'] for n in G.nodes(data=True)]
    coords = np.vstack(coords)
    # print(coords)
    edges_list = [] 
    for edge in list(G.edges()):
        edges_list.append(G.nodes(data=True)[edge[0]]['node_coords'])
        edges_list.append(G.nodes(data=True)[edge[1]]['node_coords'])
        edges_list.append([None, None, None])

    nodes_trace = scatter3d(np.array(coords), name='nodes', showlegend=True, legendgroup='pore_graph', legendgrouptitle_text='pore graph')
    edges_trace = scatter3d(np.array(edges_list), mode='lines', name='edges', line_width=8, showlegend=True, legendgroup='pore_graph', legendgrouptitle_text='pore graph')
    
    return nodes_trace, edges_trace

    

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    # Do something if this file is invoked on its own
    print("You probably shouldn't call this by itself")



def draw_framework2(ase_atoms, bonds=None, bond_orders=None, show_aromatics=True, show_atom_labels=False, bond_data=None, bond_color=None, bond_opacity=0.9, colorscheme=get_default_colors(), atoms_opacity=0.9, dia_scale=1):
        """
        It takes an ASE atoms object and returns a list of plotly traces that can be used to plot the
        framework
        
        :param ase_atoms: ASE atoms object
        :param bonds: list of bonds in the form of [[0,1],[1,2],...]
        :param bond_orders: A list of bond orders for each bond in the bonds list
        :param show_rings: If True, the function will detect and plot aromatic rings, defaults to True
        (optional)
        :param bond_data: A dictionary of bond data. The keys are the names of the bond data, and the values
        are the bond data. The bond data should be a list of the same length as the number of bonds in the
        structure
        :param colorscheme: A dictionary of colors for each element. The default is the ASE CPK color scheme
        :param atoms_opacity: Opacity of the atoms
        :param dia_scale: This is a scaling factor for the atom sizes, defaults to 1 (optional)
        :param default_ring_sizes: list of integers, default = [5,6,7,8]
        :return: A list of plotly.graph_objects.Scatter3d objects.
        """

        # from ase.io import read
        from ase import neighborlist
        from .analyze import update_with_covalent_bonds, find_aromatic_cycles_in_sgraph
        # from ase.data import covalent_radii 
        # data = read(path_to_cif)
        # A_unit = data.get_cell().T
        coord = ase_atoms.get_positions()
        names = ase_atoms.get_chemical_symbols()



        import numpy as np
        from pymatgen.core import Molecule
        # print("1. Creating Molecule ...")
        # mol = Molecule(names, coord)


        # print(bonds_atoms_list)
        print("Creating Plotly figure traces")
        # * Figure out how to plot the framework atoms

        import plotly.graph_objects as go

        # * Create a color scale using ASE cpk colorlist.
        from ase.data import colors
        from sklearn.preprocessing import minmax_scale
        cols = np.round(colorscheme * 255, decimals=0).astype(int)
        colorlist = ["rgb(" + ",".join(c.astype(str)) + ")" for c in cols]
        # atomic_numbers        = np.vstack(data.get_atomic_numbers())
        # atom_colors           =   minmax_scale(np.linspace(0,103, 104), feature_range=(0,1))
        atom_colors = np.linspace(0, 1, len(cols))
        cscale = [[atom_colors[i], colorlist[i]] for i in range(len(colorlist))]

        # * Hovertext for each point
        hovertext = ['Atom No: ' + str(i) + ', Atom type: ' + names[i] for i in range(len(names))]

        # * Scale the atom sizes beautifully.
        # atom_dia              = f1.sigma_array*2-f1.sigma
        atom_dia = 2 * neighborlist.natural_cutoffs(ase_atoms)

        # # * Create a color scale using ASE cpk colorlist.
        # from ase.data import colors
        # from sklearn.preprocessing import minmax_scale
        # colors_in_mof = colorscheme[ase_atoms.get_atomic_numbers()]
        # cols = np.round(colors_in_mof * 255, decimals=0).astype(np.int)
        # colorlist = ["rgb(" + ",".join(c.astype(np.str)) + ")" for c in cols]
        
        # # * Hovertext for each point
        # hovertext = ['Atom No: ' + str(i) + ', Atom type: ' + names[i] for i in range(len(names))]
        
        # # * Scale the atom sizes beautifully.
        # atom_dia = 2 * neighborlist.natural_cutoffs(ase_atoms)
        # print(atom_dia)
        scatter_size = np.round(minmax_scale(atom_dia, feature_range=(12, 30))*dia_scale)


        traces=[]
        print("1. Creating atoms trace\n")

        trace_scatter = go.Scatter3d(x=coord[:, 0], y=coord[:, 1], z=coord[:, 2], mode='markers', hovertext=hovertext,
                                    marker=dict(opacity=atoms_opacity, size=scatter_size,
                                                color=atom_colors[ase_atoms.get_atomic_numbers()-1],
                                                sizemode='diameter',colorscale=cscale), name='framework atoms', text=names, textfont=dict(color=colorlist, size=20),surfaceaxis=-1)
        traces.append(trace_scatter)
        
        if bonds is None:
            ase_atoms_updated = update_with_covalent_bonds(ase_atoms)
            bonds = ase_atoms_updated.info['covalent_bonds']
            if bond_orders is None:
                bond_orders = ase_atoms_updated.info['covalent_bond_orders'] 
            else:
                raise Exception("You can't provide bond orders without providing a list of bonds!")
        else:
            ase_atoms_updated = ase_atoms

        print("2. Creating bonds' trace\n")
        edges_master = []   
        for b in bonds:
            edges_master.append(coord[b[0]])
            edges_master.append(coord[b[1]])
            edges_master.append([None, None, None])
        
        if bond_color is None:
            bond_color = np.repeat(bond_orders,3)
        elif type(bond_color)==str:
            bond_color=bond_color
        else:
            bond_color = 'rgba(0,0,0,0.35)'
        edges_master = np.vstack(edges_master)
        print(bond_color)
        try:
            trace_edges = go.Scatter3d(x=edges_master[:, 0], y=edges_master[:, 1], z=edges_master[:, 2],
                                line=dict(width=15, color=bond_color, colorscale='aggrnyl', colorbar=dict(thickness=20, len=0.6, title='bond order')), opacity=bond_opacity, mode='lines',
                                name='framework bonds')
        except:
            trace_edges = go.Scatter3d(x=edges_master[:, 0], y=edges_master[:, 1], z=edges_master[:, 2],
                                line=dict(width=15, color='rgba(150,150,150,0.35)', colorscale='aggrnyl' ),opacity=bond_opacity, mode='lines',
                                name='framework bonds')
        traces.append(trace_edges)


        if show_aromatics:
            print('\t2-a. Detecting and plotting aromatic rings\n')
            # show aromatic rings 
            # if bond order data is provided

        # Detect and show requested rings
            # bond_types = nx.get_edge_attributes(G, 'bond_type')
            # Detect aromatic rings
            G = sgraph_from_atoms(ase_atoms_updated)
            isomorphic_aromatic_cycles = find_aromatic_cycles_in_sgraph(G)                 
            # Highlight and circle all aromatic rings

            ring_traces = plot_aromatic_rings(sgraph=G, isomorphic_aromatic_cycles=isomorphic_aromatic_cycles)
            traces = traces + list(ring_traces)
            # circle_data=[]
            # from numpy.core.numeric import cross
            # theta = np.linspace(0, 2*np.pi, 100)
            # xc = 0.5*np.cos(theta)
            # yc = 0.5*np.sin(theta) 
            # zc = np.ones(np.shape(xc))
            # circle_coord = np.array([xc, yc, zc]).T

            # from scipy.spatial.transform import Rotation as R

            # # ring circles
            # for ring in isomorphic_aromatic_cycles:
            #     ring_coord = coord[ring]
            #     ring_center = np.mean(ring_coord, axis=0)
            #     ring_normal = np.cross(ring_coord[3]-ring_coord[2],ring_coord[2]-ring_coord[1])
            #     ring_normal = ring_normal/np.linalg.norm(ring_normal)
            #     r = R.align_vectors(ring_normal.reshape(-1,3),np.array([0,0,1]).reshape(-1,3))
            #     new_circle = r[0].apply(circle_coord)
            #     new_circle= new_circle + ring_center-np.mean(new_circle, axis=0) 
            #     circle_data.append(new_circle)
            #     circle_data.append([None, None, None])
            # circle_data = np.vstack(circle_data)
            # data_circle = go.Scatter3d(x=circle_data[:,0], y=circle_data[:,1], z=circle_data[:,2], mode='lines', line=dict(color='white',width=2), showlegend=False, name='Aromatic Cloud', legendgroup='aromatic circles')
            
            # traces.append(data_circle)

            # # ring highlight
            # ring_centers = []
            # for ring in isomorphic_aromatic_cycles:
            #     ring_coord = coord[ring]
            #     ring_centers.append(np.mean(ring_coord, axis=0))
            #     ring_highlight_trace = go.Mesh3d(x=ring_coord[:,0],y=ring_coord[:,1],z=ring_coord[:,2], color='lightskyblue', opacity=0.1, alphahull=0, name='Benzene ring highlight')
            #     traces.append(ring_highlight_trace)


        if bool(bond_data):
            print('4. Creating traces for additional bond data provided')
            for key in bond_data.keys():

                true_edges = edges_master[edges_master != None ].reshape(-2,3)
                true_edges = 0.5*(true_edges[0::2]+true_edges[1::2])
                trace_bond_data = go.Scatter3d(x=true_edges[:, 0], y=true_edges[:, 1], z=true_edges[:, 2], mode='text',  textfont=dict(color='white', size=9), text=list(map(str, bond_data[key])),
                                        name=key )
                traces.append(trace_bond_data)
        if show_atom_labels:
            trace_atom_labels = go.Scatter3d(x=coord[:, 0], y=coord[:, 1], z=coord[:, 2], mode='text',  textfont=dict(color=colorlist, size=20), text=names,
                                    name='atom labels' )
            traces.append(trace_atom_labels)
        return traces


def plot_aromatic_rings(sgraph, isomorphic_aromatic_cycles=None, highlight_color='lightskyblue'):
    """
    It takes in a graph and a list of aromatic cycles, and returns a tuple of plotly traces that can be
    used to plot the aromatic cycles
    
    :param sgraph: the sgraph object
    :param isomorphic_aromatic_cycles: a list of lists of node indices that are aromatic rings. If you
    don't have this, you can set it to None and the function will try to find it for you
    :param highlight_color: The color of the highlighted aromatic rings, defaults to lightskyblue
    (optional)
    :return: A tuple of the data_circle and ring_highlight_traces
    """

    from .analyze import find_aromatic_cycles_in_sgraph
    import plotly.graph_objects as go
    import numpy as np 

    if not (isomorphic_aromatic_cycles is None):
        try:
            isomorphic_aromatic_cycles = find_aromatic_cycles_in_sgraph(sgraph)
        except:
            Exception('Cannot find aromatic cycles from sgraph! Please provide as input.')
    circle_data=[]
    from numpy.core.numeric import cross
    theta = np.linspace(0, 2*np.pi, 100)
    xc = 0.5*np.cos(theta)
    yc = 0.5*np.sin(theta) 
    zc = np.ones(np.shape(xc))
    circle_coord = np.array([xc, yc, zc]).T

    from scipy.spatial.transform import Rotation as R
    import networkx as nx 
    # coord  = nx.get_node_attributes(sgraph, "position")# ring circles
    try :
        coord = np.vstack(list(nx.get_node_attributes(sgraph, "position").values()))
    except:
        raise Exception("Cannot find position data in your sgraph")
    for ring in isomorphic_aromatic_cycles:
        ring_coord = coord[ring]
        ring_center = np.mean(ring_coord, axis=0)
        ring_normal = np.cross(ring_coord[3]-ring_coord[2],ring_coord[2]-ring_coord[1])
        ring_normal = ring_normal/np.linalg.norm(ring_normal)
        r = R.align_vectors(ring_normal.reshape(-1,3),np.array([0,0,1]).reshape(-1,3))
        new_circle = r[0].apply(circle_coord)
        new_circle= new_circle + ring_center-np.mean(new_circle, axis=0) 
        circle_data.append(new_circle)
        circle_data.append([None, None, None])
    circle_data = np.vstack(circle_data)
    data_circle = go.Scatter3d(x=circle_data[:,0], y=circle_data[:,1], z=circle_data[:,2], mode='lines', line=dict(color='white',width=2), showlegend=False, name='Aromatic Cloud', legendgroup='aromatic circles')
    

    # ring highlight
    ring_centers = []
    ring_highlight_traces = []
    for ring in isomorphic_aromatic_cycles:
        ring_coord = coord[ring]
        ring_centers.append(np.mean(ring_coord, axis=0))
        ring_highlight_trace = go.Mesh3d(x=ring_coord[:,0],y=ring_coord[:,1],z=ring_coord[:,2], color=highlight_color, opacity=0.1, alphahull=0, name='Aromatic ring highlight')
        ring_highlight_traces.append(ring_highlight_trace)

    return tuple([data_circle] + ring_highlight_traces)

#* Functions for ease
def trimesh_from_pyvista(pv_mesh):
    """
    `trimesh_from_pyvista(pv_mesh)` takes a PyVista mesh and returns a Trimesh mesh
    
    :param pv_mesh: a pyvista mesh object
    :return: A trimesh object.
    """
    import trimesh
    faces_array = pv_mesh.faces.reshape((pv_mesh.n_faces,4))[:,1:]
    tmesh = trimesh.Trimesh(pv_mesh.points, faces_array)
    return tmesh

def spheres_for_atoms(ase_atoms,radii=None, scale=0.8):
    """
    It takes an ASE atoms object and returns a list of trimesh objects, each of which is a sphere with
    radius and center corresponding to the atom's covalent radius and position.
    
    :param ase_atoms: an ase atoms object
    :param radii: the radii of the spheres. If None, the covalent radii from ase.data.covalent_radii are
    used
    :param scale: the scale of the sphere radius
    :return: A list of trimesh objects.
    """
    import numpy as np
    import pyvista as pv 
    from ase.data import covalent_radii 
    natoms = ase_atoms.get_global_number_of_atoms()
    if radii is None:
        radii = covalent_radii[ase_atoms.get_atomic_numbers()]
    positions =ase_atoms.get_positions()
    spheres = [pv.Sphere(radius=radii[i]*scale, center=positions[i],theta_resolution=20, phi_resolution=20) for i in range(natoms)]
    spheres_tmesh= [trimesh_from_pyvista(s) for s in spheres]
    return spheres_tmesh

def cylinders_for_bonds(bonds,ase_atoms, radius=0.25):
    """
    It takes a list of bonds (each bond is a list of two indices) and an ASE atoms object, and returns a
    list of trimesh objects, each of which is a cylinder centered on the bond with radius 0.25 and
    height equal to half the bond length
    
    :param bonds: a list of pairs of atom indices
    :param ase_atoms: an ASE atoms object
    :return: A list of trimesh objects.
    """
    import numpy as np
    import pyvista as pv
    positions = ase_atoms.get_positions()
    bond_coords = np.array([positions[b] for b in bonds])
    directions  = bond_coords[:,1] - bond_coords[:,0]
    centers = (bond_coords[:,0] + bond_coords[:,1])/2
    heights = np.linalg.norm(directions, axis=1)
    cylinders = [pv.Cylinder(center=c, direction=d, resolution=25, radius=radius, height=h).triangulate() for c,d,h in zip(centers, directions, heights)]
    cylinders_tmesh = [trimesh_from_pyvista(c) for c in cylinders]
    return cylinders_tmesh

    """
    > This function takes an ASE atoms object and returns a trimesh object that contains spheres and
    cylinders for the atoms and bonds, respectively. 
    
    The function is a bit more complicated than it needs to be because it has to handle the case where
    the bonds are not provided. In that case, it uses the `mofography` package to calculate the bonds. 
    
    The function also has a few optional arguments that allow you to customize the appearance of the
    spheres and cylinders. 
    
    Let's try it out on a simple example: 
    
    # Example
    
    :param ase_atoms: an ase.Atoms object
    :param bonds: a list of pairs of atom indices. If None, then the bonds are calculated from the atoms
    :param colormap: a dictionary of colors for each atomic number
    :param bond_rgba_array: the color of the bonds
    :param kwargs_spheres: 
    :param kwargs_cylinders: radius of the cylinder
    :return: A trimesh object
    """
def mesh_for_atoms(ase_atoms, bonds=None, colormap=None, bond_rgba_array=[128,128,128,255], kwargs_spheres=dict(scale=0.8, radii=None), kwargs_cylinders=dict(radius=0.25)):
    import trimesh 
    import numpy as np 
    if colormap is None:
        from ase.data.colors import cpk_colors 
        colormap = cpk_colors
    if bonds is None:
        import mofography as mgr
        ase_atoms_updated= mgr.update_with_covalent_bonds(ase_atoms)
        ase_atoms = ase_atoms_updated
        import mofography as mgr
        atoms_graph = mgr.sgraph_from_atoms(ase_atoms_updated)
        bonds = np.vstack(list(atoms_graph.edges()))
    
    spheres = spheres_for_atoms(ase_atoms, **kwargs_spheres)
    atomic_numbers =ase_atoms.get_atomic_numbers()
    for i,s in enumerate(spheres):
        s.visual.vertex_colors = cpk_colors[atomic_numbers[i]]
    cylinders = cylinders_for_bonds(bonds,ase_atoms, **kwargs_cylinders)
    for i,c in enumerate(cylinders):
        c.visual.vertex_colors = bond_rgba_array #cpk_colors[atomic_numbers[i]]
    combined = trimesh.util.concatenate([spheres+cylinders])

    return combined


def plot_framework_mesh(combined, stylize=False):
    """
    It takes a combined framework object and returns a plotly mesh trace
    
    :param combined: the combined mesh object
    :param stylize: if True, the mesh will be rendered with a nice lighting effect, defaults to False
    (optional)
    :return: A mesh object
    """

    import plotly.graph_objects as go 
    verts = combined.vertices
    faces= combined.faces
    colors = ['rgba('+ ','.join(c.astype(str))+')' for c in combined.visual.vertex_colors]
    framework_mesh_trace= go.Mesh3d(x=verts[:,0], y=verts[:,1], z=verts[:,2], flatshading=False, i=faces[:,0], j=faces[:,1], k=faces[:,2], vertexcolor=colors, showlegend=True, showscale=False, name='framework')
    if stylize:
        framework_mesh_trace.update(               lighting=dict(ambient=0.18,
                                            diffuse=1,
                                            fresnel=0.1,
                                            specular=1,
                                            roughness=0.05,
                                            facenormalsepsilon=1e-15,
                                            vertexnormalsepsilon=1e-15),
                                lightposition=dict(x=10,
                                            y=10,
                                            z=30
                                                ),showscale=False
                )
    return framework_mesh_trace


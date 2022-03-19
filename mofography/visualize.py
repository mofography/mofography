"""
 
 Visualization subroutines for mofography
 ██╗   ██╗██╗███████╗██╗   ██╗ █████╗ ██╗     ██╗███████╗███████╗
 ██║   ██║██║██╔════╝██║   ██║██╔══██╗██║     ██║╚══███╔╝██╔════╝
 ██║   ██║██║███████╗██║   ██║███████║██║     ██║  ███╔╝ █████╗  
 ╚██╗ ██╔╝██║╚════██║██║   ██║██╔══██║██║     ██║ ███╔╝  ██╔══╝  
  ╚████╔╝ ██║███████║╚██████╔╝██║  ██║███████╗██║███████╗███████╗
   ╚═══╝  ╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝╚══════╝╚══════╝
                                                                 
 
"""
def plot_surface_dgrid(dgrid, isovalue,  ase_atoms,region_labels=None, select='all', opacity=1.0, cscale=None, color='blue',  stepsize=1, stylize=False):
    """Creates a plolty mesh3D trace of an isosurface of a given value for a given 3D distance/energy grid. The surface can then
    be colored by the regions or any other property

    :param dgrid: 3D numpy array of grid values
    :type dgrid: numpy.array
    :param isovalue: value at which the isosurface is calculated using marching cubes.
    :type isovalue: float
    :param region_labels: 3D nmpy array of region labels or any other coloring function
    :type region_labels: numpy.array
    :param ase_atoms: ASE atoms object containing crystallographic information.
    :type ase_atoms: ase.atoms
    :param select: regions to select, useful for showing just one or more regions, defaults to 'all'
    :type select: int or same type as the coloring function, optional
    :param opacity: opacity of the plotted isosurface, defaults to 1.0
    :type opacity: float, optional
    :param cscale: colorscale used to color the plotting surface, defaults to a custom random colorscale everytime you plot
    :type cscale: can be one of plotly.express.colors.sequential/continuous/discrete or a custom color scale as specified by plotly, optional, defaults to .
    if colorscale is none, use the color argument to plot in a single color. Color names can be found in plotly docs or see this discussion
    https://community.plotly.com/t/plotly-colours-list/11730/3
    :param color: color of the isosurface used only when colorscale is not specified, defaults to 'blue'
    :type color: str, optional
    :param stepsize: whether to use all the grid points for the isosurface or to choose every nth point. Can be used to save visualization load on
    less capable machines.  , defaults to 1
    :type stepsize: int, optional
    :param stylize: secret set of styling options to make your surface look cool, purely personal preference. Refer to plotly
    documentation for detailed options, defaults to False
    :type stylize: bool, optional
    """        
    import numpy as np
    import copy
    import plotly.graph_objects as go
    import ase 
    dgrid_copy = copy.deepcopy(dgrid) # deepcopy of the grid
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()) #get cell matrix
    # Exclude the regions that are not selected
    if select != 'all':
        dgrid_copy[np.isin(region_labels, select, assume_unique=False, invert=True)] = -10

    # Compute the isosurface using Marching Cubes Lewiner
    from skimage.measure import marching_cubes as mc
    verts, faces, normals, values = mc(dgrid_copy, isovalue, step_size=stepsize)  

    if np.any(region_labels) != None: # if region labels are provided

        #  Make a colorscale if no colorscale is provided
        if cscale == None:
            cvals = np.linspace(0, 1, (region_labels.max()).astype(int))
            cscale = [[cvals[i],
                    "rgb(" + str(np.random.randint(0, 256)) + ',' + str(np.random.randint(0, 256)) + ',' + str(
                        np.random.randint(0, 256)) + ')'] for i in range((region_labels.max()).astype(int))]

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
                        j=faces[:, 1], k=faces[:, 2], opacity=opacity, name='Isosurface',showlegend=True)

    else: #if no region labels are provided color by a flat color

        true_verts = np.dot(cell, (verts /dgrid.shape).T).T
        data_plot = go.Mesh3d(x=true_verts[:, 0], y=true_verts[:, 1], z=true_verts[:, 2],
                        hoverinfo='all', 
                        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2], opacity=opacity, color=color, name='Isosurface', showlegend=True)
    
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

    A_unit =ase_atoms.get_cell()
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
def draw_framework(data, colorscheme=get_default_colors(), opacity=[1.0, 1.0],dia_scale=1):
      """ Take a grid3D object and plot the atoms and bonds after detecting them using PyMatGen and return
      the Plotly scatter trace to be plotted alongside isosurfaces, windows and diffusion paths.  
      :param t1: [grid3D object that needs to be plotted]
      :type t1: [grid3D object]
      """
      from ase.io import read
      from ase import neighborlist
      # data = read(path_to_cif)
      A_unit = data.get_cell().T
      coord = data.get_positions()
      names = data.get_chemical_symbols()

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
      atom_dia = 2 * neighborlist.natural_cutoffs(data)
    #   print(atom_dia)
      scatter_size = np.round(minmax_scale(atom_dia, feature_range=(12, 30))*dia_scale)

      bond_atoms_diameter = np.array(atom_dia)[np.hstack(bonds)]
      bond_atoms_Z = data.get_atomic_numbers()[np.hstack(bonds)] - 1  # The indices in the color scale start at 0, not 1.
      bond_scatter_size = np.round(minmax_scale(bond_atoms_diameter, feature_range=(12, 35)))

      data_scatter = go.Scatter3d(x=coord[:, 0], y=coord[:, 1], z=coord[:, 2], mode='markers', hovertext=hovertext,
                                  marker=dict(opacity=opacity[0], size=scatter_size,
                                              color=atom_colors[data.get_atomic_numbers() - 1], colorscale=cscale,
                                              sizemode='diameter'), name='framework_atoms')
      data_edges = go.Scatter3d(x=edges_master[:, 0], y=edges_master[:, 1], z=edges_master[:, 2],
                                line=dict(width=15, color='grey'), mode='lines', opacity=opacity[1],
                                name='framework_bonds')
      # data_framework          = go.Scatter3d(x=edges_master[:,0] , y=edges_master[:,1], z=edges_master[:,2] , line=dict(width=15, color='grey'), mode='markers+lines', hovertext=hovertext, marker=dict(opacity=opacity[0], size=bond_scatter_size, color=bond_atoms_Z, colorscale=cscale), opacity=opacity[1], name='Framework')
      return data_scatter, data_edges

        
if __name__ == "__main__":
    # Do something if this file is invoked on its own
    # Do something if this file is invoked on its own
    print("You probably shouldn't call this by itself")

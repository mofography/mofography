"""
 Compute functions for mofography
  ██████╗ ██████╗ ███╗   ███╗██████╗ ██╗   ██╗████████╗███████╗
 ██╔════╝██╔═══██╗████╗ ████║██╔══██╗██║   ██║╚══██╔══╝██╔════╝
 ██║     ██║   ██║██╔████╔██║██████╔╝██║   ██║   ██║   █████╗  
 ██║     ██║   ██║██║╚██╔╝██║██╔═══╝ ██║   ██║   ██║   ██╔══╝  
 ╚██████╗╚██████╔╝██║ ╚═╝ ██║██║     ╚██████╔╝   ██║   ███████╗
  ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝      ╚═════╝    ╚═╝   ╚══════╝
                                                               
 
"""

def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format).

    Replace this function and doc string for your own project.

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from.

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution.
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


def gpd_aabb(grid_points, ase_atoms, n_neighbors=10, probe_radius=0.0 ):
    """This subroutine computes the shortest distance to the framework surface for any grid point using the Axis Aligned
    Bounding Box method  as implemented by the Glotzer group from Univ. of Michigan.

    :param grid_points: (N,3) array of the cartesian coordinates of the grid points
    :type grid_points: numpy array of floats
    :param ase_atoms: ase atoms object containing the framework information
    :type ase_atoms: ase.atoms obj
    :param n_neighbors: number of neighbors of each grid point checked for the shortest dist to surface, defaults to 10
    :type n_neighbors: int, optional
    :param probe_radius: radius of the distance grid calculation probe, defaults to 0.0
    :type probe_radius: float, optional
    :return: (N,1) array of the shortest distance corresponding to the N grid points.
    :rtype: numpy array of floats
    """    
    from freud.locality import AABBQuery
    from freud.box import Box
    import ase 
    from ase.data import vdw_radii

    import numpy as np


    framework_atom_positions = ase_atoms.get_positions()
    cell = ase.geometry.complete_cell(ase_atoms.get_cell())
    radii = vdw_radii[ase_atoms.get_atomic_numbers()]

    box = Box.from_matrix(cell, dimensions =3) # Define the box for Freud
    ab  = AABBQuery(box=box, points=framework_atom_positions) # Compute the AABB Tree
    # Perform the nearest neighbor query with n-neighbors (10 by default)
    ab_r = ab.query(grid_points, query_args=dict(mode='nearest', num_neighbors=n_neighbors, exclude_ii=True))
    ab_r_nl = ab_r.toNeighborList(sort_by_distance=True) #convert query to neighborlist
    nls = ab_r_nl.point_indices.reshape(-1,n_neighbors)  # get the indices of the neighbors for each grid pt. to lookup the radii later
    dists =ab_r_nl.distances.reshape(-1,n_neighbors) # get the distances to the neighbors of each grid pt.

    # return the minimum (shortest distance to the surface) in each row
    return np.min(dists - [radii[nl] for nl in nls],axis=1)

def dask_grid_over_atoms(ase_atoms, spacing=0.1, chunksize=50000):
    """Generate a dask array of grid points over an ase atoms object with a specified grid spacing

    :param ase_atoms: ASE atoms object
    :type ase_atoms: ase.atoms obj
    :param spacing: grid fineness in angstroms, defaults to 0.1
    :type spacing: float, optional
    :param chunksize: number of grid points in each block of the dask array, defaults to 50000
    :type chunksize: int, optional
    :return: dask array of size (N_grid_points, 3) with each block of size 'chunksize'. 
    :rtype: dask.array
    """    
    import ase
    import dask.array as da
    import numpy as np
    #number of grid points in each direction 
    [nx, ny, nz] = (ase_atoms.get_cell_lengths_and_angles()[0:3] / spacing).astype(int) + 1
    gpoints = (da.stack(da.meshgrid(np.linspace(0, 1, nx), np.linspace(0,1, ny), np.linspace(0,1, nz), indexing='ij'),-1).reshape(-1, 3)).rechunk(chunksize,3)
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()) # cell matrix
    return da.dot(cell, gpoints.T).T, (nx,ny,nz) #return the actual coordinates

def grid_from_atoms(ase_atoms, spacing=0.25, chunksize=50000, n_neighbors=10, probe_radius=0):
    """Compute distance grid from ase atoms objects

    :param ase_atoms: ASE atoms object
    :type ase_atoms: ase.atoms
    :param spacing: fineness of the grid in angstrom, defaults to 0.25
    :type spacing: float, optional
    :param chunksize: number of grid points per block of the dask grid. Usually a tradeoff between fewer blocks 
    and data transfer limitation from large blocksizes, defaults to 50000
    :type chunksize: int, optional
    :param n_neighbors: number of neighbors to check for distance, fewer the faster, but risk getting wrong distance 
    since the neighbor sorting by default uses the center to center distance and does not account for the atomic radii, defaults to 10
    :type n_neighbors: int, optional
    :param probe_radius: size of the probe for computing distace to framework surface, defaults to 0.0 corresponding to a point probe
    :type probe_radius: float, optional
    :return: Grid of shape (nx, ny, nz), wher nx, ny, nz are the number of grid points corresponding to the spacing specified earlier
    :rtype: 3D numpy array
    """    
    import numpy as np
    import dask.array as da
    gpoints, shape = dask_grid_over_atoms(ase_atoms, spacing=spacing, chunksize=chunksize) #mesh grid
    dgrid = []
    for block in gpoints.blocks:
        dgrid.append(gpd_aabb(block.compute(), ase_atoms, n_neighbors=n_neighbors,probe_radius=probe_radius )) #compute grid
    return da.hstack(dgrid).rechunk(chunks=chunksize).reshape(shape) #return in correct shape



def compute_dgrid_gpu(data, spacing=0.2, chunk_size=5000):
    """Code to compute the distance grid using the GPU

    :param data: ASE style atoms object
    :type data: ase.atoms
    :param spacing: grid spacing in A, defaults to 0.2
    :type spacing: float, optional
    :param chunk_size: block size of the grid points array, number of grid points to be prcessed by each core of the CPU while talking to the GPU, defaults to 5000
    :type chunk_size: int, optional
    :return: cupy 3D array with the distance values
    :rtype: cupy.array
    """
    from ase.io import read
    import numpy as np
    import dask.array as da
    import cupy as cp
    import pandas as pd
    from tqdm import tqdm
    import ase

    # data    = read(path_to_cif)
    fpoints = cp.asarray(data.get_scaled_positions() - 0.5)
    cell = cp.asarray(ase.geometry.complete_cell(data.get_cell()))
    from ase.data import vdw_radii
    radius_series = pd.Series({'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96,
                            'B': 0.84, 'C': 0.73, 'N': 0.71, 'O': 0.66,
                            'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41,
                            'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
                            'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76,
                            'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39,
                            'Mn': 1.50, 'Fe': 1.42, 'Co': 1.38, 'Ni': 1.24,
                            'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.20,
                            'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16,
                            'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75,
                            'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46,
                            'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
                            'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38,
                            'I': 1.39, 'Xe': 1.40, 'Cs': 2.44, 'Ba': 2.15,
                            'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01,
                            'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96,
                            'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89,
                            'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75,
                            'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44,
                            'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32,
                            'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40,
                            'At': 1.50, 'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21,
                            'Ac': 2.15, 'Th': 2.06, 'Pa': 2.00, 'U': 1.96,
                            'Np': 1.90, 'Pu': 1.87, 'Am': 1.80, 'Cm': 1.69})
    # radii = cp.asarray(radius_series[data.get_chemical_symbols()].values).reshape(-1, 1)
    radii = cp.asarray(vdw_radii[data.get_atomic_numbers()]).reshape(-1, 1)
    [nx, ny, nz] = (data.get_cell_lengths_and_angles()[0:3] / spacing).astype(int) + 1
    #     print([nx,ny,nz])
    gpoints = da.stack(
        da.meshgrid(np.linspace(-0.5, 0.5, nx), np.linspace(-0.5, 0.5, ny), np.linspace(-0.5, 0.5, nz), indexing='ij'),
        -1).reshape(-1, 3).rechunk(chunk_size, 3)

    #     def gpd(points, fpoints=fpoints, cell=cell, radii=radii):
    #         points = cp.asarray(points)
    #         return cp.min(cp.linalg.norm(cp.dot(cell, (cp.expand_dims(points, axis=1)-cp.expand_dims(fpoints, axis=0)-cp.around(cp.expand_dims(points, axis=1)-cp.expand_dims(fpoints, axis=0))).reshape(-1,3).T).T, axis=1).reshape(fpoints.shape[0],-1)-radii, axis=0)

    def gpd(points, fpoints=fpoints, cell=cell, radii=radii):
        points = cp.asarray(points)

        diff = cp.expand_dims(points, axis=1) - cp.expand_dims(fpoints, axis=0)
        diff = (diff - cp.around(diff)).reshape(-1, 3)
        diff = cp.dot(cell, diff.T).T
        diff = cp.linalg.norm(diff, axis=1).reshape(-1, fpoints.shape[0]) - radii.T
        return cp.min(diff, axis=1)

    gpoints_da = gpoints.map_blocks(gpd, chunks=(chunk_size, 1))

    distance_grid = []
    from rich.progress import track
    for block in track(gpoints_da.blocks, total=gpoints_da.npartitions):
        # print("Working on block {} \n".format(i))
        distance_grid.append(block.compute())

    #  Free the gpu memeory before returning the output
    del data, cell, fpoints, radius_series, radii
    cp._default_memory_pool.free_all_blocks()
    cp._default_pinned_memory_pool.free_all_blocks()
    return np.hstack(distance_grid).reshape(nx, ny, nz)



def regions_from_grid(dgrid, peak_min=2.0, mask_thickness=1):
    """Computes the region segmentaion from the distance grid

    :param dgrid: 3D distance grid
    :type dgrid: numpy.array
    :param peak_min:minimum height of the peak in angstrom, defaults to 2.0
    :type peak_min: float, optional, defaults to 2.0 angstroms
    :param mask_thickness:thickness of the framework cover. Points within this distance from the framework surface are not segmented, defaults to 1
    :type mask_thickness: float, optional
    :return:region labels, same shape as input grid. 0 refers to the back ground
    :rtype: 3D numpy.array same shape and input grid
    """    
    from scipy import ndimage as ndi
    import numpy as np
    import sparse
    from skimage import segmentation
    from skimage.morphology import extrema
    peak_min =2
    max_indices = np.where(extrema.h_maxima(dgrid, peak_min))
    markers_sparse = sparse.COO(max_indices, 1, shape=dgrid.shape)

    probe = mask_thickness
    import dask.array as da
    reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int))
    return reg_labels, max_indices


def connections_from_regions(region_labels, ase_atoms):
    """Find the windows in a structure using the region_labels

    :param ase_atoms: ASE atoms object
    :type ase_atoms: ase.atoms
    :param rlpbc: Periodic/Non periodic region labels output from the make_labels_grid function
    :type rlpbc: 3D numpy.array
    :return: A tuple with two elements, a 2D array of connected regions and the (N_windows, 3) 2D array of coordinates of the calculated window centers     
    :rtype:  tuple
    """    
    from skimage.segmentation import find_boundaries
    import numpy as np
    import itertools
    from tqdm import tqdm
    from ase.io import read
    import ase
    # from sklearn.cluster import DBSCAN 
    # from sklearn.preprocessing import StandardScaler
    from icecream import ic
    from sklearn.cluster import AgglomerativeClustering

    A_unit = ase.geometry.complete_cell(ase_atoms.get_cell()).T
    
    def check_connected(regs):
        import numpy as np
        i = regs[0]
        j = regs[1]
    
        check_flag = np.logical_or(np.logical_and(region_labels == i, outer[j - 1]), np.logical_and(region_labels == j, outer[i - 1]))
        if np.sum(check_flag) > 0:
            # print("Edges found between "+str(i)+' and '+ str(j) + ' sum: '+str(np.sum(check_flag)))
            windices = np.vstack(np.where(check_flag)).T  # * (N,3)
            window_points = np.dot(A_unit, (windices / region_labels.shape).T).T  # * (N,3)
            # ic(window_points)
            # * No clustering
            # window_centers = np.mean(window_points, axis=0)
            
            # * DBScan
            # X_scaled = StandardScaler().fit_transform(window_points)
            # db = DBSCAN(eps=0.5, min_samples=10).fit(X_scaled)
            # window_centers = np.vstack([window_points[db.labels_==i].mean(axis=0) for i in np.unique(db.labels_) if i != -1])
            
            # * Agglomerative
            ag = AgglomerativeClustering(n_clusters=None, linkage='single',distance_threshold=1.0, compute_full_tree=True).fit(window_points)
            window_centers = np.vstack([window_points[ag.labels_==i].mean(axis=0) for i in np.unique(ag.labels_) if i != -1])

            return np.array([regs, True, window_centers], dtype='object') # returns a list of [[r1,r2], connected or not, array of window_centers]
        else:
            return np.array([regs, False, None])
    

    # Find the outer boundary of each individual region
    outer = np.stack([find_boundaries(region_labels == reg, mode='outer') for reg in range(1, np.max(region_labels) + 1)])
    
    # Check for connectedness between regions
    check_list = np.vstack(list(itertools.combinations(range(1, np.max(region_labels) + 1), 2))) # list of pairs of region labels
    connections_list = np.array([check_connected(c) for c in check_list], dtype='object') # now this step could cause an overflow
    # connected_regions = [[c[0], [c[2]] for c in connections_list if c[1] == True] # Make a new list with all the 
    return connections_list # for only the true connections use connections_list[:,[0,2]][np.where(connections_list[:,1]==True)] #np.vstack(check_list[np.where(flag_list[:, 0])]), np.round(np.vstack(flag_list[:, 1][np.where(flag_list[:, 0])]), decimals=2)



if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())

# """
#  Compute functions for mofography
#   ██████╗ ██████╗ ███╗   ███╗██████╗ ██╗   ██╗████████╗███████╗
#  ██╔════╝██╔═══██╗████╗ ████║██╔══██╗██║   ██║╚══██╔══╝██╔════╝
#  ██║     ██║   ██║██╔████╔██║██████╔╝██║   ██║   ██║   █████╗  
#  ██║     ██║   ██║██║╚██╔╝██║██╔═══╝ ██║   ██║   ██║   ██╔══╝  
#  ╚██████╗╚██████╔╝██║ ╚═╝ ██║██║     ╚██████╔╝   ██║   ███████╗
#   ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝      ╚═════╝    ╚═╝   ╚══════╝
                                                               
 
# """
# TODO Identify the different unique regions
# TODO Fix the issue with  the border regions getting cut off
# def canvas(with_attribution=True):
#     """
#     Placeholder function to show example docstring (NumPy format).

#     Replace this function and doc string for your own project.

#     Parameters
#     ----------
#     with_attribution : bool, Optional, default: True
#         Set whether or not to display who the quote is from.

#     Returns
#     -------
#     quote : str
#         Compiled string including quote and optional attribution.
#     """

#     quote = "The code is but a canvas to our imagination."
#     if with_attribution:
#         quote += "\n\t- Adapted from Henry David Thoreau"
#     return quote








def gpd_aabb(grid_points, ase_atoms, radii, n_neighbors=10, probe_radius=0.0):
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
    # from ase.data import vdw_radii
    from icecream import ic
    import pandas as pd
    import numpy as np

    framework_atom_positions = ase_atoms.get_positions()
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()).T # this needs to be an upper triangular matrix, essential column vectors, not row vectors
    
    box = Box.from_matrix(cell, dimensions =3) # Define the box for Freud
    ab  = AABBQuery(box=box, points=framework_atom_positions) # Compute the AABB Tree
    # Perform the nearest neighbor query with n-neighbors (10 by default)
    ab_r = ab.query(grid_points, query_args=dict(mode='nearest', num_neighbors=n_neighbors, exclude_ii=True))
    ab_r_nl = ab_r.toNeighborList(sort_by_distance=True) #convert query to neighborlist
    nls = ab_r_nl.point_indices.reshape(-1,n_neighbors)  # get the indices of the neighbors for each grid pt. to lookup the radii later
    dists =ab_r_nl.distances.reshape(-1,n_neighbors) # get the distances to the neighbors of each grid pt.
    # print(ab, ab_r, ab_r_nl, nls, [radii[nl] for nl in nls])#, np.min(dists - [radii[nl] for nl in nls],axis=1))
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
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()).T # cell matrix
    return da.dot(cell, gpoints.T).T, (nx,ny,nz) #return the actual coordinates

def dgrid_from_atoms(ase_atoms, radii=None, spacing=0.25, block_size=50000, n_neighbors=10, probe_radius=0):
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
    gpoints, shape = dask_grid_over_atoms(ase_atoms, spacing=spacing, chunksize=block_size) #mesh grid
    dgrid = []
    import pandas as  pd 
    if radii is None:
        vdw_radii = pd.Series({'H': 110.00000000000001,'He': 140.0,'Li': 182.0,'Be': 153.0,'B': 192.0,'C': 170.0,'N': 155.0,'O': 152.0,'F': 147.0,'Ne': 154.0,'Na': 227.0,'Mg': 173.0,'Al': 184.0,'Si': 210.0,'P': 180.0,'S': 180.0,'Cl': 175.0,'Ar': 188.0,'K': 275.0,'Ca': 231.0,'Sc': 215.0,'Ti': 211.0,'V': 206.99999999999997,'Cr': 206.0,'Mn': 204.99999999999997,'Fe': 204.0,'Co': 200.0,'Ni': 197.0,'Cu': 196.0,'Zn': 200.99999999999997,'Ga': 187.0,'Ge': 211.0,'As': 185.0,'Se': 190.0,'Br': 185.0,'Kr': 202.0,'Rb': 303.0,'Sr': 249.00000000000003,'Y': 231.99999999999997,'Zr': 223.0,'Nb': 218.00000000000003,'Mo': 217.0,'Tc': 216.0,'Ru': 213.0,'Rh': 210.0,'Pd': 210.0,'Ag': 211.0,'Cd': 218.00000000000003,'In': 193.0,'Sn': 217.0,'Sb': 206.0,'Te': 206.0,'I': 198.0,'Xe': 216.0,'Cs': 343.0,'Ba': 268.0,'La': 243.00000000000003,'Ce': 242.0,'Pr': 240.0,'Nd': 239.0,'Pm': 238.0,'Sm': 236.0,'Eu': 235.0,'Gd': 234.0,'Tb': 233.0,'Dy': 231.0,'Ho': 229.99999999999997,'Er': 229.0,'Tm': 227.0,'Yb': 225.99999999999997,'Lu': 224.00000000000003,'Hf': 223.0,'Ta': 222.00000000000003,'W': 218.00000000000003,'Re': 216.0,'Os': 216.0,'Ir': 213.0,'Pt': 213.0,'Au': 214.0,'Hg': 223.0,'Tl': 196.0,'Pb': 202.0,'Bi': 206.99999999999997,'Po': 197.0,'At': 202.0,'Rn': 220.00000000000003,'Fr': 348.0,'Ra': 283.0,'Ac': 247.00000000000003,'Th': 245.00000000000003,'Pa': 243.00000000000003,'U': 241.0,'Np': 239.0,'Pu': 243.00000000000003,'Am': 244.0,'Cm': 245.00000000000003,'Bk': 244.0,'Cf': 245.00000000000003,'Es': 245.00000000000003,'Fm': 245.00000000000003,'Md': 246.0,'No': 246.0,'Lr': 246.0,'Rf': np.nan,'Db': np.nan,'Sg': np.nan,'Bh': np.nan,'Hs': np.nan,'Mt': np.nan,'Ds': np.nan,'Rg': np.nan,'Cn': np.nan,'Nh': np.nan,'Fl': np.nan,'Mc': np.nan,'Lv': np.nan,'Ts': np.nan,'Og': np.nan})
        radii = vdw_radii[ase_atoms.get_chemical_symbols()].values/100.0
    
    for block in gpoints.blocks:
        dgrid.append(gpd_aabb(block.compute(), ase_atoms, radii, n_neighbors=n_neighbors,probe_radius=probe_radius )) #compute grid
    return da.hstack(dgrid).rechunk(chunks=block_size).reshape(shape) #return in correct shape



def dgrid_from_atoms_gpu(ase_atoms, radii=None, spacing=0.2, block_size=10000):
    """Code to compute the distance grid using the GPU

    :param data: ASE style atoms object
    :type data: ase.atoms
    :param spacing: grid spacing in A, defaults to 0.2
    :type spacing: float, optional
    :param block_size: block size of the grid points array, number of grid points to be prcessed by each core of the CPU while talking to the GPU, defaults to 5000
    :type block_size: int, optional
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
    fpoints = cp.asarray(ase_atoms.get_scaled_positions()  - 0.5)
    cell = cp.asarray(ase.geometry.complete_cell(ase_atoms.get_cell())).T
    # from ase.data import vdw_radii
    # radius_series = pd.Series({'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96,
    #                         'B': 0.84, 'C': 0.73, 'N': 0.71, 'O': 0.66,
    #                         'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41,
    #                         'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
    #                         'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76,
    #                         'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39,
    #                         'Mn': 1.50, 'Fe': 1.42, 'Co': 1.38, 'Ni': 1.24,
    #                         'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.20,
    #                         'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16,
    #                         'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75,
    #                         'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46,
    #                         'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
    #                         'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38,
    #                         'I': 1.39, 'Xe': 1.40, 'Cs': 2.44, 'Ba': 2.15,
    #                         'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01,
    #                         'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96,
    #                         'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89,
    #                         'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75,
    #                         'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44,
    #                         'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32,
    #                         'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40,
    #                         'At': 1.50, 'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21,
    #                         'Ac': 2.15, 'Th': 2.06, 'Pa': 2.00, 'U': 1.96,
    #                         'Np': 1.90, 'Pu': 1.87, 'Am': 1.80, 'Cm': 1.69})
    # radii = cp.asarray(radius_series[data.get_chemical_symbols()].values).reshape(-1, 1)
    if radii is None:
        vdw_radii = pd.Series({'H': 110.00000000000001,'He': 140.0,'Li': 182.0,'Be': 153.0,'B': 192.0,'C': 170.0,'N': 155.0,'O': 152.0,'F': 147.0,'Ne': 154.0,'Na': 227.0,'Mg': 173.0,'Al': 184.0,'Si': 210.0,'P': 180.0,'S': 180.0,'Cl': 175.0,'Ar': 188.0,'K': 275.0,'Ca': 231.0,'Sc': 215.0,'Ti': 211.0,'V': 206.99999999999997,'Cr': 206.0,'Mn': 204.99999999999997,'Fe': 204.0,'Co': 200.0,'Ni': 197.0,'Cu': 196.0,'Zn': 200.99999999999997,'Ga': 187.0,'Ge': 211.0,'As': 185.0,'Se': 190.0,'Br': 185.0,'Kr': 202.0,'Rb': 303.0,'Sr': 249.00000000000003,'Y': 231.99999999999997,'Zr': 223.0,'Nb': 218.00000000000003,'Mo': 217.0,'Tc': 216.0,'Ru': 213.0,'Rh': 210.0,'Pd': 210.0,'Ag': 211.0,'Cd': 218.00000000000003,'In': 193.0,'Sn': 217.0,'Sb': 206.0,'Te': 206.0,'I': 198.0,'Xe': 216.0,'Cs': 343.0,'Ba': 268.0,'La': 243.00000000000003,'Ce': 242.0,'Pr': 240.0,'Nd': 239.0,'Pm': 238.0,'Sm': 236.0,'Eu': 235.0,'Gd': 234.0,'Tb': 233.0,'Dy': 231.0,'Ho': 229.99999999999997,'Er': 229.0,'Tm': 227.0,'Yb': 225.99999999999997,'Lu': 224.00000000000003,'Hf': 223.0,'Ta': 222.00000000000003,'W': 218.00000000000003,'Re': 216.0,'Os': 216.0,'Ir': 213.0,'Pt': 213.0,'Au': 214.0,'Hg': 223.0,'Tl': 196.0,'Pb': 202.0,'Bi': 206.99999999999997,'Po': 197.0,'At': 202.0,'Rn': 220.00000000000003,'Fr': 348.0,'Ra': 283.0,'Ac': 247.00000000000003,'Th': 245.00000000000003,'Pa': 243.00000000000003,'U': 241.0,'Np': 239.0,'Pu': 243.00000000000003,'Am': 244.0,'Cm': 245.00000000000003,'Bk': 244.0,'Cf': 245.00000000000003,'Es': 245.00000000000003,'Fm': 245.00000000000003,'Md': 246.0,'No': 246.0,'Lr': 246.0,'Rf': np.nan,'Db': np.nan,'Sg': np.nan,'Bh': np.nan,'Hs': np.nan,'Mt': np.nan,'Ds': np.nan,'Rg': np.nan,'Cn': np.nan,'Nh': np.nan,'Fl': np.nan,'Mc': np.nan,'Lv': np.nan,'Ts': np.nan,'Og': np.nan})
        radii = cp.asarray(vdw_radii[ase_atoms.get_chemical_symbols()].values/100.0).reshape(-1,1)
    
    # * the vdw radii mendeleev is in picometers
        # radii = cp.asarray(vdw_radii[ase_atoms.get_atomic_numbers()].values/100).reshape(-1, 1)
    [nx, ny, nz] = (ase_atoms.get_cell_lengths_and_angles()[0:3] / spacing).astype(int) + 1
    #     print([nx,ny,nz])
    gpoints = da.stack(
        da.meshgrid(np.linspace(-0.5, 0.5, nx), np.linspace(-0.5, 0.5, ny), np.linspace(-0.5, 0.5, nz), indexing='ij'),
        -1).reshape(-1, 3).rechunk(block_size, 3)

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

    gpoints_da = gpoints.map_blocks(gpd, chunks=(block_size, 1))

    distance_grid = []
    from rich.progress import track
    for block in track(gpoints_da.blocks, total=gpoints_da.npartitions):
        # print("Working on block {} \n".format(i))
        distance_grid.append(block.compute())

    #  Free the gpu memeory before returning the output
    del cell, fpoints, radii
    cp._default_memory_pool.free_all_blocks()
    cp._default_pinned_memory_pool.free_all_blocks()
    return np.hstack(distance_grid).reshape(nx, ny, nz)



# def regions_from_dgrid(dgrid, mask_thickness=1, **kwargs):
#     """Computes the region segmentaion from the distance grid

#     :param dgrid: 3D distance grid
#     :type dgrid: numpy.array
#     # :param peak_min:minimum height of the peak in angstrom, defaults to 2.0
#     # :type peak_min: float, optional, defaults to 2.0 angstroms
#     :param relative_threshold: Minimum intensity of peaks, calculated as max(image) * relative_threshold, defaults to 0.1
#     :type relative_threshold: float, 
#     :param mask_thickness:thickness of the framework cover. Points within this distance from the framework surface are not segmented, defaults to 1
#     :type mask_thickness: float, optional
#     :return:region labels, same shape as input grid. 0 refers to the back ground
#     :rtype: 3D numpy.array same shape and input grid
#     """    
#     from scipy import ndimage as ndi
#     import numpy as np
#     import sparse
#     from .helper import interpolate_labels, find_coord_from_indices
#     from skimage import segmentation
#     from skimage.morphology import extrema
#     from skimage.feature import peak_local_max, corner_peaks
#     from sklearn.cluster import AgglomerativeClustering
#     # peak_min =2
#     # max_indices = np.vstack(np.where(extrema.h_maxima(dgrid, h=0.5)))
#     # max_indices = np.vstack(np.where(extrema.local_maxima(dgrid)))
#     # max_indices = np.array(peak_local_max(dgrid, **kwargs)).T
#     # max_indices = np.array(peak_local_max(dgrid, **kwargs)).T
#     # Finding the local maxima of the corner response function.
    
#     # ag = AgglomerativeClustering(n_clusters=None, linkage='average', distance_threshold=np.max(np.array(dgrid.shape))/5, compute_full_tree=True).fit(max_indices)
#     # max_indices= np.vstack([max_indices[ag.labels_==i].mean(axis=0) for i in np.unique(ag.labels_) if i != -1]).astype(int)
#     # print(max_indices)
#     # markers_sparse = sparse.COO(max_indices, 1, shape=dgrid.shape)
#     from icecream import ic 
#     probe = mask_thickness
#     import dask.array as da
#     import warnings
#     # reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int))

#     # max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1, **kwargs)).T
#     max_indices = np.array(peak_local_max(dgrid,  **kwargs)).T
#     markers_sparse = sparse.COO(max_indices, 1, shape=dgrid.shape)

#     probe = mask_thickness
#     import dask.array as da
#     reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int))


#     new_max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1)).T
#     if len(new_max_indices.T) < len(max_indices.T):
#         warnings.warn('Some maxima excluded as they fell inside the mask_thickness')
    
#     markers_sparse = sparse.COO(new_max_indices, 1, shape=dgrid.shape)

#     probe = mask_thickness
#     import dask.array as da
#     reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int))

#     #* Now interpolate the region labels and sort the maxima in the correct order    
#     sort_index = np.argsort([reg_labels[tuple(m)] for m in new_max_indices.T])
#     new_max_indices = (new_max_indices.T[sort_index]).T
#     # ir = interpolate_labels(reg_labels)
#     # new_max_indices = (new_max_indices.T[np.argsort(ir(new_max_indices.T/dgrid.shape))]).T
    
#     # # new_max_indices = []
#     # old_labels = np.unique(reg_labels)  # This is sorted aalreadylready
#     # number_of_regions = len(old_labels)
#     # new_labels = range(number_of_regions)

#     # # * Replace the old-labels with new
#     # for i in range(number_of_regions):
#     #     reg_labels[reg_labels== old_labels[i]] = new_labels[i]
#     return reg_labels, new_max_indices.T
def regions_from_dgrid(dgrid, mask_thickness=1, h=0.2, min_distance=2, compactness=0):
    """Computes the region segmentaion from the distance grid

    :param dgrid: 3D distance grid
    :type dgrid: numpy.array
    # :param peak_min:minimum height of the peak in angstrom, defaults to 2.0
    # :type peak_min: float, optional, defaults to 2.0 angstroms
    :param relative_threshold: Minimum intensity of peaks, calculated as max(image) * relative_threshold, defaults to 0.1
    :type relative_threshold: float, 
    :param mask_thickness:thickness of the framework cover. Points within this distance from the framework surface are not segmented, defaults to 1
    :type mask_thickness: float, optional
    :return:region labels, same shape as input grid. 0 refers to the back ground
    :rtype: 3D numpy.array same shape and input grid
    """    
    from scipy import ndimage as ndi
    import numpy as np
    import sparse
    # from .helper import interpolate_labels, find_coord_from_indices
    from skimage import segmentation
    from skimage.morphology import extrema
    from skimage.feature import peak_local_max, corner_peaks
    from sklearn.cluster import AgglomerativeClustering
    # peak_min =2
    # max_indices = np.vstack(np.where(extrema.h_maxima(dgrid, h=0.5)))
    # max_indices = np.vstack(np.where(extrema.local_maxima(dgrid)))
    # max_indices = np.array(peak_local_max(dgrid, **kwargs)).T
    # max_indices = np.array(peak_local_max(dgrid, **kwargs)).T
    # Finding the local maxima of the corner response function.
    
    # ag = AgglomerativeClustering(n_clusters=None, linkage='average', distance_threshold=np.max(np.array(dgrid.shape))/5, compute_full_tree=True).fit(max_indices)
    # max_indices= np.vstack([max_indices[ag.labels_==i].mean(axis=0) for i in np.unique(ag.labels_) if i != -1]).astype(int)
    # print(max_indices)
    # markers_sparse = sparse.COO(max_indices, 1, shape=dgrid.shape)
    from icecream import ic 
    probe = mask_thickness
    import dask.array as da
    import warnings
    # reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int))

    # max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1, **kwargs)).T
    # max_indices = np.array(peak_local_max(dgrid,  **kwargs)).T
    max_indices = np.vstack(np.where(extrema.h_maxima(dgrid,  h=h)))
    markers_sparse = sparse.COO(max_indices, 1, shape=dgrid.shape)

    probe = mask_thickness
    import dask.array as da
    reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int), compactness=compactness)



    # print(keep_these)
    # new_max_indices = max_indices[keep_these-1].T

    new_max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1)).T
    if len(new_max_indices.T) < len(max_indices.T):
        warnings.warn('Some maxima excluded as they fell inside the mask_thickness')
        # warnings.warn('Some maxima excluded as they the regions were less that 0.001 of the box')
    
    markers_sparse = sparse.COO(new_max_indices, 1, shape=dgrid.shape)

    probe = mask_thickness
    import dask.array as da
    reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int), compactness=compactness)

    new_max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1)).T
    #* Now interpolate the region labels and sort the maxima in the correct order    
    sort_index = np.argsort([reg_labels[tuple(m)] for m in new_max_indices.T])
    new_max_indices = (new_max_indices.T[sort_index]).T
    # ir = interpolate_labels(reg_labels)
    # new_max_indices = (new_max_indices.T[np.argsort(ir(new_max_indices.T/dgrid.shape))]).T
    
    # # new_max_indices = []
    # old_labels = np.unique(reg_labels)  # This is sorted aalreadylready
    # number_of_regions = len(old_labels)
    # new_labels = range(number_of_regions)

    # # * Replace the old-labels with new
    # for i in range(number_of_regions):
    #     reg_labels[reg_labels== old_labels[i]] = new_labels[i]

    # Lets clean out the regions that are too small 
    # props = regionprops(reg_labels)
    # print([(p.label, p.area/dgrid.size) for p in props])
    # len([(p.label, p.area/dgrid.size) for p in props if p.area > dgrid.size/1000])
    # keep_these = [p.label for p in props if p.area > dgrid.size/1000]
    # maxima2 = np.vstack([new_max_indices.T[k-1] for k in keep_these])
    # markers_sparse = sparse.COO(maxima2.T, 1, shape=dgrid.shape)

    # probe = mask_thickness
    # import dask.array as da
    # reg_labels = segmentation.watershed(-dgrid, ndi.label(markers_sparse.todense())[0], mask=(dgrid>probe).astype(int), compactness=False)

    # new_max_indices = np.array(peak_local_max(dgrid, labels=reg_labels, num_peaks_per_label=1)).T
 

    return reg_labels, new_max_indices.T

def connections_from_regions():

    """
    It takes a 3D array of region labels, a list of maxima, an ASE atoms object, and a distance
    threshold, and returns a list of connections between regions, and the centers of the windows between
    them
    
    :param region_labels: The labels of the regions, as returned by the skimage.measure.label function
    :param maxima: the maxima of the image
    :param ase_atoms: ASE atoms object
    :param wall_windows: If True, then the windows on the wall are also calculated, defaults to True
    (optional)
    :param distance_threshold: The distance threshold for the agglomerative clustering
    :return: A list of connections between regions.
    """


    from skimage.segmentation import find_boundaries
    import numpy as np
    import itertools
    import ase
    # from sklearn.cluster import DBSCAN 
    # from sklearn.preprocessing import StandardScaler
    from icecream import ic
    from sklearn.cluster import AgglomerativeClustering
    from .helper import interpolate_labels, get_fractional_coordinates
    A_unit = ase.geometry.complete_cell(ase_atoms.get_cell()).T
    from icecream import ic
    def check_connected(regs):
        import numpy as np
        from sklearn.cluster import AgglomerativeClustering
        i = regs[0]
        j = regs[1]
        ic.disable()
        ic(regs)
    
        check_flag = np.logical_or(np.logical_and(region_labels == i, outer.get(j)), np.logical_and(region_labels == j, outer.get(i)))
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
            
            # if len(window_points) >1:
            # * Agglomerative
            ag = AgglomerativeClustering(n_clusters=None, linkage='single',distance_threshold=minimum_window_separation, compute_full_tree=True).fit(window_points)
            window_centers = np.vstack([window_points[ag.labels_==i].mean(axis=0) for i in np.unique(ag.labels_) if i != -1])
        # else:
            # window_centers = window_points
            window_fractional = get_fractional_coordinates(window_centers,ase_atoms)
            window_radii = np.array(dinterp(window_fractional))
            return np.array([regs, True,'internal', window_centers, window_radii], dtype='object') # returns a list of [[r1,r2], connected or not, array of window_centers]
        else:
            return np.array([regs, False,'internal', None, None], dtype='object')
    
    dinterp = interpolate_labels(dgrid)
    # * Find the outer boundary of each individual region
    regions_to_check = np.delete(np.unique(region_labels),0)
    ic(regions_to_check)
    outer = {reg:find_boundaries(region_labels == reg, mode='outer') for reg in regions_to_check}
    ic(len(outer))
    # * Check for connectedness between regions
    check_list = np.vstack(list(itertools.combinations(regions_to_check, 2))) # list of pairs of region labels
    connections_list = [check_connected(c) for c in check_list] # now this step could cause an overflow
    # connected_regions = [[c[0], [c[2]] for c in connections_list if c[1] == True] # Make a new list with all the 
    
    #* Windows on the wall are to be calculated
    if wall_windows:
        wwl = get_wall_windows2(regions=region_labels, maxima=maxima, dgrid=dgrid, ase_atoms=ase_atoms, wall_thickness=wall_thickness)
        if len(wwl)>0:
            connections_list.append(wwl)
    return np.vstack(connections_list)
    # return np.vstack([c for c in np.vstack(connections_list) if c[1]])# for only the true connections use connections_list[:,[0,2]][np.where(connections_list[:,1]==True)] #np.vstack(check_list[np.where(flag_list[:, 0])]), np.round(np.vstack(flag_list[:, 1][np.where(flag_list[:, 0])]), decimals=2))

# def get_wall_windows(regions, maxima, ase_atoms):
#     """
#     > It finds the windows in the walls of the simulation box, and returns the indices of the regions
#     that are in the window
    
#     :param regions: the regions array
#     :param maxima: the maxima of the regions
#     :param ase_atoms: ASE atoms object
#     """
#     import numpy as np
#     from .helper import find_coord_from_indices
#     wall_windows = []
#     shape =regions.shape
#     for a in [0,1,2]:
        
#         llr = regions.take(indices=0 ,axis=a) #* left wall
#         rrr = regions.take(indices=shape[a]-1 ,axis=a)  #* right wall

#         ll = [l for l in np.delete(np.unique(llr),0) if maxima[l-1][a]!=0] #* unique regions in left wall whose maxima are not on the wall
#         rr = [ r for r in np.delete(np.unique(rrr),0) if maxima[r-1][a]!=shape[a]-1]#* unique regions in the right wall whose maxima are not on the wall
#         left_right_map = np.vectorize(dict(zip(llr.flatten(), rrr.flatten())).get)
#         right_left_map = np.vectorize(dict(zip(rrr.flatten(), llr.flatten())).get)
#         ll_2 = [l for l in ll if maxima[left_right_map(l)-1][a]!=shape[a]-1]
#         rr_2 = [r for r in rr if maxima[right_left_map(r)-1][a]!=0]
        
#         if len(ll_2)>0:
#             #* Mean of the indices of each windowed region in the left nad right wall, only if the opposite wall maxima is not on the wall
#             wls = np.vstack([np.vstack(np.where(llr==l)).mean(axis=1) for l in ll_2 ])
#             #* Add the third index value, in the correct position, 0 for left wall, shape[axis] - 1 for the right wall
#             wls = np.insert(wls,[a], np.zeros((len(wls),1)), axis=1)
#             wls = find_coord_from_indices(wls, shape=shape, ase_atoms=ase_atoms)
#             # * use the map to come up with the correct touching regions
#             ll_p =  left_right_map(ll_2) #[np.unique(rrr[np.where(llr==l)])[0] for l in ll ]
#             left_windows  = [[np.array([r1, r2]), True,'periodic', np.array(wc).reshape(-1,3)] for r1, r2,wc in zip(ll, ll_p, wls)]
#             wall_windows.append(left_windows) 
    
#         if len(rr_2)>0:    
#             wrs = np.vstack([np.vstack(np.where(rrr==r)).mean(axis=1) for r in rr_2 ])
#             wrs = np.insert(wrs,[a], (shape[a]-1)*np.ones((len(wrs),1)), axis=1)
#             wrs = find_coord_from_indices(wrs, shape=shape, ase_atoms=ase_atoms)
#             rr_p = right_left_map(rr_2) # [np.unique(llr[np.where(rrr==r)])[0] for r in rr ]
#             right_windows = [[np.array([r1, r2]), True,'periodic', np.array(wc).reshape(-1,3)] for r1, r2,wc in zip(rr, rr_p, wrs)]
#             #* Append the newly found window indices
#             wall_windows.append(right_windows)
    
#     if len(wall_windows)>0:
#         wall_windows =  np.vstack(wall_windows)
#     return wall_windows


def get_region_props_for_regions(regions):
    """
    It takes a list of regions and returns a list of region properties
    
    :param regions: The labeled image
    """
    from skimage.measure import regionprops
    return regionprops(regions)

def add_maxima_to_rag(G, maxima,  maxima_radii, shape, ase_atoms):
    """
    > This function takes a graph, a list of maxima, a list of maxima radii, a shape, and an ase_atoms
    object and adds the maxima and maxima radii to the graph. 
    
    Let's look at the function in more detail. 
    
    The first thing we do is import the `find_coord_from_indices` function from the `helper` module. 
    
    Next, we check to make sure that the number of nodes in the graph is equal to the number of maxima.
    If not, we raise a `ValueError`. 
    
    If the number of nodes in the graph is equal to the number of maxima, we loop through the nodes and
    maxima and add the maxima indices, cartesian coordinates and maxima radii to the graph. 
    

    
    :param G: the graph
    :param maxima: the indices of the maxima (2nd output of regions_from_grid routine)
    :param maxima_radii: the radius of the maxima
    :param shape: the shape of the image
    :param ase_atoms: The atoms object from ASE
    :return: The graph with the maxima added to the nodes
    """
    

    from .helper import find_coord_from_indices
    if len(G.nodes()) != len(maxima):
        raise ValueError('The number of nodes in the graph and the number of maxima do not match')
    else:
        # regionprops = get_region_props_for_regions(G.nodes())
        for n, m in zip(G.nodes(), maxima):
            G.nodes[n]['maxima'] = find_coord_from_indices(maxima[G.nodes[n]['labels'][0]-1], shape=shape, ase_atoms=ase_atoms) #* add the coordinates of the maxima to the node
            G.nodes[n]['maxima_indices'] = maxima[G.nodes[n]['labels'][0]-1] #* add the indcies of the maxima 
            G.nodes[n]['maxima_radii'] =maxima_radii[G.nodes[n]['labels'][0]-1] #* add the indcies of the maxima 
    return G
# def add_size_dist_to_graph(G, dgrid):

def add_psd_to_rag(G, region_props, bins, ase_atoms, shape):
    """
    > This function takes a graph, a list of region properties, a list of bins, an ase atoms object, and
    a shape, and adds the psd to the graph. 
    
    For each region, we calculate the distances between the centroid and the coordinates. 
    
    Then we calculate the histogram of the distances. 
    
    Then we calculate the bin centers. 
    
    Then we add the bin centers and the histogram to the graph. 
    
    Then we return the graph. 
    
    Note that these distance histograms or centroids do not account for periodic boundary conditions.
    To include periodic boundary conditions, use the add_psd_to_rag_with_pbc function.
    
    :param G: the graph
    :param region_props: a list of dictionaries, each dictionary contains the properties of a region
    :param bins: The number of bins to use for the histogram
    :param ase_atoms: The atoms object from ASE
    :param shape: shape of the image
    :return: The graph with the psd_x and psd_y attributes added to the nodes.
    """


    import numpy as np 
    from .helper import find_coord_from_indices 
    for region_prop in region_props:    
        distances= np.linalg.norm(find_coord_from_indices(region_prop['centroid'] - region_prop['coords'], shape=shape, ase_atoms=ase_atoms), axis=1)

        G.nodes[region_prop['label']]['psd_y'], bin_edges = np.histogram(distances, bins=bins, density=True, )
        # bin_edges = np.histogram(region_prop['coords']-region_prop['centroid'], bins=bins, normed=True, density=1)[1]
        bin_centers = (bin_edges[1:]+bin_edges[:-1])/2
        G.nodes[region_prop['label']]['psd_x'] = bin_centers
    return G

def get_rag_from_regions_and_grid(regions, dgrid, ase_atoms, maxima=None):
    """
    > It takes a grid of regions and a grid of distances and returns a region adjacency graph with the centroid of each
    region as a node property
    
    :param regions: the regions that we computed from the distance grid
    :param dgrid: the distance grid
    :param ase_atoms: the atoms object from ase
    :return: A region adjacency graph
    """
    
    import numpy as np
    # * compute a region adjacency graph from skimage
    from skimage.future import graph 
    from .helper import find_coord_from_indices
    gr_mean = graph.rag_mean_color(dgrid, regions)
    
    from skimage.measure import regionprops 
    # * add the centroid property to the rag nodes
    region_props = regionprops(regions)
    for region_prop in region_props:
        #add the coordinates of the centroid to the node
        gr_mean.nodes[region_prop['label']]['centroid'] = find_coord_from_indices( indices = np.array(region_prop['centroid']), ase_atoms=ase_atoms,shape=dgrid.shape)
        gr_mean.nodes[region_prop['label']]['centroid_indices'] = region_prop['centroid'] # * add the indices of the centroid to the node
    gr_mean.remove_node(0)
    if maxima is not None:
        maxima_radii = dgrid[tuple(maxima.T)] # * radii of the maxima
        gr_mean = add_maxima_to_rag(gr_mean, maxima, maxima_radii, regions.shape, ase_atoms)
        for edge in gr_mean.edges():
            gr_mean.edges[edge]['weight'] = np.linalg.norm(gr_mean.nodes[edge[0]]['maxima'] - gr_mean.nodes[edge[1]]['maxima'])
    else:
        for edge in gr_mean.edges():
            gr_mean.edges[edge]['weight'] = np.linalg.norm(gr_mean.nodes[edge[0]]['centroid'] - gr_mean.nodes[edge[1]]['centroid'])
    
    return gr_mean

def get_periodic_map(regions, maxima, maxima_radii, return_conn=False, **kwargs):
    """
    It takes a 3D array of integers representing region indices that is obtained from region segmentation, a list of maxima, and a list of maxima radii, and returns a dictionary
    that maps each region to the nearest maximum.
    
    :param regions: a 3D array of integers, where each integer is a region index, obtained from regions_from_dgrid function
    :param maxima: the 3D indices of the maxima, obtained from regions_from_dgrid function
    :param maxima_radii: the list of radii of the maxima (in Angstrom)
    :param return_conn: if True, returns the connectivity groups of the regions, defaults to False
    (optional)
    :return: a dictionary with the keys being the starting region indices and the values being the periodic group of
    the regions.

    then to apply this periodic map to regions, use: 
        np.vectorize(periodic_map.get)(regions)
    """
    
    
    # r1 = apply_pbc(regions, maxima, **kwargs)
    if return_conn:
        return apply_pbc(regions, maxima, maxima_radii, return_conn=return_conn, **kwargs)
    else:
        return dict(zip(regions.flatten(), apply_pbc(regions, maxima, maxima_radii, **kwargs).flatten()))



def apply_periodic_map(input_array, periodic_map):
    """
    It takes an array of integers and a dictionary that maps integers to integers, and returns an array
    of integers that is the result of applying the dictionary to the input array
    
    :param input_array: a numpy array of integers
    :param periodic_map: a dictionary that maps the input values to the output values
    :return: the value of the key in the dictionary.
    """
    import numpy as np
    return np.vectorize(periodic_map.get)(input_array)
# if __name__ == "__main__":
#     # Do something if this file is invoked on its own
#     print(canvas())
# def graph_from_connections(connections, maxima, only_use_internal=True):
    
#     """ Create a basic networkx graph (with no properties annotated) from the array of connections
#     :type connections: np.array of shape (Nregions-Choose-2 ,1) where Nregions is the number of unique regions. Each row has 
#     three elements, the region indices [r1, r2], a boolean flag showing whether they are connected, and window coordinate(s) if the regions are connected, and None is they are not.
#     :param connections: numpy array

#     :raises:

#     :rtype: a graph object of the class networkx.Graph
#     """
#     import networkx as nx 
#     import numpy as np
#     G = nx.Graph()
#     # number_of_regions = np.max([c[0] for c  in connections])
#     number_of_nodes = len(maxima) #np.max([c[0] for c  in connections])
#     # print(number_of_nodes)
#     G.add_nodes_from(range(1, number_of_nodes+1)) 
#     if only_use_internal:# range stops one short
#         G.add_edges_from([c[0] for c in connections if np.logical_and(c[1], c[2] =='internal')])
#     else:
#         G.add_edges_from([c[0] for c in connections if c[1]])

#     return G

def add_props_to_graph(G, node_props=None, node_prop_names=None, edge_props=None, edge_prop_names=None):
    """
    The function takes in a graph, and a list of node properties and edge properties, and adds them to the graph as node and
    edge attributes.


    :param G: the graph
    :param node_props: a list of lists, where each list contains the node properties for a given node property
    :param node_prop_names: list of names of node properties
    :param edge_props: a list of lists, where each list contains the edge properties for a given edge property
    :param edge_prop_names: list of strings, names of edge properties
    :return: A graph with the node and edge properties added to it.
    """


    import numpy as np
    import networkx as nx 
    if np.any(node_props) != None:
        for index in  range(len(node_props)):
            d = {i+1:k for i,k in enumerate(node_props[index])}
            nx.set_node_attributes(G, d, name=node_prop_names[index])
    
    if np.any(edge_props) != None:
        for index in  range(len(edge_props)):
            d = {i:k for i,k in zip(list(G.edges()), edge_props[index]) }
            nx.set_edge_attributes(G, d, name=edge_prop_names[index])
    
    return G


def AABB_on_atoms(ase_atoms, points=None):
    """
    > This function takes an ASE atoms object and returns a freud AABBQuery object
    
    :param ase_atoms: ASE atoms object
    :param points: the positions of the atoms
    :return: The AABB tree
    """
    from freud.locality import AABBQuery
    from freud.box import Box
    import ase 
    # from ase.data import vdw_radii


    if points is None:
        points = ase_atoms.get_positions()
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()).T # this needs to be an upper triangular matrix, essential column vectors, not row vectors
    
    box = Box.from_matrix(cell, dimensions =3) # Define the box for Freud
    ab  = AABBQuery(box=box, points=points) # Compute the AABB Tree
    # Perform the nearest neighbor query with n-neighbors (10 by default)
        
    return ab
def get_wall_windows2(regions, maxima, dgrid, ase_atoms, wall_thickness=2):

    import numpy as np
    from icecream import ic
    from .helper import find_coord_from_indices, interpolate_labels
    wall_windows = []
    shape =regions.shape
    dinterp = interpolate_labels(dgrid)
    for a in [0,1,2]:
        ic.disable()
        ic(a)
        llr = regions.take(indices=0 ,axis=a) #* left wall
        rrr = regions.take(indices=shape[a]-1 ,axis=a)  #* right wall

        map_list = np.unique(list(zip(llr.flatten(), rrr.flatten())), axis=0)
        map_list = np.vstack([m for m in map_list if np.all(m!=0)])
        flags = np.vstack([[maxima[k - 1][a]>wall_thickness, maxima[v-1][a]<(shape[a]-wall_thickness)] for k,v in map_list])    
        flags = np.logical_and(flags[:,0],flags[:,1])
        # window_indices = [np.vstack(np.where(llr==m[0])).mean(axis=1) for i,m in enumerate(map_list) if flags[i]]
        wregs_l  = [(np.array([m[0], m[1]]), True,'periodic') for i, m in enumerate(map_list) if flags[i]]
        wregs_r  = [(np.array([m[1], m[0]]), True,'periodic') for i, m in enumerate(map_list) if flags[i]]
        window_indices = [np.vstack(np.where(llr==m[0])).mean(axis=1).T for i, m in enumerate(map_list) if flags[i]]
        ic(window_indices)
        if np.any(flags):
            window_indices_l = np.insert(np.vstack(window_indices), [a], np.zeros((len(window_indices),1)), axis=1  )
            window_indices_r = np.insert(np.vstack(window_indices), [a], (shape[a] -1)*np.ones((len(window_indices),1)), axis=1 )
            window_radii_l = np.vstack(dinterp(window_indices_l/regions.shape))
            window_radii_r = np.vstack(dinterp(window_indices_r/regions.shape))
            window_coords_l = find_coord_from_indices(window_indices_l, shape=shape, ase_atoms=ase_atoms)
            window_coords_r = find_coord_from_indices(window_indices_r, shape=shape, ase_atoms=ase_atoms)
            ic(window_coords_l)
            ic(window_coords_r)
            window_rows_l = np.vstack([[*w, np.array([ window_coords_l[i] ]), window_radii_l[i]] for i,w in enumerate(wregs_l)])
            window_rows_r = np.vstack([[*w, np.array([ window_coords_r[i] ]), window_radii_r[i]] for i,w in enumerate(wregs_r)])
            wall_windows.append(window_rows_l) 
            wall_windows.append(window_rows_r) 
            
    if len(wall_windows)>0:
        wall_windows =  np.vstack(wall_windows)
    return wall_windows

# def apply_pbc(regions, maxima, maxima_radii, wall_thickness=2, return_conn =False):

#     import copy
#     # a = 0 
#     import numpy as np
#     from icecream import ic
#     ic.disable()
#     regions_pbc = copy.deepcopy(regions)
#     shape = regions.shape
#     import networkx as nx
#     g1= nx.Graph()
#     g1.add_nodes_from(range(len(np.unique(regions))))

#     for a in [0,1,2]:
        
#         ic(a)
#         llr = regions_pbc.take(indices=0 ,axis=a) #* left wall
#         rrr = regions_pbc.take(indices=shape[a]-1 ,axis=a)  #* right wall

#         # ic(dict(zip(llr.flatten(), rrr.flatten())))
#         # ic(dict(zip(rrr.flatten(),llr.flatten())))
#         # dlist = [dict(sorted(zip(llr.flatten(), rrr.flatten()))),dict(sorted(zip(rrr.flatten(),llr.flatten())))]

#         #* Check left vs right wall
#         # map_dict = dlist[0]
#         # map_dict.pop(0)
#         # ic(map_dict)


#         map_list = np.unique(list(zip(llr.flatten(), rrr.flatten())),axis=0)
#         map_list = map_list[np.all(map_list!=0,axis=1)]
#         ic(map_list)
#         if len(map_list)> 0:
#             reg_numbers = np.vstack([(k, v) for k,v in map_list])
#             flags = np.vstack([[maxima[k-1][a]<wall_thickness, maxima[v-1][a]>(shape[a]-wall_thickness)] for k,v in map_list])
#             # flags = np.logical_or(flags[:,0],flags[:,1,])
#             # flags = np.logical_or(flllags[:,0],flags[:,1,])
#             flags1 = np.logical_and(flags[:,0], flags[:,1])
#             flags2 = np.logical_xor(flags[:,0], flags[:,1])
#             reg_radii = np.vstack([[maxima_radii[r[0]-1], maxima_radii[r[1]-1]] for r in reg_numbers])
#             flags3 = np.vstack([reg_radii[i][np.where(f)]<reg_radii[i][np.where(np.invert(f))][0] if flags2[i] else False for i,f in enumerate(flags)])
#             # flags3 = np.vstack([ True if flags2[i] else False for i,f,r  in enumerate(list(zip(flags, reg_numbers)))])
#             # reg_radii = [maxima_radii[r[0]-1], maxima_radii[r[1]-1] for r in reg_numbers]
            
#             ic(flags)
#             flags = [flags1[i] if flags1[i] else flags2[i] if flags3[i] else False for i,f in enumerate(flags)]
#             ic(flags1)    
#             ic(flags2)    
#             ic(flags3)    
#             ic(flags)       # flags = [True if flags1[i] else True if np.logical_and(flags2[i], flags3[i] ) else False for i,f in enumerate(flags)]
#             g1.add_edges_from(reg_numbers[np.where(flags)])
#             # ic(flags) 

#             # map_dict = dlist[1]
#             # map_dict.pop(0)
#             # ic(map_dict)
#             # reg_numbers = np.vstack([(k, v) for k,v in map_dict.items()])
#             # flags = np.vstack([[maxima[v-1][a]<wall_thickness, maxima[k-1][a]>(shape[a]-wall_thickness)] for k,v in map_dict.items()])
#             # # flags = np.logical_or(flags[:,0],flags[:,1,])
#             # ic(flags)
#             # flags1 = np.logical_and(flags[:,0], flags[:,1])
#             # flags2 = np.logical_xor(flags[:,0], flags[:,1])
#             # reg_radii = np.vstack([[maxima_radii[r[0]-1], maxima_radii[r[1]-1]] for r in reg_numbers])
#             # flags3 = np.vstack([reg_radii[i][np.where(f)]<reg_radii[i][np.where(np.invert(f))] if flags2[i] else False for i,f in enumerate(flags)])
#             # # flags3 = np.vstack([maxima_radii[r[f]-1] < maxima_radii[r[np.invert(f)]-1] if flags2[i] else False for i,f,r  in enumerate(list(zip(flags, reg_numbers)))])
#             # # reg_radii = [maxima_radii[r[0]-1], maxima_radii[r[1]-1] for r in reg_numbers]
            
#             # flags = [flags1[i] if flags1[i] else flags2[i] if flags3[i] else False for i,f in enumerate(flags)]
#             # g1.add_edges_from(reg_numbers[np.where(flags)])

#             # ic(flags1)    
#             # ic(flags2)    
#             # ic(flags3)    
#             # ic(flags)    
#     conn = nx.connected_components(g1)
#     # ic(list(conn))
#     # ic(list(G.edges()))
#     # ic([c for c in list(conn)])
#     for i,c in enumerate(conn):
#         for r in c:
#             regions_pbc[np.where(regions==r)] = i
#             ic('replace {0} with {1}'.format(r, i))

#     # if reindex:
#     #     # * Re-index the regions so that it runs from zero to N-regions.
#     #     old_labels = np.unique(regions_pbc)  # This is sorted aalreadylready
#     #     number_of_regions = len(old_labels)
#     #     new_labels = range(number_of_regions)

#     #     # * Replace the old-labels with new
#     #     for i in range(number_of_regions):
#     #         regions_pbc[regions_pbc== old_labels[i]] = new_labels[i]
#     if return_conn:
#         return list(nx.connected_components(g1))[1:]
#     else:
#         return regions_pbc


def do_it_all(ase_atoms,spacing=0.25):
    """
    > It takes an ase atoms object of the mof, and returns a list of objects that are useful for analyzing the pore
    network of the material
    
    :param ase_atoms: the atoms object
    :param spacing: the spacing of the grid
    :return: updated_atoms, grid,regions, maxima, connections, pgraph
    """
    from .helper import find_coord_from_indices
    from .analyze import update_with_covalent_bonds
    updated_atoms =update_with_covalent_bonds(ase_atoms)
    grid = dgrid_from_atoms(ase_atoms=ase_atoms, spacing=spacing)
    regions, maxima =regions_from_dgrid(grid.compute(), mask_thickness=1.5, exclude_border=False, min_distance=int(2/spacing), threshold_abs=1)
    connections = connections_from_regions(maxima=maxima, region_labels=regions,ase_atoms=ase_atoms,dgrid=grid.compute(),minimum_window_separation=1,wall_windows=True)
    pgraph = graph_from_connections(connections, maxima, only_use_internal=True)
    pgraph = add_props_to_graph(pgraph, node_props = [find_coord_from_indices(maxima, shape=grid.shape, ase_atoms=updated_atoms)], node_prop_names=['node_coords'])
    return [updated_atoms, grid,regions, maxima, connections, pgraph]

def add_centroids_to_rag_with_pbc(rag, props, pbc_groups, ase_atoms, shape):
    """
    > For each group of nodes that are connected by periodic boundary conditions, find the centroid of
    the group and add it to the graph as a node attribute. 
    
    The function takes in a graph, a list of properties, a list of groups of nodes that are connected by
    periodic boundary conditions, an ASE atoms object, and the shape of the image. 
    
    The function first makes a copy of the graph. Then, it loops through each group of nodes that are
    connected by periodic boundary conditions. For each group, it finds the indices of the properties
    that correspond to the nodes in the group. If there is more than one node in the group, it finds the
    centroid of the group by taking the circular mean of the coordinates of the nodes in the group. If
    there is only one node in the group, it finds the centroid of the node by taking the mean of the
    coordinates of the node. Finally, it adds the centroid of the group
    
    :param rag: the Region Adjacency Graph object from skimage (skimage.future.graph.rag class)
    :param props: a list of dictionaries, each dictionary contains the label and the coordinates of the
    centroid of the cluster
    :param pbc_groups: a list of lists, where each list contains the labels of the atoms that are
    connected by periodic boundary conditions
    :param ase_atoms: the atoms object from ASE
    :param shape: shape of the image
    :return: The rag is being returned with the centroid_pbc added to the nodes.
    """


    from scipy import stats 
    import numpy as np 
    import copy
    from .helper import find_coord_from_indices
    G = copy.deepcopy(rag)
    prop_labels = np.array([prop['label'] for prop in props])
    for group in pbc_groups[1:]: #* 0 is always the background
        prop_indices = [np.where(prop_labels==g)[0][0] for g in group]
        # print(group, prop_indices)
        if len(group)>1:
            coords = np.vstack([props[pi]['coords'] for pi in prop_indices])/shape 
            # print(coords.shape)
            centroid = stats.circmean(coords, axis=0, low=0, high=1)
            for i, g in enumerate(group):
                G.nodes[props[prop_indices[i]]['label']]['centroid_pbc'] = find_coord_from_indices(indices = centroid*shape, ase_atoms=ase_atoms, shape=shape)
        else:
            pi = np.where(prop_labels==list( group )[0])[0][0]
            # print(pi)
            coords = props[pi]['coords']/shape
            centroid = np.mean(coords, axis=0)
            G.nodes[props[pi]['label']]['centroid_pbc'] = find_coord_from_indices(indices = centroid*shape, ase_atoms=ase_atoms, shape=shape)
    return G

def add_psd_to_rag_with_pbc(rag, props, pbc_groups, bins, ase_atoms, shape):
    """
    > It takes a rag, a list of properties, a list of pbc groups, a list of bins, an ase atoms object
    and a shape and returns a region adjacency graph with pbc properties added to it. 
    
    Let's see how we can use this function. 
    
    First, we need to get the pbc groups. 
    
    We can do that using the `get_pbc_groups` function and setting return_conn option to True (it is False by default)
    
    Let's see how we can use it.
    
    :param rag: the region adjacency graph
    :param props: list of dictionaries containing the properties of each region
    :param pbc_groups: list of lists of labels of regions that are connected by periodic boundary
    conditions
    :param bins: the bins for the histogram
    :param ase_atoms: ASE atoms object
    :param shape: shape of the image
    :return: The graph with the pbc properties added to the nodes.
    """


    def get_pbc_distances(coords, centroid_pbc, ase_atoms):
        """
        > The function takes the coordinates of the grid points in the region, the centroid of the region, and
        the atoms object, and returns the distances of the grid points from the periodic region centroid
        
        :param coords: the coordinates of the atoms in the cluster
        :param centroid_pbc: the centroid of the cluster in the PBC frame
        :param ase_atoms: ASE atoms object
        :return: The distances between the centroid and the atoms in the cluster.
        """
        
        from ase.geometry import complete_cell 
        cell = complete_cell(ase_atoms.get_cell()).T
        # cell_inv = np.linalg.inv(cell) 
        df = coords - centroid_pbc
        df = (df - np.round(df))
        distances = np.linalg.norm(np.dot(cell, df.T).T, axis=1)
        return distances

    from scipy import stats 
    import numpy as np 
    import copy
    from .helper import find_coord_from_indices
    # shape = regions.shape
    G = copy.deepcopy(rag)
    prop_labels = np.array([prop['label'] for prop in props])
    for group in pbc_groups: 
        prop_indices = [np.where(prop_labels==g)[0][0] for g in group]
        # print(group, prop_indices)
        if len(group)>1:
            #* fractional coordinates of all the points present in this group
            coords = np.vstack([props[pi]['coords'] for pi in prop_indices])/shape 
            # print(coords.shape)
            #* pbc centroid of the group in fractional coordinates
            centroid = stats.circmean(coords, axis=0, low=0, high=1)
            #* distances for the group in cartesian units
            distances = get_pbc_distances(coords, centroid, ase_atoms)
            # * histogram of the distances and bin centers
            psd_y_pbc, bine = np.histogram(distances, bins=bins, density=True)
            psd_x_pbc = (bine[1:] + bine[:-1])/2
            # * Now update the property of all the nodes in the graph with the pbc psd and cetnroid properties
            for i, g in enumerate(group):
                G.nodes[props[prop_indices[i]]['label']]['centroid_pbc'] = find_coord_from_indices(indices = centroid*shape, ase_atoms=ase_atoms, shape=shape)
                G.nodes[props[prop_indices[i]]['label']]['psd_x_pbc'] = psd_x_pbc
                G.nodes[props[prop_indices[i]]['label']]['psd_y_pbc'] = psd_y_pbc
        else:
            pi = np.where(prop_labels==list( group )[0])[0][0]
            # print(pi)
            #* get coordinates and centroids
            coords = props[pi]['coords']/shape
            centroid = np.mean(coords, axis=0)
            # * histogram of the distances and bin centers
            distances = get_pbc_distances(coords, centroid, ase_atoms)
            psd_y_pbc, bine = np.histogram(distances, bins=bins, density=True)
            psd_x_pbc = (bine[1:] + bine[:-1])/2

            # * Now update the property of all the nodes in the graph with the pbc psd and cetnroid properties
            G.nodes[props[pi]['label']]['centroid_pbc'] = find_coord_from_indices(indices = centroid*shape, ase_atoms=ase_atoms, shape=shape)
            G.nodes[props[pi]['label']]['psd_x_pbc'] = psd_x_pbc
            G.nodes[props[pi]['label']]['psd_y_pbc'] = psd_y_pbc
    
    return G
def add_clustering_to_rag(rag2,  return_cluster_obj = False, no_pbc =False, **kwargs):
    """
    > We take the RAG, get the PSD data for each node, and then use agglomerative clustering to cluster
    the nodes.
    
    :param rag2: the RAG object
    :param distance_threshold: The linkage distance at which to cut the dendrogram to determine the
    clusters
    :param return_cluster_obj: if True, returns the clustering object as well as the rag, defaults to
    False (optional)
    :return: The rag2 object with the cluster labels added to the nodes.
    """
    
    import numpy as np 
    import networkx as nx 

    # *  get the data array, while reordering to match the node order   
    if no_pbc:
        X = np.vstack(dict(sorted(nx.get_node_attributes(rag2, 'psd_y').items())).values())
    else:
        X = np.vstack(dict(sorted(nx.get_node_attributes(rag2, 'psd_y_pbc').items())).values())

    # * do the clustering using agglomerative clustering 
    from sklearn.cluster import AgglomerativeClustering 

    # get kwargs that are passed, if not set defaults
    affinity = kwargs.get('affinity', 'euclidean')
    linkage = kwargs.get('linkage', 'single')
    n_clusters = kwargs.get('n_clusters', None)
    connectivity = kwargs.get('connectivity', None)
    compute_full_tree = kwargs.get('compute_full_tree', True)
    distance_threshold = kwargs.get('distance_threshold', 0.25)
    compute_distances = kwargs.get('compute_distances', False)

    if len(X) > 1:
        clustering = AgglomerativeClustering(compute_full_tree=compute_full_tree, connectivity=connectivity, compute_distances=compute_distances, linkage=linkage, distance_threshold = distance_threshold, n_clusters=n_clusters, affinity=affinity).fit(X)
        
        #* Add the clustering info to rag
        for i in rag2.nodes():
            rag2.nodes[i]['cluster'] = str(clustering.labels_[i-1])
        if return_cluster_obj:
            return (rag2, clustering)
        else:
            return rag2
    elif len(X) == 1:
        # rag2.nodes[1]['cluster'] = '1'
        nx.set_node_attributes(rag2, '1', 'cluster')
        print('Only one node in the graph')
        return rag2
    else: 
        raise ValueError('No PSD with PBC data in the RAG (pore graph). Call add_pbc_psd_to_rag() first.')

def add_clustering_to_rag_pca(rag2,  return_cluster_obj = False, no_pbc =False, use_kpca=True, kwargs_cluster={}, kwargs_pca={}):
    """
    > We take the RAG, get the PSD data for each node, and then use agglomerative clustering to cluster
    the nodes.
    
    :param rag2: the RAG object
    :param distance_threshold: The linkage distance at which to cut the dendrogram to determine the
    clusters
    :param return_cluster_obj: if True, returns the clustering object as well as the rag, defaults to
    False (optional)
    :return: The rag2 object with the cluster labels added to the nodes.
    """
    
    import numpy as np 
    import networkx as nx 

    # *  get the data array, while reordering to match the node order   
    if no_pbc:
        X = np.vstack(dict(sorted(nx.get_node_attributes(rag2, 'psd_y').items())).values())
    else:
        X = np.vstack(dict(sorted(nx.get_node_attributes(rag2, 'psd_y_pbc').items())).values())

    if use_kpca:
        # get kwargs that are passed, if not set defaults
        kernel = kwargs_pca.get('kernel', 'cosine')

        # do the kernal PCA on X 
        from sklearn.decomposition import KernelPCA
        kpc = KernelPCA(**kwargs_pca).fit(X)
        X = kpc.transform(X)

        #* if you cannot solve for eigen values use regular PCA
        # * usually happens when there is only one pocket in the framework
        if X.shape[1]==0:
            X = kpc.X_fit_



    # * do the clustering using agglomerative clustering 
    from sklearn.cluster import AgglomerativeClustering 

    # get kwargs that are passed, if not set defaults
    affinity = kwargs_cluster.get('affinity', 'euclidean')
    linkage = kwargs_cluster.get('linkage', 'single')
    n_clusters = kwargs_cluster.get('n_clusters', None)
    connectivity = kwargs_cluster.get('connectivity', None)
    compute_full_tree = kwargs_cluster.get('compute_full_tree', True)
    distance_threshold = kwargs_cluster.get('distance_threshold', 0.25)
    compute_distances = kwargs_cluster.get('compute_distances', False)

    if len(X) > 1:
        clustering = AgglomerativeClustering(compute_full_tree=compute_full_tree, connectivity=connectivity, compute_distances=compute_distances, linkage=linkage, distance_threshold = distance_threshold, n_clusters=n_clusters, affinity=affinity).fit(X)
        
        #* Add the clustering info to rag
        for i in rag2.nodes():
            rag2.nodes[i]['cluster'] = str(clustering.labels_[i-1])
        if return_cluster_obj:
            return (rag2, clustering)
        else:
            return rag2
    elif len(X) == 1:
        # rag2.nodes[1]['cluster'] = '1'
        nx.set_node_attributes(rag2, '1', 'cluster')
        print('Only one node in the graph')
        return rag2
    else: 
        raise ValueError('No PSD with PBC data in the RAG (pore graph). Call add_pbc_psd_to_rag() first.')


def compare_two_paths(p1, p2, rag):
    """
    It takes two paths and a region adjacency graph, and returns a list of tuples, where each tuple
    contains a node in the first path and a node in the second path
    
    :param p1: a list of tuples, each tuple is a coordinate in the path
    :param p2: the path to be compared to p1
    :param rag: the region adjacency graph
    """
    import networkx.algorithms.isomorphism as iso
    g1 = rag.subgraph(p1)
    g2 = rag.subgraph(p2)
    nm = iso.categorical_node_match('cluster', default='None')
    return iso.is_isomorphic(g1, g2, node_match=nm)

def find_isomorphic_path_groups(list_of_paths,rag, cluster=True):
    import itertools 
    import networkx as nx
    gg = nx.Graph()
    if cluster:
        rag = add_clustering_to_rag(rag)
    for i in range(len(list_of_paths)):
        gg.add_node(i)
    indices = list(itertools.combinations(range(len(list_of_paths)), 2))
    for i,p in enumerate(itertools.combinations(list_of_paths,2)):
        # /else:
        if compare_two_paths(p[0], p[1], rag):
            # print(p)
            gg.add_edge(indices[i][0], indices[i][1])
    # * get a list tof groups of paths that are isomorphic 
    ccgg = list(nx.connected_components(gg))
    grouped_paths = [[list_of_paths[i] for i in cg] for cg in ccgg]
    return grouped_paths

def calc_path_components(rag, path: list):
    """
    Calculate the PCA components of a path in a RAG.

    Parameters
    ----------
    rag : nx.Graph
        The RAG to calculate the PCA components for.
    path : list
        The path to calculate the PCA components for.

    Returns
    -------
    np.ndarray
        The PCA components of the path.
    """
    import numpy as np
    import networkx as nx
    from sklearn.decomposition import PCA 
    from sklearn.preprocessing import StandardScaler
    pca_cart = PCA(n_components=np.min([len(path)-1,3]))
    pca_frac= PCA(n_components=np.min([len(path)-1,3]))
    maxima_coordinates = np.vstack(nx.get_node_attributes(rag.subgraph(path), 'maxima_indices').values())
    # maxima_fractional = mgr.get_fractional_coordinates(maxima_cartesian, ase_atoms)
    maxima_coordinates= StandardScaler().fit_transform(maxima_coordinates)
    pca_cart.fit(maxima_coordinates)
    # pca_frac.fit(maxima_fractional)
    # flag = np.logical_and(pca_cart.explained_variance_ratio_[0] > threshold, pca_frac.explained_variance_ratio_[0] > threshold)
    # flag = pca_cart.explained_variance_ratio_[0] > threshold
    
    # print(path, pca.explained_variance_ratio_)
    # return (pca_cart.explained_variance_ratio_, pca_frac.explained_variance_ratio_, flag )
    return pca_cart.explained_variance_ratio_

def find_paths(regions, rag2):

    import numpy as np 
    #* Now we are going to go face by face, and find all the paths 
    #* that have source on starting face and destination on one of the other 5 faces
    import networkx as nx 
    all_paths = list(nx.all_pairs_bellman_ford_path(rag2, weight='weight'))
    # all_paths = list(nx.all_pairs_dijkstra_path(rag2, weight='weight'))
    # all_paths = list(nx.all_pairs_dijkstra_path(rag2, weight='weight'))
        #* make a subset of these paths to check later for channels
    check_these_paths =[]

    shape = regions.shape 
    map_list=[]
    for a  in [0,1,2]:
            llr = regions.take(indices=0 ,axis=a) #* left wall
            rrr = regions.take(indices=shape[a]-1 ,axis=a)  #* right wall
            map_list.append(np.unique(list(zip(llr.flatten(), rrr.flatten())),axis=0))
    # map_list = np.vstack(map_list)


    #* check in the forward direction
    for a in [0,1,2]: #* along the three faces (directions) of the cube
        for thing in all_paths: #* loop over the tuples in all_paths
            paths = list(thing[1].items())[1:] #* all paths starting from source
            for dest, path in paths: #* the first entry is the source node, second is a dictionary of paths
                if(len(path)>1):
                    # print(len(path))
                    for node_i,node in list(enumerate(path))[1:]:  
                        #* check in the forward direction
                        # source_nodes = np.unique(wall_list[a*2+0])
                        source_nodes = map_list[a][:,0]
                        # print(node_i, node, path, source_nodes)
                        #* only the starting node should be on the starting face
                        if np.isin(node, source_nodes):
                            break 
                        # destination_nodes  =   np.unique([v for i,w in enumerate(wall_list) for v in w  if i!=a*2+0])
                        # destination_nodes= np.unique(wall_list[a*2+1])
                        destination_nodes= map_list[a][:,1]
                        if np.logical_and(*[np.isin(thing[0],source_nodes),np.isin(node,destination_nodes)]): #* thing = (source, {dest1:path1,dest2:path2,... })
                            check_these_paths.append(path[:node_i+1]) #* grab that path up to that node
                            # * add a path only once, once a destination is reached, move on to the next path.
                            break
                #* Add the single node paths to the list of paths to check
                # * where a single node occupies the whole path
                # if (len(path)==1):
                #     source_nodes = np.unique(wall_list[a*2+0])
                #     if np.isin(path[0], source_nodes):
                #         check_these_paths.append(path)

    #* Check in the backward direction
    for a in [0,1,2]: #* along the three faces (directions) of the cube
        for thing in all_paths: #* loop over the tuples in all_paths
            paths = list(thing[1].items())[1:] #* all paths starting from source
            for dest, path in paths: #* the first entry is the source node, second is a dictionary of paths
                if(len(path)>1):
                    # print(len(path))
                    for node_i,node in list(enumerate(path))[1:]:  
                        #* check in the backward direction
                        # source_nodes = np.unique(wall_list[a*2+1])
                        source_nodes = map_list[a][:,1]
                        # print(node_i,node, path, source_nodes)
                        #* only the starting node should be on the starting face
                        if np.isin(node, source_nodes):
                            break 
                        # destination_nodes  =   np.unique([v for i,w in enumerate(wall_list) for v in w  if i!=a*2+1])
                        # destination_nodes= np.unique(wall_list[a*2+0])
                        destination_nodes= map_list[a][:,0]
                        if np.logical_and(*[np.isin(thing[0],source_nodes),np.isin(node,destination_nodes)]): #* thing = (source, {dest1:path1,dest2:path2,... })
                            check_these_paths.append(path[:node_i+1]) #* grab that path up to that node
                            # * add a path only once, once a destination is reached, move on to the next path.
                            break
                            #* Add the single node paths to the list of paths to check


    #region #* Remove paths that non unique, in terms of repetiion indepent of order (assuming that no relevant paths are cyclic)
    frozen_paths = [] 
    check_these_paths_new=[]
    for i, p in enumerate(check_these_paths):
        if frozenset(p) not in frozen_paths:
            frozen_paths.append(frozenset(p))
            check_these_paths_new.append(p)
    checks = [calc_path_components(rag2,p) for p in check_these_paths_new]
    return checks, check_these_paths_new

def get_connections_for_rag(rag,region_labels, maxima,ase_atoms,dgrid, minimum_window_separation=1.0, wall_windows=True, wall_thickness=2):

    

    def find_windows(regs):
        

            import numpy as np
            from sklearn.cluster import AgglomerativeClustering
            i = regs[0]
            j = regs[1]
            check_flag = np.logical_or(np.logical_and(region_labels== i, outer.get(j)), np.logical_and(region_labels== j, outer.get(i)))
            if np.sum(check_flag) > 0:
                # print("Edges found between "+str(i)+' and '+ str(j) + ' sum: '+str(np.sum(check_flag)))
                windices = np.vstack(np.where(check_flag)).T  # * (N,3)
                window_points = np.dot(A_unit, (windices /region_labels.shape).T).T  # * (N,3)
                # ic(window_points)
                # * No clustering
                # window_centers = np.mean(window_points, axis=0)
                
                # * DBScan
                # X_scaled = StandardScaler().fit_transform(window_points)
                # db = DBSCAN(eps=0.5, min_samples=10).fit(X_scaled)
                # window_centers = np.vstack([window_points[db.labels_==i].mean(axis=0) for i in np.unique(db.labels_) if i != -1])
                
                # if len(window_points) >1:
                # * Agglomerative
                ag = AgglomerativeClustering(n_clusters=None, linkage='single',distance_threshold=minimum_window_separation, compute_full_tree=True).fit(window_points)
                window_centers = np.vstack([window_points[ag.labels_==i].mean(axis=0) for i in np.unique(ag.labels_) if i != -1])
            # else:
                # window_centers = window_points
                window_fractional = get_fractional_coordinates(window_centers,ase_atoms)
                window_radii = np.array(dinterp(window_fractional))
                return np.array([regs, True,'internal', window_centers, window_radii], dtype='object') # returns a list of [[r1,r2], connected or not, array of window_centers]
            else:
                return np.array([regs, False,'internal', None, None], dtype='object')

    #* find the internal windows
    from skimage.segmentation import find_boundaries
    from .helper import interpolate_labels, get_fractional_coordinates
    import ase
    import numpy as np
    connected_regions = np.vstack(list(rag.edges()))
    list_of_region_indices = np.delete(np.unique(region_labels),0)
    outer = {reg:find_boundaries(region_labels== reg, mode='outer') for reg in list_of_region_indices}
    dinterp = interpolate_labels(dgrid)
    minimum_window_separation = 1.0
    A_unit = ase.geometry.complete_cell(ase_atoms.get_cell()).T    

    connections= [find_windows(regs) for regs in rag.edges]

    #* Windows on the wall are to be calculated (periodic windows)
    if wall_windows:
        wwl = get_wall_windows2(regions=region_labels, maxima=maxima, dgrid=dgrid, ase_atoms=ase_atoms, wall_thickness=wall_thickness)
        if len(wwl)>0:
            connections = connections+wwl.tolist()
    
    return connections


def get_sampling_paths_for_region(region_index, connections, dgrid, maxima,  shape, ase_atoms):
    """
    > It finds the paths from the maxima of a region to the windows of that region
    
    :param region_index: the index of the region you want to find the paths for
    :param connections: a list of tuples, each tuple contains the indices of the connected regions and
    the window coordinates
    :param dgrid: the distance grid
    :param maxima: the indices of the maxima of the region
    :param shape: shape of the grid
    :param ase_atoms: ASE atoms object
    :return: The paths are being returned.
    """

    import numpy as np
    from .helper import get_fractional_coordinates
    #* maxima of the regiuon 
    region_maxima_index = maxima[region_index-1]

    #* let's find the indices ofthe windows of that region
    connected_regions = np.vstack([c[0] for c in connections])
    region_window_coords = np.vstack(np.array(connections)[np.any(connected_regions==region_index, axis=1)][:,-2])
    region_window_frac_coords = get_fractional_coordinates(region_window_coords, ase_atoms)
    region_window_indices = (region_window_frac_coords*shape).astype(int)

    #* dijkstra path finding
    import dijkstra3d 

    # * use exp field this works fine.
    field = np.exp(-(dgrid - dgrid.min()))

    # * find the paths
    paths = []
    for wc in region_window_indices:
        paths.append(dijkstra3d.dijkstra(field, region_maxima_index, wc, connectivity=26))

    return paths


def sample_path(path_indices, points_per_path=100):
    """
    It takes a list of points and returns a list of points that are evenly spaced along the path defined
    by the input points using spline interpolation
    
    :param path_indices: a numpy array of shape (n,3) where n is the number of points in the path
    :param points_per_path: The number of points to sample along the path, defaults to 100 (optional)
    :return: a list of points that are evenly spaced along the path.
    """

    from scipy.interpolate import splprep, splev  
    import numpy as np

    tck, u = splprep([path_indices[:,0],path_indices[:,1],path_indices[:,2]], s=0)
    x  = np.linspace(0,1,points_per_path)
    new_points = np.array(splev(x, tck)).T
    return new_points



def get_distances_for_path(path_indices, shape, ase_atoms, n_neighbors=10, radii=None):
    """
    > For each grid point, find the nearest neighbors, and return the minimum distance to the nearest
    neighbor minus the covalent radius of the neighbor
    
    :param path_indices: the indices of the path
    :param shape: the shape of the grid
    :param ase_atoms: the ase atoms object for the framework
    :param n_neighbors: the number of neighbors to consider for each grid point, defaults to 10
    (optional)
    :param radii: the covalent radii of the atoms in the structure. If not provided, they will be
    calculated from the atomic numbers using ase.data.covalent_radii
    :return: The distances between the atoms and the grid points.
    """

    from .helper import find_coord_from_indices
    import numpy as np
    from ase.data import covalent_radii
    if radii is None:
        radii = covalent_radii[ase_atoms.get_atomic_numbers()]
    ab = AABB_on_atoms(ase_atoms)
    ab_r = ab.query(find_coord_from_indices(path_indices, shape, ase_atoms), query_args=dict(mode='nearest', num_neighbors=n_neighbors, exclude_ii=True))

    ab_r_nl = ab_r.toNeighborList(sort_by_distance=True) #convert query to neighborlist
    nls = ab_r_nl.point_indices.reshape(-1,n_neighbors)  # get the indices of the neighbors for each grid pt. to lookup the radii later
    dists =ab_r_nl.distances.reshape(-1,n_neighbors) # get the distances to the neighbors of each grid pt.
    values =  np.min(dists - [radii[nl] for nl in nls],axis=1)
    return values



def get_psd_for_region(region_index, connections, dgrid, maxima, shape, ase_atoms, n_neighbors=10, radii=None, points_per_path=1000, set_bins=50, return_bins=False):
    """
    This function takes in a region index, the connections between regions, the distance grid, the
    maxima, the shape of the distance grid, the ase atoms object, the number of neighbors to consider,
    the radii of the atoms, the number of points to sample along the path, the number of bins to use for
    the histogram, and a boolean to return the bins or not. 
    
    It returns the histogram values and the bins if return_bins is True, otherwise it returns just the
    histogram values.
    
    :param region_index: The index of the region for which you want to get the psd
    :param connections: list of tuples of the form (region_index, region_index)
    :param dgrid: distance grid
    :param maxima: the maxima of the distance map
    :param shape: shape of the grid
    :param ase_atoms: ASE atoms object
    :param n_neighbors: number of neighbors to consider for each point in the path, defaults to 10
    (optional)
    :param radii: list of radii for which you want to calculate the psd. If you want to calculate the
    psd for a single radius, then you can pass a single value
    :param points_per_path: number of points to sample along the path, defaults to 1000 (optional)
    :param set_bins: The number of bins to use for the histogram, defaults to 50 (optional)
    :param return_bins: If True, returns the bin centers as well as the histogram values, defaults to
    False (optional)
    :return: the psd for a given region.
    """

    import numpy as np
    #* get the paths 
    connected_regions = np.vstack([c[0] for c in connections])
    if np.any(connected_regions == region_index):


        #* get the paths 
        paths = get_sampling_paths_for_region(region_index, connections, dgrid, maxima, shape, ase_atoms)


        #* resample the paths to make them finer
        paths = [sample_path(path, points_per_path=points_per_path) for path in paths]

        #* get the distances for each path
        distances = []
        for path in paths:
            distances.append(get_distances_for_path(path, dgrid.shape, ase_atoms, n_neighbors=n_neighbors,radii=radii))

        #* get the psd for all the paths
        vals, bins = np.histogram(np.concatenate(distances), density=True, bins=set_bins)
        # binc = (bins[1:]+bins[:-1])/2
        if return_bins:
            binc = (bins[1:]+bins[:-1])/2
            return vals, binc
        else:
            return vals
    else:
        print('The region {} is not connected to any other region. Hence unable to detect cavity axis for sampling'.format(region_index))
        return 'NaN', 'NaN'

def apply_pbc(regions, maxima, maxima_radii, ase_atoms, wall_thickness=1, return_conn =False, minimum_overlap_fraction=0.05):
    """
    It takes a 3D array of regions and returns a 3D array of regions with periodic boundary conditions
    applied
    
    :param regions: the regions array
    :param maxima: the list of maxima
    :param maxima_radii: the radius of each region
    :param ase_atoms: The atoms object from ASE
    :param wall_thickness: The thickness of the wall in Angstroms, defaults to 1 (optional)
    :param return_conn: if True, returns a list of lists of connected regions, defaults to False
    (optional)
    :param minimum_overlap_fraction: This is the minimum fraction of overlap at the wall region between two regions for
    them to be considered for connected via PBC. Then we check the maxima locations also to make sure, they are not separated by periodic windows.
    :return: a list of connected components. Each connected component is a list of regions a.k.a the PBC group.
    """



    import copy
    # a = 0 
    import numpy as np
    from icecream import ic
    ic.disable()
    regions_pbc = copy.deepcopy(regions)
    shape = regions.shape
    import networkx as nx
    g1= nx.Graph()
    g1.add_nodes_from(range(len(np.unique(regions))))



    #* we need to recompute the wall thickness for non orthorhombic cells 
    #* for this we shall use the perpendicualr length , which is the norm of the reciprocal lattice vector in each direction
    from ase import geometry


    #* let's compute the perpendicular length of each voxel
    #* we will use the reciprocal lattice vectors to do this
    #* we need this to determine how many voxels we need to select as the wall of the unit cell in each direction
    la = ase_atoms.get_cell_lengths_and_angles()[0]/regions.shape[0]
    lb = ase_atoms.get_cell_lengths_and_angles()[1]/regions.shape[1]
    lc = ase_atoms.get_cell_lengths_and_angles()[2]/regions.shape[2]
    alpha = ase_atoms.get_cell_lengths_and_angles()[3] * (np.pi/180.0)
    beta = ase_atoms.get_cell_lengths_and_angles()[4] * (np.pi/180.0)
    gamma = ase_atoms.get_cell_lengths_and_angles()[5] * (np.pi/180.0)
    vol = ase_atoms.get_volume()/regions.size
    eA = [la, 0, 0]
    eB = [lb*np.cos(gamma), lb*np.sin(gamma), 0]
    eC = [lc*np.cos(beta), lc*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma), vol/(la*lb*np.sin(gamma))]

    la_perp =  vol/np.linalg.norm(np.cross(eB, eC))
    lb_perp =  vol/np.linalg.norm(np.cross(eC, eA))
    lc_perp =  vol/np.linalg.norm(np.cross(eA, eB))


    wall_cells = np.ceil([wall_thickness/la_perp, wall_thickness/lb_perp, wall_thickness/lc_perp]).astype(int)

    for a in [0,1,2]:
        
        # ic(a)
        # print('a', a)

        #* Here we are comparing the left wall and the right wall, so to get th ecorrespoding regions, 
        #* we have to mirror the right wall
        llr = regions_pbc.take(indices=np.arange(0, wall_cells[a]) ,axis=a) #* left wall
        rrr = regions_pbc.take(indices=np.arange(shape[a]-1, shape[a]-1-wall_cells[a], -1) ,axis=a)  #* right wall
        llf = llr.flatten()
        rrf = rrr.flatten()
        # ic(dict(zip(llr.flatten(), rrr.flatten())))
        # ic(dict(zip(rrr.flatten(),llr.flatten())))
        # dlist = [dict(sorted(zip(llr.flatten(), rrr.flatten()))),dict(sorted(zip(rrr.flatten(),llr.flatten())))]

        #* Check left vs right wall
        # map_dict = dlist[0]
        # map_dict.pop(0)
        # ic(map_dict)


        map_list = np.unique(list(zip(llr.flatten(), rrr.flatten())),axis=0)
        map_list = map_list[np.all(map_list!=0,axis=1)]
        ic(map_list)
        # print('map_list', map_list)
        if len(map_list)> 0:

            # * forwards check
            reg_numbers = np.vstack([(k, v) for k,v in map_list])


            #* let's compute and prnt percentage overlap between pairs of mapping regions
            #* let's filter those connections with less thant 10% overlap
            sumlist = [np.sum(llf==r[0]) for r in reg_numbers]
            matching = [np.sum(np.logical_and(llf==r[0], rrf==r[1])) for r in reg_numbers] 
            fractional_overlap = [m/s for m,s in zip(matching, sumlist)]

            #* flags where there should be a left to right connection 
            left_right = [f> minimum_overlap_fraction for f in fractional_overlap]

            
            # reg_numbers = reg_numbers[np.where([f > minimum_overlap_fraction for f in fractional_overlap])]
            # print('left-right\n')
            # print(list(zip(reg_numbers, fractional_overlap)))

            # * backwards check

            # #* let's compute and prnt percentage overlap between pairs of mapping regions
            # #* let's filter those connections with less thant minimum overlap fraction
            sumlist = [np.sum(rrf==r[1]) for r in reg_numbers]
            matching = [np.sum(np.logical_and(llf==r[0], rrf==r[1])) for r in reg_numbers] 
            fractional_overlap = [m/s for m,s in zip(matching, sumlist)]
            
            #* flags where there should be a right to left connection
            right_left = [f> minimum_overlap_fraction for f in fractional_overlap]
            
            # print('right-left\n')
            # print(list(zip(reg_numbers, fractional_overlap)))
            # now let's combine the two flags
            pick_these = np.logical_or(left_right, right_left)

            reg_numbers = reg_numbers[pick_these]
            # print(pick_these)
            # print()
            # print(reg_numbers, 'reg numbers after picking')
            if len(reg_numbers)>0:
                flags = np.vstack([[maxima[r[0]-1][a]<wall_cells[a], maxima[r[1]-1][a]>(shape[a]-wall_cells[a])] for r in reg_numbers])
                # flags = np.logical_or(flags[:,0],flags[:,1,])
                # flags = np.logical_or(flllags[:,0],flags[:,1,])
                flags1 = np.logical_and(flags[:,0], flags[:,1])
                # print('flags1\n',list(zip(reg_numbers,flags1)))
                
                flags2 = np.logical_xor(flags[:,0], flags[:,1])
                # print('flags2\n',list(zip(reg_numbers,flags2)))
                reg_radii = np.vstack([[maxima_radii[r[0]-1], maxima_radii[r[1]-1]] for r in reg_numbers])
                flags3 = np.vstack([reg_radii[i][np.where(f)]<reg_radii[i][np.where(np.invert(f))][0] if flags2[i] else False for i,f in enumerate(flags)])
                # print('flags3\n',list(zip(reg_numbers,flags3)))
                # flags3 = np.vstack([ True if flags2[i] else False for i,f,r  in enumerate(list(zip(flags, reg_numbers)))])
                # reg_radii = [maxima_radii[r[0]-1], maxima_radii[r[1]-1] for r in reg_numbers]
                
                flags = [flags1[i] if flags1[i] else flags2[i] if flags3[i] else False for i,f in enumerate(flags)]
                # print('after checks\n',list(zip(reg_numbers,flags)))
                # ic(flags1)    
                # ic(flags2)    
                # ic(flags3)    
                # ic(flags)       # flags = [True if flags1[i] else True if np.logical_and(flags2[i], flags3[i] ) else False for i,f in enumerate(flags)]
                
                #* filter based on region numbers
                reg_numbers = reg_numbers[np.where(flags)]
                
                # g1.add_edges_from(reg_numbers[np.where(flags)])
                g1.add_edges_from(reg_numbers)


    conn = nx.connected_components(g1)
    # ic(list(conn))
    # ic(list(G.edges()))
    # ic([c for c in list(conn)])
    for i,c in enumerate(conn):
        for r in c:
            regions_pbc[np.where(regions==r)] = i
            ic('replace {0} with {1}'.format(r, i))

    # if reindex:
    #     # * Re-index the regions so that it runs from zero to N-regions.
    #     old_labels = np.unique(regions_pbc)  # This is sorted aalreadylready
    #     number_of_regions = len(old_labels)
    #     new_labels = range(number_of_regions)

    #     # * Replace the old-labels with new
    #     for i in range(number_of_regions):
    #         regions_pbc[regions_pbc== old_labels[i]] = new_labels[i]
    if return_conn:
        return list(nx.connected_components(g1))[1:]
    else:
        return regions_pbc

def add_pbc_to_rag(rag, pbc_groups):
    """
    It takes a RAG and a list of PBC groups, and adds a 'pbc_group' attribute to each node in the RAG

    :param rag: the region adjacency graph
    :param pbc_groups: a list of lists of labels that are in the same periodic boundary condition group
    :return: A rag with the pbc_group attribute added to each node.
    """
    import numpy as np
    for node in rag.nodes:
        pbc_index = np.where([rag.nodes[node]['labels'][0] in list(g) for g in pbc_groups] )[0][0]
        rag.nodes[node]['pbc_group'] = pbc_index
    return rag
#region  # * Energy grid functions
#######################################################
############################
# * Let's put the energy grid functions in this section
############################


def get_repetitions_to_meet_cutoff(ase_atoms, cutoff=12.5):
    """
    It takes an ASE atoms object and a cutoff distance, and returns the number of unit cells in each
    direction that are needed to ensure that the cutoff distance is satisfied
    
    :param ase_atoms: ASE atoms object
    :param cutoff: the cutoff distance for the pairwise interactions
    :return: The number of repetitions in each direction to meet the cutoff.
    """

    import numpy as np
    la, lb, lc, alpha, beta, gamma = ase_atoms.get_cell_lengths_and_angles()
    alpha, beta, gamma = np.radians(alpha), np.radians(beta), np.radians(gamma)
    vol = ase_atoms.get_volume()
    eA = [la, 0, 0]
    eB = [lb*np.cos(gamma), lb*np.sin(gamma), 0]
    eC = [lc*np.cos(beta), lc*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma), vol/(la*lb*np.sin(gamma))]
    A=[eA, eB, eC]
    A = np.array(A)
    A = A.T
    A_inv = np.linalg.inv(A)
    lx_unit, ly_unit, lz_unit = vol/np.linalg.norm(np.cross(eB, eC)), vol/np.linalg.norm(np.cross(eC, eA)), vol/np.linalg.norm(np.cross(eA, eB))
    nx_cells, ny_cells, nz_cells = int(np.ceil(2*cutoff/lx_unit)), int(np.ceil(2*cutoff/ly_unit)), int(np.ceil(2*cutoff/lz_unit))
    return nx_cells, ny_cells, nz_cells


def repeat_to_meet_cutoff(ase_atoms, forcefield_mixed, cutoff=12.5):
    """
    It takes an ASE atoms object and a forcefield, and returns a new ASE atoms object and a new
    forcefield that are repeated in such a way that the vdw cutoff is met
    
    :param ase_atoms: the atoms object from ASE
    :param forcefield_mixed: the forcefield parameters for the mixed system
    :param cutoff: the cutoff distance for the forcefield
    :return: ase_atoms_repeat is a supercell of ase_atoms, forcefield_mixed_repeat is a tiled array of
    forcefield_mixed
    """


    # * let's find the repetitions needed to meet cutoff
    import numpy as np
    repeats =get_repetitions_to_meet_cutoff(ase_atoms, cutoff=cutoff)
    nrepeats = np.prod(repeats)

    
    # * make the P matrix for the supercell
    P = np.eye(3)
    from ase.build import make_supercell 
    nrepeats = np.prod(repeats)
    np.fill_diagonal(P, val=list(repeats))


    # * let's repeat everything
    ase_atoms_repeat = make_supercell(ase_atoms, P)

    forcefield_mixed_repeat = np.tile(forcefield_mixed, (nrepeats,1))
    return ase_atoms_repeat, forcefield_mixed_repeat

def gpe_aabb2(grid_points,positions, cell,forcefield_mixed,cutoff):

    import numpy as np

    from freud.locality import AABBQuery
    ab_r=AABBQuery(cell,positions).query(grid_points,query_args=dict(mode='ball',r_max=cutoff,exclude_ii=True))
    # ab_r=aabb_object.query(grid_points,query_args=dict(mode='ball',r_max=cutoff,exclude_ii=True))
    ab_r_nl=ab_r.toNeighborList(sort_by_distance=True)#convertquerytoneighborlist
    #nls=ab_r_nl.point_indices.reshape(-1,n_neighbors)#gettheindicesoftheneighborsforeachgridpt.tolookuptheradiilater
    #dists=ab_r_nl.distances.reshape(-1,n_neighbors)#getthedistancestotheneighborsofeachgridpt.
    #print(ab,ab_r,ab_r_nl,nls,[radii[nl]fornlinnls])#,np.min(dists-[radii[nl]fornlinnls],axis=1))
    #returntheminimum(shortestdistancetothesurface)ineachrow


    indices = ab_r_nl.point_indices # * frame indices within cutoff at each grid point
    distances = ab_r_nl.distances # * distances to frame indices at each grid point
    segments = ab_r_nl.segments # * number of frame indices at each grid point
    inter_params = np.take(forcefield_mixed, indices, axis=0) # * get the parameters for each frame index at each grid point
    sig  = inter_params[:,1] # * get the sigma values for each frame index at each grid point
    eps  = inter_params[:,0]    # * get the epsilon values for each frame index at each grid point
    sd = np.power(np.divide(sig, distances), 6) # * (sigma/distance)^6
    # sig = np.array([forcefield[i][1] for i in indices])
    # eps = np.array([forcefield[i][0] for i in indices])

    # print(len(indices),len(distances),len(segments),len(sig),len(eps)) #, np.take(forcefield, indices))
    v = np.add.reduceat(np.multiply(4*eps, ((sd)**2-(sd))), segments) # * sum the energy contributions for each grid point
    # E = np.add.reduceat(4*eps*((sig/distances)**12-(sig/distances)**6), segments, axis=0)


    # E = [lj_sum(d, s, e) for d, s, e in zip(distances, sig, eps)]

    # return np.vstack(E).reshape(-1,1)
    return v


def gpe_aabb(grid_points,aabb_object,forcefield_mixed,cutoff):
    """
    It takes a list of grid points,an AABB query object from freud, a list of forcefield parameters, and a cutoff
    distance, and returns the energy of each grid point
    
    :param grid_points: the points you want to evaluate the energy at
    :param aabb_object: aabb object from the aabbquery class from freud
    :param forcefield_mixed: a numpy array of shape (n,2) where n is the number of atoms in the system.
    The first column is the epsilon parameter and the second column is the sigma parameter
    :param cutoff: the cutoff distance for the neighbor search
    :return: The energy of each grid point.
    """

    # def lj_sum(dists, sig, eps):
    #     return np.sum(4*eps*((sig/dists)**12-(sig/dists)**6))

    #Performthenearestneighborquerywithn-neighbors(10bydefault)


    import numpy as np
    ab_r=aabb_object.query(grid_points,query_args=dict(mode='ball',r_max=cutoff,exclude_ii=True))


    ab_r_nl=ab_r.toNeighborList(sort_by_distance=True)#convertquerytoneighborlist
    #nls=ab_r_nl.point_indices.reshape(-1,n_neighbors)#gettheindicesoftheneighborsforeachgridpt.tolookuptheradiilater
    #dists=ab_r_nl.distances.reshape(-1,n_neighbors)#getthedistancestotheneighborsofeachgridpt.
    #print(ab,ab_r,ab_r_nl,nls,[radii[nl]fornlinnls])#,np.min(dists-[radii[nl]fornlinnls],axis=1))
    #returntheminimum(shortestdistancetothesurface)ineachrow
    indices = ab_r_nl.point_indices
    distances = ab_r_nl.distances
    segments = ab_r_nl.segments
    inter_params = np.take(forcefield_mixed, indices, axis=0)
    sig  = inter_params[:,1]
    eps  = inter_params[:,0]
    sd = np.power(np.divide(sig, distances), 6)
    # sig = np.array([forcefield[i][1] for i in indices])
    # eps = np.array([forcefield[i][0] for i in indices])

    # print(len(indices),len(distances),len(segments),len(sig),len(eps)) #, np.take(forcefield, indices))
    v = np.add.reduceat(np.multiply(4*eps, ((sd)**2-(sd))), segments)
    # E = np.add.reduceat(4*eps*((sig/distances)**12-(sig/distances)**6), segments, axis=0)


    # E = [lj_sum(d, s, e) for d, s, e in zip(distances, sig, eps)]

    # return np.vstack(E).reshape(-1,1)
    return v

# def aabb_for_atoms(ase_atoms):

#     from freud.locality import AABBQuery
#     from freud.box import Box
#     import ase
#     #fromase.dataimportvdw_radii
#     import pandas as pd
#     import numpy as np

#     framework_atom_positions=ase_atoms.get_positions()
#     cell=ase.geometry.complete_cell(ase_atoms.get_cell()).T#thisneedstobeanuppertriangularmatrix,essentialcolumnvectors,notrowvectors

#     box=Box.from_matrix(cell,dimensions=3)#DefinetheboxforFreud
#     aabb=AABBQuery(box=box,points=framework_atom_positions)#ComputetheAABBTree

#     return aabb

def egrid_from_atoms2(ase_atoms, spacing=0.5, blocksize=20000, forcefield_mixed= None, cutoff=12.5, probe_symbol='He', return_ncells=False, precompute=True):
    """
    > It takes an ASE atoms object, and returns a dask array of the vdw energy on a grid. 
    
    The function is a bit long, but it's not too complicated. 
    
    The first part of the function is just to get the grid points. 
    
    The second part is to get the forcefield parameters for the atoms. 
    
    The third part is to get the positions of the atoms in the unit cell, and to get the unit cell
    vectors. 
    
    The fourth part is to repeat the unit cell to get the positions of the atoms in the repeated unit
    cell. 
    
    The fifth part is to map the grid point energy function onto the grid points. 
    
    The sixth part is to reshape the grid points into a 3D array. 
    
    The seventh part is to return the grid points. 
    
    The eighth part is to return the number of repetitions of the unit cell if that return_ncells parameter is turned on.
    
    :param ase_atoms: an ase atoms object
    :param spacing: the spacing of the grid points
    :param blocksize: the size of the chunks that the grid is split into. This is important for memory
    management, defaults to 20000 (optional)
    :param forcefield_mixed: This is a list of the forcefield parameters for each atom in the unit cell
    :param cutoff: the cutoff distance for the Lennard-Jones potential
    :param probe_symbol: the symbol of the probe atom used to calculate the van der Waals radii,
    defaults to He (optional)
    :param return_ncells: if True, returns the number of cells needed to meet the cutoff, defaults to
    False (optional)
    :param precompute: if True, then the AABB is precomputed and passed to the function. If False, then
    the positions and cell are passed to the function, defaults to True (optional)
    :return: The grid2 is a dask array, which is a lazy array. It is not computed until you call
    compute() on it.
    """
    
    import dask.array as da
    import numpy as np
    from .data import mix_lorentz_berthelot
    from .data import forcefields
    import ase 
    # import mofography as mgr 
    grid, shape = dask_grid_over_atoms(ase_atoms, spacing=spacing, chunksize=blocksize)
    
    # ab = aabb_for_atoms(ase_atoms)
    if forcefield_mixed is None:
        forcefield_mixed = mix_lorentz_berthelot(ase_atoms, forcefield = forcefields['uff'], probe_symbol=probe_symbol)
    # grid.rechunk(chunks=(10000,3))


    # * Now here use the the repeated box to get the framework atom positions
    # * also the grid point coordinates now only span one unit cell, since others are simply copies

    from ase.build import make_supercell 
    P = np.eye(3)
    repeats =get_repetitions_to_meet_cutoff(ase_atoms, cutoff=cutoff)
    nrepeats = np.prod(repeats)
    np.fill_diagonal(P, val=list(repeats))


    # * let's repeat everything
    ase_atoms_repeat = make_supercell(ase_atoms, P)
    positions_repeat =ase_atoms_repeat.get_positions()
    cell_repeat = ase.geometry.complete_cell(ase_atoms_repeat.get_cell()).T
    forcefield_mixed_repeat  = np.tile(forcefield_mixed, (nrepeats,1))
    
    # mix_lorentz_berthelot(ase_atoms_repeat, forcefield = forcefields['uff'], probe_symbol=probe_symbol)
    if precompute:
        ab_repeat = AABB_on_atoms(ase_atoms_repeat)
        # * map the grid point energy function that uses position and cell, on to the grid points
        grid2 = da.map_blocks(gpe_aabb, grid, ab_repeat, forcefield_mixed_repeat, cutoff=cutoff, dtype=np.float64,drop_axis=1)
        grid2 = grid2.reshape(shape)
        grid2 = da.from_array(grid2.compute(scheduler='threads')).rechunk(chunks=(10,10,10)).reshape(shape)
    
    else:
        # * map the grid point energy function that uses position and cell, on to the grid points
        grid2 = da.map_blocks(gpe_aabb2, grid,positions_repeat, cell_repeat, forcefield_mixed_repeat, cutoff=cutoff, dtype=np.float64,drop_axis=1)
        grid2 = grid2.reshape(shape)

    if return_ncells:
        return grid2, get_repetitions_to_meet_cutoff(ase_atoms, cutoff=cutoff)
    else:
        return grid2

def get_energies_for_group(group, props, ase_atoms, forcefield_mixed, shape, cutoff=12.5, return_coords=False):
    
    from .helper import find_coord_from_indices
    import numpy as np
    # * repeat the ase atoms and forcefield to meet cutoff    
    # nrepeats = np.prod(mgr.get_repetitions_to_meet_cutoff(ase_atoms, cutoff))
    ase_atoms_repeat, forcefield_mixed_repeat = repeat_to_meet_cutoff(ase_atoms, forcefield_mixed, cutoff)
    # import numpy as np  

    # * get aabb object on repeated system
    ab_repeat = AABB_on_atoms(ase_atoms_repeat)
    # nrepeats = np.prod(mgr.get_repetitions_to_meet_cutoff(ase))

    # ase_atom_repeat = mgr.get_aseatoms_to_meet_cutoff(ase_atoms, cutoff)
    # forcefield_mixed_repeat = np.tile()
    # * get region indices for the regions in the group
    prop_labels = np.array([prop['label'] for prop in props])
    prop_indices = [np.where(prop_labels==g)[0][0] for g in group]
    
    # * this group has more than one entry, assuming it is a pbc group
    if len(group)>1:

        #* grid indices of all the points present in this group
        coords_indices = np.vstack([props[pi]['coords'] for pi in prop_indices])


        # * get the actual grid point coordinates from indices for the pbc_group
        coords =find_coord_from_indices(coords_indices, shape, ase_atoms) 

        # * get the grid point energies for the points in the group
        energies = gpe_aabb(coords, ab_repeat, forcefield_mixed_repeat, cutoff)

    else:

        pi = np.where(prop_labels==list( group )[0])[0][0]
        # print(pi)
        #* get coordinates and centroids
        coords_indices = props[pi]['coords']
        coords = find_coord_from_indices(coords_indices, shape, ase_atoms)
        # centroid = np.mean(coords, axis=0)
        # * histogram of the distances and bin centers
        energies = gpe_aabb(coords, ab_repeat, forcefield_mixed_repeat, cutoff)

    if return_coords:
        return energies, coords
    else:
        return energies

def add_vdw_hist_to_rag_with_pbc(rag, props, pbc_groups, ase_atoms, forcefield_mixed, shape, bins_energy=None , cutoff=12.5, minimum_npockets_per_type=3, add_energies=False):

    import numpy as np 
    import networkx as nx
    if 'cluster' not in rag.nodes[1].keys():
        raise ValueError('RAG must have a cluster attribute. Follow examples in docs to add clustering to RAG')
    else:
        types_list = [rag.nodes[list(g)[0]]['cluster'] for g in pbc_groups]

        if bins_energy is None:
            bins_energy = np.linspace(-3000, 100, 100)
        
        # list(nx.get_node_attributes(rag2, 'cluster').values())
        ctypes, cindices = np.unique(types_list,return_inverse=True)
        # for c in ctypes:
        types_indices = [np.where([t==c for t in types_list])[0] for c in ctypes]
        types_counts = [len(t) for t in types_indices]
        type_selection = [np.random.choice(ti , minimum_npockets_per_type, replace=False) if types_counts[i]>minimum_npockets_per_type else ti  for i,ti in enumerate(types_indices)]

        import itertools  
        merged_type_selection = [list(itertools.chain(*[list(pbc_groups[t]) for t in ts])) for ts in type_selection]
        energies_list = [get_energies_for_group(mts, props, ase_atoms, forcefield_mixed, shape) for mts in merged_type_selection]
        histograms_bins= [np.histogram(e, bins=bins_energy, density=False) for e in energies_list]
        ehistograms = [hb[0] for hb in histograms_bins]
        ebincenters = [(hb[1][:-1]+hb[1][1:])/2 for hb in histograms_bins]
        ehists_data = [ehistograms[cindices[i]] for i,g in enumerate(pbc_groups) for gg in list(g)]
        ebins_data =  [ebincenters[cindices[i]] for i,g in enumerate(pbc_groups) for gg in list(g)] 


        for node in rag.nodes():
            rag.nodes[node]['vdw_hist_y_pbc']=ehistograms[cindices[rag.nodes[node]['pbc_group']]]
            rag.nodes[node]['vdw_hist_x_pbc']=ebincenters[cindices[rag.nodes[node]['pbc_group']]]
        # nx.set_node_attributes(rag, dict(zip(rag.nodes,ehists_data)), 'vdw_hist_y_pbc')
        # nx.set_node_attributes(rag, dict(zip(rag.nodes,ebins_data)), 'vdw_hist_x_pbc')
        if add_energies:
            for node in rag.nodes():
                rag.nodes[node]['energies']=energies_list[cindices[rag.nodes[node]['pbc_group']]]
            
            # energies_nodes = [energies_list[c] for c in cindices]
            # nx.set_node_attributes(rag, dict(zip(rag.nodes, energies_nodes)), 'vdw_energies_pbc')
        
    return rag

#endregion

# """

# Helper functions for mofography
# ██╗  ██╗███████╗██╗     ██████╗ ███████╗██████╗ 
# ██║  ██║██╔════╝██║     ██╔══██╗██╔════╝██╔══██╗
# ███████║█████╗  ██║     ██████╔╝█████╗  ██████╔╝
# ██╔══██║██╔══╝  ██║     ██╔═══╝ ██╔══╝  ██╔══██╗
# ██║  ██║███████╗███████╗██║     ███████╗██║  ██║
# ╚═╝  ╚═╝╚══════╝╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝
                                                

# """                         
 
def find_coord_from_indices(indices, shape, ase_atoms):
    """Get the cartesian coordinates from the indces of a grid mapped over an ase atoms object.

    :param indices: indices of interest, can be a list of a (N,3) numpy array, where N is number of points
    :type indices: list or np.array (N,3)
    :param shape: overall shape of the 3D grid
    :type shape: np.shape output, so tuple or array of shape (1,3)
    :param ase_atoms: ASE atoms object corressponding to the grid
    :type ase_atoms: ase.atoms object
    :return: cartesian coordinates of the grid points
    :rtype: numpy arry of size (N,3) where N is the number of grid points
    """
    import ase
    import numpy as np
    cell = ase.geometry.complete_cell(ase_atoms.get_cell()).T
    return np.dot(cell, (indices/shape).T).T


def interpolate_labels(regions):
    """Create a scipy.interpolator.RegularGridInterpolator, using 'nearest' method for the region labels (or any 3D numpy grid) 
    over a unit box ([0,1]) in all axes.

    :param regions: 3D numpy array of region labels
    :type regions: numpy.array 
    :return: a scipy RegularGridInterpolator object that works on the fractional coordinates of atoms
    :rtype: scipy.interpolate.RegularGridInterpolator
    """    
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator as RGI
    xx = np.linspace(0,  1,regions.shape[0])
    yy = np.linspace(0,  1,regions.shape[1])
    zz = np.linspace(0,  1,regions.shape[2])
    rl = RGI((xx, yy, zz), regions, method='nearest')

    return rl


def get_fractional_coordinates(points, ase_atoms):
    """Find the fractional coordinates for a bunch of cartesian corrdinates in a box corresponding to an ase.atoms object

    :param points: numpy array of shape (Npoints, 3)
    :type points: numpy.array
    :param ase_atoms: ASE atoms object for the box
    :type ase_atoms: ase.atoms
    :return: numpy array of fractional coordinates of shape (Npoints,3)
    :rtype: numpy.array
    """    
    import numpy as np
    import ase
    cell= ase.geometry.complete_cell(ase_atoms.get_cell()).T 
    cell_inv = np.linalg.inv(cell)
    frac_coords = np.dot(cell_inv, points.T).T
    return frac_coords

def zoom_grid(grid, zoom_factor):
    """
    It takes a 2D array and a zoom factor, and returns a 2D array that is the original array zoomed by
    the zoom factor
    
    :param grid: the grid to be zoomed
    :param zoom_factor: The zoom factor along the axes. If a float, zoom is the same for each axis. If
    an array, zoom should contain one value for each axis
    :return: A dask array
    """


    import numpy as np
    import dask.array as da
    from scipy.ndimage import zoom
    return da.from_array(zoom(grid, zoom_factor, mode = 'grid-wrap', order=1))

def find_psd_peaks(psd, **kwargs):


    from scipy.signal import find_peaks_cwt, cwt, ricker

    # * set default parameters
    height  = kwargs.get('height', 0.2)
    width = kwargs.get('width', 15)

    # print('height: {}, width: {}'.format(height, width))
    # * compute the continuous wavelet transform
    cwt_psd= cwt(psd, wavelet=ricker, widths=[width])[0]

    #* compute peaks
    peaks = find_peaks_cwt(psd, min_snr = height,widths=width)

    #* filter out peaks that are below the height threshold
    valid_peaks = [p for p in peaks if cwt_psd[p]>height ]

    return valid_peaks



def compare_psd_peaks(peak_values1, peak_values2):



    import numpy as np 

    # * sort by length of peaks
    peak_values_array = np.array([peak_values1, peak_values2])
    sorted_indices  = np.argsort([len(peak_values1), len(peak_values2)])
    peak_values_array = peak_values_array[sorted_indices]

    # * assign to right order
    x2 = peak_values_array[0] # shorter one
    x1 = peak_values_array[1] # longer one

    # Find the indices of the closest values in x1 for each value in x2
    indices = np.interp(x2, x1, np.arange(len(x1)))

    # Round the indices to the nearest integer
    indices = np.round(indices).astype(int)

    # Get the closest values in x1 for each value in x2
    closest_values = x1[indices]

    sorted_and_matched = np.array([x2, closest_values])

    reconstructed_list = [[],[]]
    for i in sorted_indices:
        reconstructed_list[i] = sorted_and_matched[i]
    # Print the closest values
    # print(closest_values)
    return reconstructed_list
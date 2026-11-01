def make_raspa_input_for_cif(folder_path, raspa_settings_dict, other_files_to_copy, templates_dir, template_name):
    """Render simulation.input from a Jinja template and copy files into a run folder."""
    
    
    import jinja2 as jj 
    import shutil
    import os 
    ej = jj.Environment(loader = jj.FileSystemLoader(templates_dir))
    tlt = ej.get_template(template_name)
    simulation_input_string = tlt.render(**raspa_settings_dict)
    with open(os.path.join(os.path.abspath(folder_path), 'simulation.input'), 'w') as f: 
        f.write(simulation_input_string)
    for f in other_files_to_copy:
        shutil.copy(os.path.abspath(f), os.path.abspath(folder_path))

def get_isotherm_df_for_raspa_folder(path_to_folder):
    """Parse RASPA output *.data files in a folder into a sorted isotherm DataFrame.

    :param path_to_folder: the path to the folder containing the raspa output files
    :return: A dictionary with the isotherm data
    """

    import numpy as np
    import os
    import glob
    import pandas as pd

    current_path = os.getcwd()
    os.chdir(path_to_folder)
    files = glob.glob('**/*.data', recursive=True)
    data = np.zeros((len(files), 8))
    for i in range(len(files)):
            file = files[i]
            # print(file)
            data[i,0] = float(file.split('/')[0]) # Pressure
            lines = open(file).readlines()
            if np.sum(['Average loading excess [molecules/unit cell]' in ll for ll in lines])>0:
                # * Read the average loading and isosteric heats
                pfl = np.where(['Partial fugacity' in ll for ll in lines])[0][0]
                mucl =  np.where(['Average loading excess [molecules/unit cell]' in ll for ll in lines])[0][0]
                ccpgl =  np.where(['Average loading excess [cm^3 (STP)/gr framework]' in ll for ll in lines])[0][0]
                ugl=  np.where(['Ug should be subtracted' in ll for ll in lines])[0][0]
                data[i, 1] = float(lines[pfl].split()[2])
                data[i, 2] = float(lines[mucl].split()[5])
                data[i, 3] = float(lines[mucl].split()[7])
                data[i, 4] = float(lines[ccpgl].split()[6])
                data[i, 5] = float(lines[ccpgl].split()[8])
                data[i, 6] = float(lines[ugl-1].split()[0])
                data[i, 7] = float(lines[ugl-1].split()[2])
            else:
                data[i,:] = np.NaN*np.ones(np.size(data[i,:]))
    data = data[np.argsort(data[:, 0])]
    df = pd.DataFrame(data,
                      columns=['Pressure', 'Partial fugacity', 'Average loading excess [molecules/unit cell]', 'Error',
                               'Average loading excess [cm^3 (STP)/gr framework]', 'Error', 'Isos. Heat (kJ/mol)',
                               'Error'])
    # np.savetxt(filename,out,fmt='%.4f',delimiter='\t\t')
    os.chdir(current_path)
    return df
def get_isotherms_dict_for_multiple_raspa_folders(list_of_folder_paths):
    """Build a dict of isotherm DataFrames keyed by folder basename for each run folder.
    :param list_of_folder_paths: list of paths to the folders containing the raspa output files
    :return: a dictionary with the isotherm data
    """
    import numpy as np
    from tqdm import tqdm
    isotherms_dict = {}
    for folder_path in tqdm(list_of_folder_paths):
        df = get_isotherm_df_for_raspa_folder(folder_path)
        isotherms_dict[folder_path.split('/')[-1]] = df
    return isotherms_dict



#TODO: make a function to do clustered regions grid for a give output folder

# TODO: make a function to decompose the isotherm into component isotherms

# TODO: make a function to do the isotherm decomposition for a given output folder that wraps around all the above

def get_isos_fit_for_raspa_folder(path_to_folder, T, order):
    """Fit Legendre polynomials to RASPA energies and return isosteric heat vs loading.
    :type path_to_folder: string
    :param path_to_folder: Path to the RASPA folder with all the output_*.data files
    :type T: float
    :param T: Temperature
    :type order: integer
    :param order:Order of the Legendre polynomials to be used in regression
    :raises:
    :rtype: Number in box array, RASPA P.E array [K], Fitted P.E array[K], number in mol/kg, Isosteric heat (J/mol)
    """
    #* ###  #####  #######  #####  ####### ####### ######  ###  #####     #     # #######    #    #######
    #*  #  #     # #     # #     #    #    #       #     #  #  #     #    #     # #         # #      #
    #*  #  #       #     # #          #    #       #     #  #  #          #     # #        #   #     #
    #*  #   #####  #     #  #####     #    #####   ######   #  #          ####### #####   #     #    #
    #*  #        # #     #       #    #    #       #   #    #  #          #     # #       #######    #
    #*  #  #     # #     # #     #    #    #       #    #   #  #     #    #     # #       #     #    #
    #* ###  #####  #######  #####     #    ####### #     # ###  #####     #     # ####### #     #    #
    # * Python code to read the RASPA output and calulate the isosteric heat of adsorption
    # * by regressing <U> vs <N>
    # * This code uses Legendre polynomials to do the same
    # * Takes the path to the folder as the input
    # * Written and tested by Arun Gopalan, Snurr Group, Northwestern University @ 2019
    import numpy as np
    import glob
    import pandas as pd
    import scipy as sp
    import os
    import pandas as pd
    cwd = os.getcwd()
    # self.path= path_to_folder #'/projects/p20832/ALKANES_ZIFs/ZIF-4/flexible/hybrid/pentane_ArunGopalan2018ZIF-4_TraPPE_ZIF-4_mod_FLEX_EMPTY_10NVE_UNCHARGED_LOG_303'
    # self.T = T
    os.chdir(path_to_folder)
    datafiles = glob.glob('**/*.data', recursive=True)
    # * Read  all the required energies and the number of molecules in the box
    total_energy = []
    guest_guest_energy = []
    guest_host_energy = []
    number = []
    number_molpg = []
    for d in datafiles:
        filename = d
        f1 = open(filename).readlines()
        for i in range(len(f1)):
            if "Total energy:\n" in f1[i]:
                total_energy.append(float(f1[i+8].split()[1]))
            # if "Average Adsorbate-Adsorbate energy:\n" in f1[i]:
            #     guest_guest_energy.append(float(f1[i+8].split()[1]))
            # if "Average Host-Adsorbate energy:\n" in f1[i]:
            #     guest_host_energy.append(float(f1[i+8].split()[1]))
            # ! Reads the previous line which is the average
            if "Average loading absolute [molecules/unit cell]" in f1[i]:
                # print(f1[i-5:i+5])
                number.append(float(f1[i-1].split()[3]))
            if "Average loading absolute [mol/kg framework]" in f1[i]:
                number_molpg.append(float(f1[i].split()[5]))
    total_energy = np.array(total_energy)  # in kelvin per box
    number = np.array(number)        # per box
    # sort the energy and loading
    total_energy = total_energy[np.argsort(number)]
    number = np.sort(number)
    number_molpg = np.array(number_molpg)
    number_molpg = np.sort(number_molpg)
    R = 8.314  # J/(K. mol)

    # * Do the Legendre fit and derivative
    fit_coeff = np.polynomial.legendre.legfit(
        number, total_energy, order, full=False)
    energy_fit = np.polynomial.legendre.legval(number, fit_coeff)
    der_coeff = np.polynomial.legendre.legder(fit_coeff, m=1, scl=1)
    dUdN = np.polynomial.legendre.legval(number, der_coeff)
    Isos = R*(T - dUdN)  # in J/mol
    output = np.column_stack((number, total_energy, energy_fit, number_molpg, Isos))
    df = pd.DataFrame(output, columns=['Number', 'U_total', 'U_total_fit', 'Number_molpg', 'Isos'])
    os.chdir(cwd)
    return df

def get_isos_fit_for_multiple_raspa_folders(list_of_folder_paths, T, order):
    """Run get_isos_fit_for_raspa_folder on each folder and return a path-keyed dict.
    
    :param list_of_folder_paths: a list of strings, each string being the path to a folder
    containing RASPA output files
    :param T: temperature in K
    :param order: the order of the polynomial fit
    :return: A dictionary of the form {folder_path: isos_fit_for_raspa_folder(folder_path,
    T, order)}
    """
    from tqdm import tqdm
    isos_heats_dict = {}
    for folder in tqdm(list_of_folder_paths):
        isos_heats_dict[folder] = get_isos_fit_for_raspa_folder(folder, T, order)
    return isos_heats_dict



# def get_atoms_list_for_pdb(pdb_file, index_list=None):#, natoms, bonds_in_ref=None, trace_name='adsorbates', pbc=[1, 1, 1]):
#     """
#     It reads a RASPA pdb file and returns a list of ASE atoms objects
    
#     :param pdb_file: the name of the pdb file
#     :param index_list: a list of indices of the configurations you want to extract from the pdb file. If
#     you want all of them, just leave it as None
#     :return: A list of atoms objects.
#     """

#     from ase import Atoms
#     import numpy as np
#     from .data import read_raspa_pdb
#     from rich.progress import track
#     config_data = read_raspa_pdb(pdb_file)
#     atoms_list = []
#     if index_list is None:
#         index_list = np.arange(len(config_data['coords']))
#     for index in track(index_list, description='Reading Configs'):
#         coords = config_data['coords'][index]
#         sym = config_data['symbols'][index]
#         sym = [str(c) for c in sym]
#         cell = config_data['cells'][index]
#         # nmolecules = len(coords) // natoms  # numbre of co2 molecules
#         # convert to ase atoms object
#         atoms = Atoms(symbols=sym, positions=coords, cell=cell, pbc=[1, 1, 1])
#         atoms_list.append(atoms)
#     return atoms_list
def get_atoms_list_for_pdb(pdb_file, index_list=None, symbol_map=None):
    """Read selected frames from a RASPA trajectory PDB into a list of ASE Atoms.
    
    :param pdb_file: the name of the pdb file
    :param index_list: a list of indices of the configurations you want to extract from the pdb file. If
    you want all of them, just leave it as None
    :param symbol_map: Optional dictionary mapping custom symbols to standard element symbols.
    :return: A list of atoms objects.
    """

    from ase import Atoms
    import numpy as np
    from .data import read_raspa_pdb
    from rich.progress import track

    config_data = read_raspa_pdb(pdb_file, symbol_map=symbol_map)
    atoms_list = []
    
    if index_list is None:
        index_list = np.arange(len(config_data['coords']))
    
    for index in track(index_list, description='Reading Configs'):
        coords = config_data['coords'][index]
        sym = config_data['symbols'][index]
        
        # Check if coords is empty or not
        if coords.size == 0:
            print(f"Warning: No coordinates found for index {index}. Skipping this configuration.")
            continue
        
        # Ensure symbols are in the correct format
        sym = [str(c) for c in sym]
        
        # Check if cell dimensions are available
        if index < len(config_data['cells']):
            cell = config_data['cells'][index]
        else:
            print(f"Warning: No cell dimensions found for index {index}. Skipping this configuration.")
            continue
        
        # Convert to ASE Atoms object
        try:
            atoms = Atoms(symbols=sym, positions=coords, cell=cell, pbc=[1, 1, 1])
            atoms_list.append(atoms)
        except Exception as e:
            print(f"Error creating Atoms object for index {index}: {e}")
    
    return atoms_list
def apply_region_map_to_raspa_pdb_wrapped(pdb_file, regions_cluster, unit_cell, natoms_per_molecule=1, index_list=None, symbol_map=None):

    """Count adsorbates per region after wrapping PDB coordinates into one unit cell.

    :param pdb_file: the pdb file you want to analyze
    :param regions_cluster: the region map for a single unit cell
    :param unit_cell: the unit cell matrix or ASE Atoms object containing the unit cell information
    :param natoms_per_molecule: number of atoms per molecule. This is currently only implemented for 1,
    defaults to 1 (optional)
    :param index_list: list of indices of atoms to use in the pdb file. If None, all atoms are used
    :return: A dictionary with the number of atoms in each cluster type and the total number of atoms.
    """
    import numpy as np
    from .helper import interpolate_labels
    from ase.geometry import wrap_positions
    from ase import Atoms

    atoms_list = get_atoms_list_for_pdb(pdb_file, index_list=index_list, symbol_map=symbol_map)
    unique_cluster_types = np.unique(regions_cluster)
    # remove background 0
    unique_cluster_types = unique_cluster_types[unique_cluster_types!=0]

    # * make a dictionary of cluster types to store their counts and means
    cluster_type_counts = {str(c-1):[] for c in unique_cluster_types}
    cluster_type_counts['total']= []

    # this is the interpolator function for a single unit cell
    region_interp = interpolate_labels(regions_cluster)

    if natoms_per_molecule != 1:
        print('ERROR: natoms_per_molecule is not 1, currently only implemented for 1')

    # Get the unit cell
    if isinstance(unit_cell, Atoms):
        cell = unit_cell.get_cell()
    else:
        cell = unit_cell

    # loop over the snapshots in the pdb
    for snap in atoms_list:
        # Get the positions from the snapshot
        positions = snap.get_positions()
        
        # Wrap the positions to the unit cell
        wrapped_positions = wrap_positions(positions, cell)
        
        # Convert to fractional coordinates
        scaled_positions = cell.scaled_positions(wrapped_positions)

        # shift any scaled positions if they are out of bounds
        scaled_positions = scaled_positions%1
        # Apply the region map
        cluster_indices = region_interp(scaled_positions)
        
        number_total = len(cluster_indices)
        for c in unique_cluster_types:
            cluster_type_counts[str(c-1)].append(len(np.where(cluster_indices == c)[0]))
        cluster_type_counts['total'].append(number_total)

    for c in unique_cluster_types:
        cluster_type_counts[str(c-1)+'_mean'] = np.mean(cluster_type_counts[str(c-1)])
        cluster_type_counts[str(c-1)+'_std'] = np.std(cluster_type_counts[str(c-1)])
    cluster_type_counts['total_mean'] = np.mean(cluster_type_counts['total'])

    return cluster_type_counts

# TODO: Implement for natoms > 1
def apply_region_map_to_raspa_pdb(pdb_file, regions_cluster_super, natoms_per_molecule=1, index_list=None):
    """
    Count adsorbates per segmentation region using a supercell region map on PDB frames.
    
    :param pdb_file: the pdb file you want to analyze
    :param regions_cluster_super: the supercell regions map
    :param natoms_per_molecule: number of atoms per molecule. This is currently only implemented for 1,
    defaults to 1 (optional)
    :param index_list: list of indices of atoms to use in the pdb file. If None, all atoms are used
    :return: A dictionary with the number of atoms in each cluster type and the total number of atoms.
    """

    import numpy as np
    from .helper import interpolate_labels
    atoms_list = get_atoms_list_for_pdb(pdb_file, index_list=index_list)
    unique_cluster_types = np.unique(regions_cluster_super)
    # remove background 0
    unique_cluster_types = unique_cluster_types[unique_cluster_types!=0]

    # * make a dictionary of cluster types to store their counts and means
    cluster_type_counts = {str(c-1):[] for c in unique_cluster_types}
    cluster_type_counts['total']= []

    # this is the interpolator function
    region_interp_super = interpolate_labels(regions_cluster_super)
    if natoms_per_molecule !=1:
        print('ERROR: natoms_per_molecule is not 1, currently only implemented for 1')
    # loop over the snapshots in the pdb
    for snap in atoms_list:
        cluster_indices = region_interp_super(snap.get_scaled_positions())
        number_total = len(cluster_indices)
        for c in unique_cluster_types:
            cluster_type_counts[str(c-1)].append(len(np.where(cluster_indices == c)[0]))
        cluster_type_counts['total'].append(number_total)
    for c in unique_cluster_types:
        cluster_type_counts[str(c-1)+'_mean'] = np.mean(cluster_type_counts[str(c-1)])
        cluster_type_counts[str(c-1)+'_std'] = np.std(cluster_type_counts[str(c-1)])
    cluster_type_counts['total_mean'] = np.mean(cluster_type_counts['total'])

    return cluster_type_counts


def block_pockets_at_random_and_write(file_path, rag, ase_atoms, random_seed=123):
    """
    Write RASPA blocking spheres for a random fraction of each RAG pocket type.
    
    :param file_path: the path to the file you want to write the blocking information to
    :param rag: the networkx graph object
    :param ase_atoms: ASE atoms object
    :param random_seed: the random seed to use for the random number generator, defaults to 123
    (optional)
    :return: The blocking_fractions and _ are being returned.
    """

    from .helper import get_fractional_coordinates
    import numpy as np
    import networkx as nx
    cluster_type_dict = nx.get_node_attributes(rag, 'cluster')
    centroid_dict = nx.get_node_attributes(rag, 'centroid')
    size_dict = nx.get_node_attributes(rag, 'maxima_radii')
    np.random.seed(random_seed)
    unique_pocket_types = np.unique(list(cluster_type_dict.values())).astype(int)
    nunique = len(unique_pocket_types)
    blocking_fractions = np.random.random(nunique) # what percentage of each cavity type to block
    blocking_file = open(file_path, 'w')

    # print(unique_pocket_types)
    for i,c in enumerate(unique_pocket_types):
        typec = [k for k,v in cluster_type_dict.items() if v == str(int(c))]
        pockets_to_block = np.random.choice(typec, np.max([int(blocking_fractions[i]*len(typec)),1]), replace=False)
        # print(typec, int(blocking_fractions[i]*len(typec)), pockets_to_block)
        pockets_to_block_centroids = np.vstack([centroid_dict[p] for p in pockets_to_block])
        pockets_to_block_sizes = np.vstack([size_dict[p] for p in pockets_to_block])
        scaled_centroids = get_fractional_coordinates(pockets_to_block_centroids, ase_atoms)
        scaled_centroids = np.round(scaled_centroids,2)
        for i,s in enumerate(scaled_centroids):
            blocking_file.write(str(s[0])+'\t'+str(s[1])+'\t'+str(s[2])+'\t'+str(np.round(pockets_to_block_sizes[i][0],2))+'\n')
    blocking_file.close()
    return blocking_fractions, _


def block_specific_pockets_and_write(file_path, rag, ase_atoms, pocket_ids_to_block):
    """Write RASPA blocking spheres for selected RAG pocket node IDs.
    
    :param file_path: the path to the file you want to write the blocking information to
    :param rag: the networkx graph object
    :param ase_atoms: ASE atoms object
    :param pocket_ids_to_block: list of pocket IDs to block
    """
    from .helper import get_fractional_coordinates
    import numpy as np
    import networkx as nx
    
    centroid_dict = nx.get_node_attributes(rag, 'centroid')
    size_dict = nx.get_node_attributes(rag, 'maxima_radii')
    
    # Get centroids and sizes for specified pockets
    pockets_to_block_centroids = np.vstack([centroid_dict[p] for p in pocket_ids_to_block])
    pockets_to_block_sizes = np.vstack([size_dict[p] for p in pocket_ids_to_block])
    
    # Convert to fractional coordinates
    scaled_centroids = get_fractional_coordinates(pockets_to_block_centroids, ase_atoms)
    scaled_centroids = np.round(scaled_centroids, 2)
    
    # Write to file
    with open(file_path, 'w') as blocking_file:
        for i, s in enumerate(scaled_centroids):
            blocking_file.write(f"{s[0]}\t{s[1]}\t{s[2]}\t{np.round(pockets_to_block_sizes[i][0], 2)}\n")


# --- Percolation-based blocking spheres (RASPA fractional + radius in Å) -----------------


def percolated_region_labels_winding(
    results,
    probe_radius,
    include_dangling=True,
    expand_components=True,
    graph=None,
    base_region_labels=None,
):
    """
    Segmentation region labels (ints) that belong to a percolating cluster at
    ``probe_radius`` (Å), using winding percolation on the pore graph.

    Parameters
    ----------
    results : dict
        Segmentation output (``dgrid``, ``regions``, ``connections``,
        ``maxima_info``, ...).
    probe_radius : float
        Threshold radius in Å (e.g. half kinetic diameter).
    include_dangling, expand_components : bool
        Passed to :func:`mofography.winding_percolation.percolation_all_axes_winding`.
    graph, base_region_labels : optional
        Precomputed graph / label map; if None they are built from ``results``.

    Returns
    -------
    set of int
        Region labels that are accessible (percolated) at this threshold.
    """
    from .winding_percolation import percolation_all_axes_winding

    _, _, cluster_metadata = percolation_all_axes_winding(
        results,
        threshold=probe_radius,
        graph=graph,
        base_region_labels=base_region_labels,
        include_dangling=include_dangling,
        expand_components=expand_components,
    )
    percolated = set()
    for meta in cluster_metadata.values():
        percolated.update(int(x) for x in meta.get("regions", []) if int(x) != 0)
    return percolated


def segmentation_spacing_from_results(ase_atoms, shape, results=None):
    """
    Segmentation grid spacing (Å), consistent with ``dask_grid_over_atoms``.

    **Primary:** ``results["input_params"]["spacing"]`` is set by
    :func:`mofography.rust_fun.do_periodic_segmentation_rust` and should match
    the spacing used to build ``dgrid`` / ``regions``.

    **Also:** ``results["kwargs"]["spacing"]`` on some older dict layouts.

    **Fallback:** mean spacing from ASE cell lengths and ``shape`` if those
    keys are missing (e.g. ad-hoc grids or very old caches).
    """
    import numpy as np

    if results is not None:
        ip = results.get("input_params") or {}
        if isinstance(ip, dict) and ip.get("spacing") is not None:
            return float(ip["spacing"])
        kw = results.get("kwargs") or {}
        if isinstance(kw, dict) and kw.get("spacing") is not None:
            return float(kw["spacing"])
    cell = np.asarray(ase_atoms.get_cell(), dtype=float).reshape(3, 3)
    a, b, c = (np.linalg.norm(cell[i]) for i in range(3))
    nx, ny, nz = (int(s) for s in shape)
    sx = a / max(nx - 1, 1)
    sy = b / max(ny - 1, 1)
    sz = c / max(nz - 1, 1)
    return float(np.mean([sx, sy, sz]))


_infer_spacing_from_results = segmentation_spacing_from_results


def _blocking_spheres_cart_centers_radii(spheres, ase_atoms):
    import numpy as np

    cell = np.asarray(ase_atoms.get_cell(), dtype=float).reshape(3, 3)
    cell_t = cell.T
    centers = []
    radii = []
    for row in spheres:
        sx, sy, sz, r = (float(row[0]), float(row[1]), float(row[2]), float(row[3]))
        centers.append(cell_t @ np.array([sx, sy, sz], dtype=float))
        radii.append(float(r))
    return np.asarray(centers, dtype=float), np.asarray(radii, dtype=float)


def _sphere_union_mic_raw(pts, centers, radii, ase_atoms, atol, chunk_size):
    import numpy as np
    from ase.geometry import get_distances

    n = pts.shape[0]
    covered_mic = np.zeros(n, dtype=bool)
    covered_raw = np.zeros(n, dtype=bool)
    for center, rad in zip(centers, radii):
        for start in range(0, n, chunk_size):
            end = min(start + chunk_size, n)
            chunk = pts[start:end]
            _, d_mic = get_distances(
                chunk, center.reshape(1, 3), cell=ase_atoms.cell, pbc=True
            )
            covered_mic[start:end] |= d_mic[:, 0] <= rad + atol
            d_raw = np.linalg.norm(chunk - center, axis=1)
            covered_raw[start:end] |= d_raw <= rad + atol
    return covered_mic, covered_raw


def _attach_legacy_overlap_keys(out: dict) -> None:
    """Add keys expected by older overlap helpers and scripts."""
    n = int(out.get("n_denominator_points", 0))
    res = out.get("resolution")
    dom = out.get("denominator_domain")
    proxy_vol = out.get("approx_voxel_volume_proxy_ang3")
    if dom == "accessible":
        out["n_accessible_voxels"] = n
        if res == "coarse":
            out["n_accessible_sample_points"] = n
        if res == "interface":
            out["n_accessible_sample_points"] = n
            out["n_sample_points_interface"] = n
        if res in {"coarse", "interface"} and proxy_vol is not None:
            out["approx_voxel_volume_coarse_ang3"] = float(proxy_vol)
        if res == "fine":
            out["spacing"] = out.get("fine_spacing")
            out["approx_voxel_volume_ang3"] = out.get("approx_voxel_volume_proxy_ang3")
    else:
        out["n_inaccessible_sample_points"] = n
        if res in {"coarse", "interface"} and proxy_vol is not None:
            out["approx_voxel_volume_coarse_ang3"] = float(proxy_vol)
    pl = out.get("per_label") or {}
    if pl:
        if dom == "accessible":
            out["per_accessible_label"] = pl
        else:
            out["per_inaccessible_label"] = pl


def blocking_spheres_overlap_leakage(
    results,
    ase_atoms,
    spheres,
    probe_radius,
    *,
    denominator_domain="accessible",
    resolution="coarse",
    sample_spacing=0.5,
    fine_spacing=None,
    accessible_dgrid_ge_probe=False,
    restrict_to_rag_adjacent_accessible=True,
    require_face_neighbor_inaccessible=True,
    per_label=False,
    atol=1e-4,
    chunk_size=65536,
):
    """
    How much blocking spheres intrude into the **other** void class (leakage).

    * ``denominator_domain="accessible"`` — measure overlap into percolating void
      (typical when spheres block **inaccessible** pockets, ``BlockPockets``).
    * ``denominator_domain="inaccessible"`` — measure overlap into non-percolating
      void (when spheres target **accessible** regions, ``InvertBlockPockets``).

    ``resolution``: ``"fine"`` (full segmentation grid), ``"coarse"`` (subsample
    at ``sample_spacing``), or ``"interface"`` (accessible **interface band** only;
    ``denominator_domain`` must be ``"accessible"``).

    Fine grid spacing defaults to :func:`segmentation_spacing_from_results`
    (usually ``results["input_params"]["spacing"]`` from periodic segmentation).
    """
    import logging

    import numpy as np
    from .compute import dask_grid_over_atoms

    logger = logging.getLogger(__name__)
    resolution = str(resolution).lower()
    denominator_domain = str(denominator_domain).lower()
    if resolution not in {"fine", "coarse", "interface"}:
        raise ValueError(
            f"resolution must be 'fine', 'coarse', or 'interface', got {resolution!r}"
        )
    if denominator_domain not in {"accessible", "inaccessible"}:
        raise ValueError(
            "denominator_domain must be 'accessible' or 'inaccessible', "
            f"got {denominator_domain!r}"
        )
    if resolution == "interface" and denominator_domain != "accessible":
        raise ValueError(
            "resolution='interface' only supports denominator_domain='accessible'."
        )

    regions = np.asarray(results["regions"])
    dgrid = np.asarray(results["dgrid"], dtype=float)
    shape_f = regions.shape
    vol_cell = float(ase_atoms.cell.volume)

    percolated = percolated_region_labels_winding(results, probe_radius)
    all_labels = {int(x) for x in np.unique(regions) if int(x) != 0}
    inaccessible = all_labels - percolated

    if denominator_domain == "accessible":
        denom_mask = np.isin(regions, list(percolated))
        if accessible_dgrid_ge_probe:
            denom_mask &= dgrid >= float(probe_radius)
    else:
        denom_mask = np.isin(regions, list(inaccessible))

    mask_acc_report = np.isin(regions, list(percolated))
    if accessible_dgrid_ge_probe:
        mask_acc_report &= dgrid >= float(probe_radius)
    n_fine_accessible_ctx = int(np.count_nonzero(mask_acc_report))
    n_fine_inaccessible_ctx = int(
        np.count_nonzero(np.isin(regions, list(inaccessible)))
    )

    fs = (
        float(fine_spacing)
        if fine_spacing is not None
        else segmentation_spacing_from_results(ase_atoms, shape_f, results)
    )
    grid_fine_da, grid_shape_f = dask_grid_over_atoms(ase_atoms, spacing=fs)
    if tuple(grid_shape_f) != tuple(shape_f):
        raise ValueError(
            f"Fine grid shape {grid_shape_f} != regions shape {shape_f}; "
            f"adjust fine_spacing (using {fs})."
        )

    n_sph = len(spheres)
    rag_adj_labels_sorted: list[int] = []
    n_fine_interface = 0
    n_coarse_total = None
    sample_spacing_out = None
    region_at_sample = None
    voxel_vol = float(vol_cell / max(int(np.prod(shape_f)), 1))

    if resolution == "fine":
        grid_coordinate = grid_fine_da.compute().reshape(-1, 3)
        flat_d = denom_mask.ravel()
        pts = grid_coordinate[flat_d]
        n_grid_total = int(np.prod(shape_f))
        voxel_vol = vol_cell / max(n_grid_total, 1)
        if per_label and pts.shape[0] > 0:
            region_at_sample = regions.ravel()[flat_d]

    elif resolution == "coarse":
        ss = float(sample_spacing)
        sample_spacing_out = ss
        grid_coarse_da, shape_c = dask_grid_over_atoms(ase_atoms, spacing=ss)
        Nx, Ny, Nz = (int(x) for x in shape_c)
        n_coarse_total = int(Nx * Ny * Nz)
        nx, ny, nz = (int(x) for x in shape_f)
        Ic = np.arange(Nx, dtype=np.float64)[:, None, None]
        Jc = np.arange(Ny, dtype=np.float64)[None, :, None]
        Kc = np.arange(Nz, dtype=np.float64)[None, None, :]
        fx = Ic / max(Nx - 1, 1)
        fy = Jc / max(Ny - 1, 1)
        fz = Kc / max(Nz - 1, 1)
        ii = np.clip(np.rint(fx * max(nx - 1, 1)).astype(np.intp), 0, nx - 1)
        jj = np.clip(np.rint(fy * max(ny - 1, 1)).astype(np.intp), 0, ny - 1)
        kk = np.clip(np.rint(fz * max(nz - 1, 1)).astype(np.intp), 0, nz - 1)
        submask = denom_mask[ii, jj, kk]
        grid_coordinate = grid_coarse_da.compute().reshape(-1, 3)
        pts = grid_coordinate[submask.ravel()]
        voxel_vol = vol_cell / max(n_coarse_total, 1)
        if per_label and pts.shape[0] > 0:
            region_at_sample = regions[ii, jj, kk][submask]

    else:
        ss = float(sample_spacing)
        sample_spacing_out = ss
        rag_adj = labels_accessible_adjacent_to_inaccessible_via_rag(
            results, percolated, inaccessible
        )
        rag_adj_labels_sorted = sorted(int(x) for x in rag_adj)
        interface_fine = accessible_interface_voxel_mask(
            results,
            probe_radius,
            percolated_labels=percolated,
            inaccessible_labels=inaccessible,
            restrict_to_rag_adjacent_accessible=restrict_to_rag_adjacent_accessible,
            require_face_neighbor_inaccessible=require_face_neighbor_inaccessible,
            accessible_dgrid_ge_probe=accessible_dgrid_ge_probe,
        )
        n_fine_interface = int(np.count_nonzero(interface_fine))
        grid_coarse_da, shape_c = dask_grid_over_atoms(ase_atoms, spacing=ss)
        Nx, Ny, Nz = (int(x) for x in shape_c)
        n_coarse_total = int(Nx * Ny * Nz)
        nx, ny, nz = (int(x) for x in shape_f)
        Ic = np.arange(Nx, dtype=np.float64)[:, None, None]
        Jc = np.arange(Ny, dtype=np.float64)[None, :, None]
        Kc = np.arange(Nz, dtype=np.float64)[None, None, :]
        fx = Ic / max(Nx - 1, 1)
        fy = Jc / max(Ny - 1, 1)
        fz = Kc / max(Nz - 1, 1)
        ii = np.clip(np.rint(fx * max(nx - 1, 1)).astype(np.intp), 0, nx - 1)
        jj = np.clip(np.rint(fy * max(ny - 1, 1)).astype(np.intp), 0, ny - 1)
        kk = np.clip(np.rint(fz * max(nz - 1, 1)).astype(np.intp), 0, nz - 1)
        submask = interface_fine[ii, jj, kk]
        grid_coordinate = grid_coarse_da.compute().reshape(-1, 3)
        pts = grid_coordinate[submask.ravel()]
        voxel_vol = vol_cell / max(n_coarse_total, 1)
        if per_label and pts.shape[0] > 0:
            region_at_sample = regions[ii, jj, kk][submask]

    n_pts = int(pts.shape[0])
    per_label_dict: dict = {}

    def base_out():
        return {
            "denominator_domain": denominator_domain,
            "resolution": resolution,
            "probe_radius": float(probe_radius),
            "fine_spacing": float(fs),
            "sample_spacing": float(sample_spacing_out)
            if sample_spacing_out is not None
            else None,
            "accessible_dgrid_ge_probe": bool(accessible_dgrid_ge_probe),
            "n_spheres": int(n_sph),
            "n_denominator_points": n_pts,
            "percolated_labels": sorted(int(x) for x in percolated),
            "inaccessible_labels": sorted(int(x) for x in inaccessible),
            "n_regions_percolating": len(percolated),
            "n_regions_inaccessible": len(inaccessible),
            "n_fine_voxels_accessible": n_fine_accessible_ctx,
            "n_fine_voxels_inaccessible": n_fine_inaccessible_ctx,
            "n_coarse_grid_points": n_coarse_total,
            "restrict_to_rag_adjacent_accessible": bool(
                restrict_to_rag_adjacent_accessible
            ),
            "require_face_neighbor_inaccessible": bool(
                require_face_neighbor_inaccessible
            ),
            "n_rag_adjacent_accessible_labels": len(rag_adj_labels_sorted),
            "rag_adjacent_accessible_labels": rag_adj_labels_sorted,
            "n_fine_interface_voxels": n_fine_interface,
            "approx_voxel_volume_proxy_ang3": float(voxel_vol),
        }

    if n_pts == 0 or n_sph == 0:
        logger.debug("blocking_spheres_overlap_leakage: no points or no spheres")
        out = base_out()
        out.update(
            {
                "n_overlapped_mic": 0,
                "frac_overlapped_mic": float("nan"),
                "n_overlapped_raw": 0,
                "frac_overlapped_raw": float("nan"),
                "approx_overlap_volume_mic_ang3": 0.0,
                "frac_unit_cell_overlapped_mic": 0.0,
                "per_label": per_label_dict,
            }
        )
        _attach_legacy_overlap_keys(out)
        return out

    centers, radii = _blocking_spheres_cart_centers_radii(spheres, ase_atoms)
    covered_mic, covered_raw = _sphere_union_mic_raw(
        pts, centers, radii, ase_atoms, atol, chunk_size
    )
    n_ov_mic = int(np.count_nonzero(covered_mic))
    n_ov_raw = int(np.count_nonzero(covered_raw))
    vol_ov_mic = n_ov_mic * voxel_vol

    if per_label and region_at_sample is not None:
        rl = np.asarray(region_at_sample, dtype=np.intp)
        for L in np.unique(rl):
            m = rl == L
            n_L = int(np.count_nonzero(m))
            ov_m = int(np.count_nonzero(covered_mic[m]))
            ov_r = int(np.count_nonzero(covered_raw[m]))
            per_label_dict[int(L)] = {
                "n_coarse_samples": n_L,
                "n_overlapped_mic": ov_m,
                "frac_overlapped_mic": float(ov_m / n_L) if n_L else float("nan"),
                "n_overlapped_raw": ov_r,
                "frac_overlapped_raw": float(ov_r / n_L) if n_L else float("nan"),
            }

    out = base_out()
    out.update(
        {
            "n_overlapped_mic": n_ov_mic,
            "frac_overlapped_mic": float(n_ov_mic / n_pts),
            "n_overlapped_raw": n_ov_raw,
            "frac_overlapped_raw": float(n_ov_raw / n_pts),
            "approx_overlap_volume_mic_ang3": float(vol_ov_mic),
            "frac_unit_cell_overlapped_mic": float(vol_ov_mic / vol_cell)
            if vol_cell > 0
            else 0.0,
            "per_label": per_label_dict,
        }
    )
    _attach_legacy_overlap_keys(out)
    return out


def blocking_spheres_overlap_leakage_pair(
    results,
    ase_atoms,
    spheres,
    probe_radius,
    *,
    resolution="fine",
    sample_spacing=0.5,
    fine_spacing=None,
    accessible_dgrid_ge_probe=False,
    restrict_to_rag_adjacent_accessible=True,
    require_face_neighbor_inaccessible=True,
    per_label=False,
    atol=1e-4,
    chunk_size=65536,
    inaccessible_resolution=None,
):
    """
    Pair of :func:`blocking_spheres_overlap_leakage` runs for **BlockPockets** QA.

    When spheres are built to cover **inaccessible** void, you usually want both:

    1. **Leakage** — ``denominator_domain="accessible"``: among accessible sample
       points, what fraction fall **inside** the sphere union (spill into transport
       void)?
    2. **Coverage gap** — ``denominator_domain="inaccessible"``: among inaccessible
       sample points, what fraction fall **outside** the union (missed cages)?
       Reported as ``frac_covered_inaccessible_mic`` and ``frac_missed_inaccessible_mic``.

    ``resolution`` is passed to the **accessible** leg. The **inaccessible** leg uses
    the same resolution except when ``resolution="interface"`` (interface is only
    defined for accessible); then inaccessible defaults to ``"fine"`` unless you set
    ``inaccessible_resolution``.

    Parameters
    ----------
    inaccessible_resolution : str, optional
        ``"fine"`` or ``"coarse"`` for the inaccessible leg. Default: same as
        ``resolution``, or ``"fine"`` if ``resolution=="interface"``.
    """
    resolution = str(resolution).lower()
    if inaccessible_resolution is None:
        inc_res = "fine" if resolution == "interface" else resolution
    else:
        inc_res = str(inaccessible_resolution).lower()
    if inc_res not in {"fine", "coarse"}:
        raise ValueError(
            "inaccessible_resolution must be 'fine' or 'coarse', "
            f"got {inaccessible_resolution!r}"
        )

    common = dict(
        sample_spacing=sample_spacing,
        fine_spacing=fine_spacing,
        accessible_dgrid_ge_probe=accessible_dgrid_ge_probe,
        restrict_to_rag_adjacent_accessible=restrict_to_rag_adjacent_accessible,
        require_face_neighbor_inaccessible=require_face_neighbor_inaccessible,
        per_label=per_label,
        atol=atol,
        chunk_size=chunk_size,
    )
    st_acc = blocking_spheres_overlap_leakage(
        results,
        ase_atoms,
        spheres,
        probe_radius,
        denominator_domain="accessible",
        resolution=resolution,
        **common,
    )
    st_inc = blocking_spheres_overlap_leakage(
        results,
        ase_atoms,
        spheres,
        probe_radius,
        denominator_domain="inaccessible",
        resolution=inc_res,
        **common,
    )
    n_inc = int(st_inc.get("n_denominator_points", 0))
    cov_mic = st_inc.get("frac_overlapped_mic")
    n_cov_mic = int(st_inc.get("n_overlapped_mic", 0))
    missed_n = n_inc - n_cov_mic
    missed_frac = (
        float(missed_n / n_inc)
        if n_inc > 0 and cov_mic == cov_mic
        else float("nan")
    )
    cov_frac = float(cov_mic) if cov_mic == cov_mic else float("nan")

    return {
        "leakage_into_accessible": st_acc,
        "coverage_of_inaccessible": st_inc,
        "inaccessible_resolution_used": inc_res,
        "n_inaccessible_denominator": n_inc,
        "n_inaccessible_covered_mic": n_cov_mic,
        "n_inaccessible_missed_mic": missed_n,
        "frac_covered_inaccessible_mic": cov_frac,
        "frac_missed_inaccessible_mic": missed_frac,
    }


def print_blocking_spheres_overlap_pair_report(pair: dict) -> None:
    """
    Print :func:`blocking_spheres_overlap_leakage_pair` results in plain language.
    """
    st_a = pair["leakage_into_accessible"]
    st_i = pair["coverage_of_inaccessible"]
    inc_res = pair.get("inaccessible_resolution_used", st_i.get("resolution"))

    print("\n=== Blocking spheres — paired void QA (BlockPockets-style) ===\n")
    print(
        f"Probe: {st_a.get('probe_radius')} Å | spheres: {st_a.get('n_spheres')} | "
        f"accessible leg resolution: {st_a.get('resolution')} | "
        f"inaccessible leg resolution: {inc_res}\n"
    )

    na = int(st_a.get("n_denominator_points", 0))
    oa = int(st_a.get("n_overlapped_mic", 0))
    fa = st_a.get("frac_overlapped_mic")
    fa_pct = 100.0 * fa if fa == fa else float("nan")
    print(
        "(1) **Leakage into ACCESSIBLE void** — accessible sample points inside "
        "the sphere union (unintended):\n"
        f"    MIC: {oa:,} / {na:,}  ({fa_pct:.4f}%)\n"
        f"    Volume proxy (MIC): {st_a.get('approx_overlap_volume_mic_ang3', 0):.2f} Å³ "
        f"≈ {100.0 * st_a.get('frac_unit_cell_overlapped_mic', 0):.4f}% of cell\n"
    )

    ni = int(pair.get("n_inaccessible_denominator", st_i.get("n_denominator_points", 0)))
    n_cov = int(pair.get("n_inaccessible_covered_mic", st_i.get("n_overlapped_mic", 0)))
    n_miss = int(pair.get("n_inaccessible_missed_mic", ni - n_cov))
    f_cov = pair.get("frac_covered_inaccessible_mic")
    f_miss = pair.get("frac_missed_inaccessible_mic")
    if f_cov != f_cov:
        f_cov = float("nan")
        f_miss = float("nan")
    else:
        f_cov = float(f_cov)
        f_miss = float(f_miss)
    print(
        "(2) **Coverage of INACCESSIBLE void** — inaccessible sample points that "
        "should be blocked:\n"
        f"    Inside sphere union (MIC, **covered**): {n_cov:,} / {ni:,}  "
        f"({100.0 * f_cov:.4f}%)\n"
        f"    **Outside** union (MIC, **missed / gap**): {n_miss:,} / {ni:,}  "
        f"({100.0 * f_miss:.4f}%)\n"
        f"    Volume proxy for covered fraction (MIC): "
        f"{st_i.get('approx_overlap_volume_mic_ang3', 0):.2f} Å³\n"
    )


def blocking_spheres_overlap_into_accessible_regions(
    results,
    ase_atoms,
    spheres,
    probe_radius,
    spacing=None,
    accessible_dgrid_ge_probe=False,
    atol=1e-4,
    chunk_size=65536,
):
    """
    Estimate how much blocking-sphere volume intrudes into **accessible** void.

    Accessible voxels = segmentation labels that **percolate** at ``probe_radius``
    (same definition as :func:`percolated_region_labels_winding`). Optionally
    restrict to voxels with ``dgrid >= probe_radius`` (``accessible_dgrid_ge_probe``).

    For each **grid point** in that accessible set, test whether its Cartesian
    position lies inside the **union** of blocking spheres (fractional center +
    radius in Å). Uses **minimum-image** distance (PBC) for the primary metric,
    and also reports **raw** Cartesian distance (no PBC) for comparison.

    The grid is the same ordering as ``dask_grid_over_atoms`` / ``dgrid.ravel()``.

    Parameters
    ----------
    results : dict
        Must contain ``regions``, ``dgrid``.
    ase_atoms : ase.Atoms
    spheres : sequence of [sx, sy, sz, r]
        Fractional centers, radii in Å (RASPA convention).
    probe_radius : float
        Percolation threshold defining accessible labels.
    spacing : float, optional
        Grid spacing; inferred from ``results`` / cell if omitted.
    accessible_dgrid_ge_probe : bool
        If True, count only accessible-label voxels with ``dgrid >= probe_radius``.
    atol : float
        Tolerance on ``distance <= r + atol``.
    chunk_size : int
        Chunk size over accessible points for memory.

    Returns
    -------
    dict
        Keys include ``n_spheres``, ``n_accessible_voxels`` (denominator),
        ``n_overlapped_mic``, ``frac_overlapped_mic``, ``n_overlapped_raw``,
        ``frac_overlapped_raw``, ``approx_voxel_volume_ang3``,
        ``approx_overlap_volume_mic_ang3``, ``frac_unit_cell_overlapped_mic``.

    Implemented via :func:`blocking_spheres_overlap_leakage` (``resolution="fine"``).
    """
    return blocking_spheres_overlap_leakage(
        results,
        ase_atoms,
        spheres,
        probe_radius,
        denominator_domain="accessible",
        resolution="fine",
        fine_spacing=spacing,
        accessible_dgrid_ge_probe=accessible_dgrid_ge_probe,
        atol=atol,
        chunk_size=chunk_size,
        per_label=False,
    )


def print_blocking_spheres_overlap_report(stats: dict) -> None:
    """
    Print a human-readable summary for :func:`blocking_spheres_overlap_leakage`
    (and thin wrappers that return the same dict shape).
    """
    res = stats.get("resolution")
    if res is None:
        # Stats from older helpers (before ``resolution`` / ``denominator_domain``).
        if "n_sample_points_interface" in stats and "n_fine_interface_voxels" in stats:
            res = "interface"
        elif stats.get("n_coarse_grid_points"):
            res = "coarse"
        else:
            res = "fine"
    dom = stats.get("denominator_domain", "accessible")
    n = int(stats.get("n_denominator_points", 0))

    if res == "interface":
        n = int(stats.get("n_sample_points_interface", stats.get("n_denominator_points", 0)))
        if not n:
            print(
                "[blocking overlap] No coarse sample points on interface mask "
                f"(fine_interface_voxels={stats.get('n_fine_interface_voxels', 0)})."
            )
            return
        print(
            "[blocking overlap] Interface band (accessible near inaccessible)\n"
            f"  Coarse sample spacing={stats.get('sample_spacing', 0):.3f} Å, "
            f"fine grid={stats.get('fine_spacing', 0):.4f} Å\n"
            f"  RAG-adjacent accessible labels={stats.get('n_rag_adjacent_accessible_labels', 0)}: "
            f"{stats.get('rag_adjacent_accessible_labels')}\n"
            f"  Fine voxels in interface mask: {stats.get('n_fine_interface_voxels', 0)}; "
            f"coarse samples: {n} / {stats.get('n_coarse_grid_points', 0)} coarse cells"
        )
        print(
            f"  Inside sphere union (MIC): {stats['n_overlapped_mic']}/{n} "
            f"({100.0 * stats['frac_overlapped_mic']:.2f}%)"
        )
        print(
            f"  Approx. overlap volume (proxy): {stats['approx_overlap_volume_mic_ang3']:.2f} Å³ "
            f"≈ {100.0 * stats['frac_unit_cell_overlapped_mic']:.4f}% of unit cell"
        )
        print(
            f"  Raw (no PBC): {stats['n_overlapped_raw']}/{n} "
            f"({100.0 * stats['frac_overlapped_raw']:.2f}%)"
        )
        return

    if res == "fine":
        if dom == "accessible":
            nacc = int(stats.get("n_accessible_voxels", stats.get("n_denominator_points", 0)))
            if not nacc:
                print(
                    "[blocking overlap] No accessible voxels in denominator "
                    f"(probe={stats.get('probe_radius')} Å, dgrid≥probe="
                    f"{stats.get('accessible_dgrid_ge_probe')})."
                )
                return
            print(
                "[blocking overlap] Fine grid — accessible void (MIC)\n"
                f"  {stats['n_overlapped_mic']}/{nacc} "
                f"({100.0 * stats['frac_overlapped_mic']:.2f}%)"
            )
            print(
                f"  Approx. overlap volume: {stats['approx_overlap_volume_mic_ang3']:.2f} Å³ "
                f"≈ {100.0 * stats['frac_unit_cell_overlapped_mic']:.4f}% of unit cell"
            )
            print(
                f"  Raw (no PBC): {stats['n_overlapped_raw']}/{nacc} "
                f"({100.0 * stats['frac_overlapped_raw']:.2f}%)"
            )
        else:
            ni = int(stats.get("n_inaccessible_sample_points", stats.get("n_denominator_points", 0)))
            if not ni:
                print("[blocking overlap] No inaccessible voxels in denominator.")
                return
            print(
                "[blocking overlap] Fine grid — **inaccessible** void (MIC)\n"
                f"  {stats['n_overlapped_mic']}/{ni} "
                f"({100.0 * stats['frac_overlapped_mic']:.2f}%)"
            )
            print(
                f"  Approx. overlap volume: {stats['approx_overlap_volume_mic_ang3']:.2f} Å³ "
                f"≈ {100.0 * stats['frac_unit_cell_overlapped_mic']:.4f}% of unit cell"
            )
            print(
                f"  Raw (no PBC): {stats['n_overlapped_raw']}/{ni} "
                f"({100.0 * stats['frac_overlapped_raw']:.2f}%)"
            )
        return

    # coarse
    if dom == "accessible":
        n = int(stats.get("n_accessible_sample_points", stats.get("n_denominator_points", 0)))
        if not n:
            print("No coarse sample points on accessible mask.")
            return
        ncg = int(stats.get("n_coarse_grid_points", 0))
        pct_cells = 100.0 * n / ncg if ncg else 0.0
        print(
            "\n=== Blocking spheres vs accessible void (coarse sample) ===\n"
            f"Probe / percolation threshold: {stats['probe_radius']:.3f} Å\n"
            f"Segmentation (fine) grid spacing: {stats['fine_spacing']:.4f} Å\n"
            f"Sample grid spacing: {stats.get('sample_spacing', 0):.3f} Å  "
            "(each sample uses **nearest** fine voxel for region / accessibility)\n"
            f"dgrid ≥ probe filter on accessible: {stats.get('accessible_dgrid_ge_probe')}\n"
        )
        print(
            "**Inaccessible** region labels (non-percolating at this probe — typical "
            "``BlockPockets`` targets):\n  "
            f"{stats.get('inaccessible_labels')}\n"
            f"  ({stats.get('n_regions_inaccessible')} regions)\n"
        )
        print(
            "**Accessible** region labels (percolating — main void / transport):\n  "
            f"{stats.get('percolated_labels')}\n"
            f"  ({stats.get('n_regions_percolating')} regions)\n"
        )
        print(
            "Fine-grid voxel counts:\n"
            f"  Accessible-mask voxels: {stats.get('n_fine_voxels_accessible')}\n"
            f"  Inaccessible-mask voxels: {stats.get('n_fine_voxels_inaccessible')}\n"
        )
        print(
            f"Coarse grid: {ncg} cells; {n} centers classified as accessible "
            f"({pct_cells:.2f}% of coarse cells).\n"
            "Of those **accessible** sample points, how many lie inside the **union** of "
            "blocking spheres:\n"
            f"  MIC (PBC): {stats['n_overlapped_mic']}/{n}  "
            f"({100.0 * stats['frac_overlapped_mic']:.2f}%)\n"
            f"  Raw (no PBC): {stats['n_overlapped_raw']}/{n}  "
            f"({100.0 * stats['frac_overlapped_raw']:.2f}%)\n"
            f"Approx. overlap volume (coarse voxel × count, MIC): "
            f"{stats['approx_overlap_volume_mic_ang3']:.2f} Å³  "
            f"(~{100.0 * stats['frac_unit_cell_overlapped_mic']:.4f}% of unit cell)\n"
        )
        pl = stats.get("per_accessible_label") or stats.get("per_label") or {}
        if pl:
            print(
                "Breakdown by **accessible** region label "
                "(coarse samples / overlapped / % MIC):\n"
            )
            for lab in sorted(pl.keys()):
                row = pl[lab]
                print(
                    f"  Label {lab:2d}:  samples={row['n_coarse_samples']:5d}  "
                    f"inside spheres={row['n_overlapped_mic']:5d}  "
                    f"({100.0 * row['frac_overlapped_mic']:.2f}%)"
                )
    else:
        n = int(stats.get("n_inaccessible_sample_points", stats.get("n_denominator_points", 0)))
        if not n:
            print("No coarse sample points on inaccessible mask.")
            return
        print(
            "\n=== Blocking overlap (coarse sample, **inaccessible** denominator) ===\n"
            f"Probe: {stats['probe_radius']:.3f} Å | sample spacing: {stats.get('sample_spacing', 0):.3f} Å\n"
            f"Coarse samples in inaccessible void: {n} / {stats.get('n_coarse_grid_points', 0)}\n"
            f"MIC: {stats['n_overlapped_mic']}/{n} ({100.0 * stats['frac_overlapped_mic']:.2f}%)\n"
            f"Raw: {stats['n_overlapped_raw']}/{n} ({100.0 * stats['frac_overlapped_raw']:.2f}%)\n"
        )


def print_blocking_overlap_stats(stats: dict) -> None:
    """Deprecated name; calls :func:`print_blocking_spheres_overlap_report`."""
    print_blocking_spheres_overlap_report(stats)


def blocking_spheres_overlap_accessible_coarse_sample(
    results,
    ase_atoms,
    spheres,
    probe_radius,
    *,
    sample_spacing=0.5,
    fine_spacing=None,
    accessible_dgrid_ge_probe=False,
    atol=1e-4,
    chunk_size=65536,
):
    """
    Same intent as :func:`blocking_spheres_overlap_into_accessible_regions`, but
    **sample** accessible void on a coarser grid (e.g. 0.5 Å) for speed.

    Each coarse cell center maps to the **nearest** fine segmentation voxel
    (same fractional mapping as :func:`mofography.compute.dask_grid_over_atoms`).
    A sample point counts as **accessible** if that fine voxel lies in a label
    that **percolates** at ``probe_radius`` (and optionally ``dgrid >= probe``).

    Returns label sets, fine-grid counts for context, and a per-label
    breakdown of how many coarse samples per percolating region fell inside the
    blocking-sphere union (MIC and raw).

    Implemented via :func:`blocking_spheres_overlap_leakage` (``resolution="coarse"``,
    ``per_label=True``).

    Parameters
    ----------
    sample_spacing : float
        Coarse probe grid in Å (e.g. 0.5).
    fine_spacing : float, optional
        Spacing of ``results['regions']``; inferred if omitted.
    """
    return blocking_spheres_overlap_leakage(
        results,
        ase_atoms,
        spheres,
        probe_radius,
        denominator_domain="accessible",
        resolution="coarse",
        sample_spacing=sample_spacing,
        fine_spacing=fine_spacing,
        accessible_dgrid_ge_probe=accessible_dgrid_ge_probe,
        atol=atol,
        chunk_size=chunk_size,
        per_label=True,
    )


def print_blocking_accessible_coarse_overlap_report(stats: dict) -> None:
    """Deprecated name; calls :func:`print_blocking_spheres_overlap_report`."""
    print_blocking_spheres_overlap_report(stats)


def print_blocking_accessible_coarse_overlap_comparison(
    stats_by_method: dict,
    *,
    title: str = "Coarse accessible overlap — method comparison",
) -> None:
    """
    Compare outputs of :func:`blocking_spheres_overlap_accessible_coarse_sample`
    for several blocking builders (e.g. largest-including vs greedy).

    Parameters
    ----------
    stats_by_method : dict
        Maps method name (str) -> stats dict from
        :func:`blocking_spheres_overlap_accessible_coarse_sample`.
    """
    print(f"\n{'=' * 60}\n{title}\n{'=' * 60}")
    print(
        "Lower **MIC %** and **volume proxy** = less blocking into percolating void "
        "(better for unintended-domain encroachment).\n"
    )
    hdr = (
        f"{'Method':<22} {'N spheres':>10}  {'MIC %':>8}  {'raw %':>8}  "
        f"{'overlap Å³':>12}  {'% unit cell':>12}"
    )
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for name, st in stats_by_method.items():
        n = int(
            st.get(
                "n_accessible_sample_points", st.get("n_denominator_points", 0)
            )
        )
        if not n:
            micp = float("nan")
            rawp = float("nan")
        else:
            micp = 100.0 * st["frac_overlapped_mic"]
            rawp = 100.0 * st["frac_overlapped_raw"]
        rows.append(
            (
                name,
                int(st.get("n_spheres", 0)),
                micp,
                rawp,
                st.get("approx_overlap_volume_mic_ang3", 0.0),
                100.0 * st.get("frac_unit_cell_overlapped_mic", 0.0),
            )
        )
    for name, ns, micp, rawp, vol, fc in rows:
        print(
            f"{name:<22} {ns:10d}  {micp:7.2f}%  {rawp:7.2f}%  "
            f"{vol:12.2f}  {fc:11.4f}%"
        )
    first = next(iter(stats_by_method.values()))
    n0 = first.get("n_accessible_sample_points")
    if n0 is None:
        n0 = first.get("n_denominator_points")
    if n0 is not None:
        print(
            f"\n(Denominator: {n0} coarse accessible sample points — same for all methods.)"
        )


def labels_accessible_adjacent_to_inaccessible_via_rag(
    results,
    percolated_labels,
    inaccessible_labels,
):
    """
    Region labels that are **percolating** and share an **internal** RAG edge
    with at least one **inaccessible** label.

    Uses :func:`mofography.compute.rag_from_connections` (same window graph as
    visualization / pocket connectivity), not 6-neighbor grid adjacency.
    """
    import numpy as np

    from .compute import rag_from_connections

    connections = results.get("connections")
    all_maxima = results.get("all_maxima")
    if connections is None or all_maxima is None:
        raise ValueError(
            "results must include 'connections' and 'all_maxima' (periodic segmentation output) "
            "to build the region adjacency graph."
        )
    maxima_arr = np.asarray(all_maxima)
    rag = rag_from_connections(connections, maxima_arr, only_use_internal=True)
    percolated = set(int(x) for x in percolated_labels)
    inaccessible = set(int(x) for x in inaccessible_labels)
    out = set()
    for inc in inaccessible:
        if inc not in rag:
            continue
        for nb in rag.neighbors(inc):
            if nb in percolated:
                out.add(nb)
    return out


def accessible_interface_voxel_mask(
    results,
    probe_radius,
    *,
    percolated_labels=None,
    inaccessible_labels=None,
    restrict_to_rag_adjacent_accessible=True,
    require_face_neighbor_inaccessible=True,
    accessible_dgrid_ge_probe=False,
):
    """
    Boolean mask (same shape as ``results['regions']``) for **accessible**
    voxels near inaccessible void, for overlap / leakage estimates.

    *Optional RAG filter:* keep only voxels whose region label lies in
    :func:`labels_accessible_adjacent_to_inaccessible_via_rag`.

    *Optional grid filter:* keep only voxels with a 6-neighbor (periodic)
    belonging to an **inaccessible** region (face contact on the segmentation
    grid).

    Percolation at ``probe_radius`` matches :func:`percolated_region_labels_winding`.
    """
    import numpy as np

    regions = np.asarray(results["regions"])
    dgrid = np.asarray(results["dgrid"], dtype=float)
    if percolated_labels is None:
        percolated_labels = percolated_region_labels_winding(results, probe_radius)
    if inaccessible_labels is None:
        all_labels = {int(x) for x in np.unique(regions) if int(x) != 0}
        inaccessible_labels = all_labels - set(percolated_labels)

    percolated = set(int(x) for x in percolated_labels)
    inaccessible = set(int(x) for x in inaccessible_labels)

    if restrict_to_rag_adjacent_accessible:
        allowed_regions = labels_accessible_adjacent_to_inaccessible_via_rag(
            results, percolated, inaccessible
        )
    else:
        allowed_regions = percolated

    mask_acc = np.isin(regions, list(allowed_regions))
    if accessible_dgrid_ge_probe:
        mask_acc &= dgrid >= float(probe_radius)

    if require_face_neighbor_inaccessible:
        inc_mask = np.isin(regions, list(inaccessible))
        touch = np.zeros_like(regions, dtype=bool)
        touch |= np.roll(inc_mask, -1, axis=0)
        touch |= np.roll(inc_mask, 1, axis=0)
        touch |= np.roll(inc_mask, -1, axis=1)
        touch |= np.roll(inc_mask, 1, axis=1)
        touch |= np.roll(inc_mask, -1, axis=2)
        touch |= np.roll(inc_mask, 1, axis=2)
        mask_acc &= touch

    return mask_acc


def blocking_spheres_overlap_into_accessible_near_interface(
    results,
    ase_atoms,
    spheres,
    probe_radius,
    *,
    sample_spacing=0.5,
    fine_spacing=None,
    accessible_dgrid_ge_probe=False,
    restrict_to_rag_adjacent_accessible=True,
    require_face_neighbor_inaccessible=True,
    atol=1e-4,
    chunk_size=65536,
):
    """
    Cheaper overlap diagnostic: sample **coarse** grid points (default 0.5 Å)
    that lie on **interface** accessible voxels (see :func:`accessible_interface_voxel_mask`).

    Each coarse voxel center is classified by **nearest** fine-grid voxel
    (same fractional mapping as :func:`mofography.compute.dask_grid_over_atoms`).

    The reported fraction is **only** over this interface sample — not over all
    accessible void. A **high** percentage here often means the blocking union
    covers the narrow gasket next to pockets (expected for large including
    spheres), not that the bulk percolating channel is blocked; compare with
    :func:`blocking_spheres_overlap_into_accessible_regions` for a global
    accessible denominator (much more expensive).

    Implemented via :func:`blocking_spheres_overlap_leakage` (``resolution="interface"``).

    Parameters
    ----------
    sample_spacing : float
        Spacing (Å) for probe points; coarser => fewer distance evaluations.
    fine_spacing : float, optional
        Spacing of ``results['regions']``; inferred if omitted.
    restrict_to_rag_adjacent_accessible, require_face_neighbor_inaccessible
        Passed to :func:`accessible_interface_voxel_mask`.
    """
    return blocking_spheres_overlap_leakage(
        results,
        ase_atoms,
        spheres,
        probe_radius,
        denominator_domain="accessible",
        resolution="interface",
        sample_spacing=sample_spacing,
        fine_spacing=fine_spacing,
        accessible_dgrid_ge_probe=accessible_dgrid_ge_probe,
        restrict_to_rag_adjacent_accessible=restrict_to_rag_adjacent_accessible,
        require_face_neighbor_inaccessible=require_face_neighbor_inaccessible,
        atol=atol,
        chunk_size=chunk_size,
        per_label=False,
    )


def print_blocking_interface_overlap_stats(stats: dict) -> None:
    """Deprecated name; calls :func:`print_blocking_spheres_overlap_report`."""
    print_blocking_spheres_overlap_report(stats)


def _grid_frac_cart_helpers(shape, cell):
    """
    Fractional grid mapping and PBC distance in Å (ASE: rows = a,b,c;
    Cartesian = cell.T @ frac).
    """
    import numpy as np

    cell = np.asarray(cell, dtype=float)
    if cell.shape != (3, 3):
        cell = np.array(cell).reshape(3, 3)
    cell_T = cell.T

    def ijks_to_fractional(i, j, k):
        sx = i / max(1, shape[0] - 1) if shape[0] > 1 else 0.5
        sy = j / max(1, shape[1] - 1) if shape[1] > 1 else 0.5
        sz = k / max(1, shape[2] - 1) if shape[2] > 1 else 0.5
        return np.array([sx, sy, sz], dtype=float)

    def pbc_distance_angstrom(frac1, frac2):
        dr_frac = np.asarray(frac1, dtype=float) - np.asarray(frac2, dtype=float)
        dr_frac = dr_frac - np.round(dr_frac)
        dr_cart = cell_T @ dr_frac
        return float(np.linalg.norm(dr_cart))

    return ijks_to_fractional, pbc_distance_angstrom


def blocking_spheres_largest_including(
    results,
    ase_atoms,
    region_labels_to_cover,
    radius_scale=1.0,
):
    """
    Build a list of RASPA blocking spheres ``[sx, sy, sz, radius_Å]`` (fractional
    center + radius) that covers every voxel in the given segmentation regions.

    For each voxel in the selected regions, the largest including sphere
    (Gelb–Gubbins style) is found: center is the voxel C with maximum
    ``dgrid[C]`` such that PBC-distance(voxel, C) <= dgrid[C]. Many voxels map to
    the same (center, radius); unique (center, radius) pairs are returned.

    Use with ``region_labels_to_cover`` = inaccessible region labels to block
    non-percolating pockets for GCMC.

    Parameters
    ----------
    results : dict
        Must include ``regions`` and ``dgrid``.
    ase_atoms : ase.Atoms
        Unit cell for PBC distances.
    region_labels_to_cover : set of int
        Region labels whose voxels must be covered by the output spheres.
    radius_scale : float
        Multiply each sphere radius by this factor (>= 0).

    Returns
    -------
    list of list
        ``[sx, sy, sz, radius_Å]`` for each unique sphere.
    """
    import numpy as np

    regions = np.asarray(results["regions"])
    dgrid = np.asarray(results["dgrid"], dtype=float)
    shape = regions.shape
    cell = ase_atoms.get_cell()
    ijks_to_fractional, pbc_distance_angstrom = _grid_frac_cart_helpers(shape, cell)

    region_voxels = {}
    region_sorted = {}
    for L in region_labels_to_cover:
        L = int(L)
        if L < 1:
            continue
        mask = regions == L
        voxels = np.argwhere(mask)
        if len(voxels) == 0:
            continue
        region_voxels[L] = voxels
        vals = dgrid[mask]
        order = np.argsort(-vals)
        sorted_ijk = voxels[order]
        sorted_r = vals[order]
        region_sorted[L] = (sorted_ijk, sorted_r)

    seen = set()
    out_spheres = []
    for L, (sorted_ijk, sorted_r) in region_sorted.items():
        voxels = region_voxels[L]
        for v in voxels:
            vi, vj, vk = int(v[0]), int(v[1]), int(v[2])
            frac_v = ijks_to_fractional(vi, vj, vk)
            center_ijk = None
            radius_val = None
            for idx in range(len(sorted_ijk)):
                ci, cj, ck = (
                    int(sorted_ijk[idx][0]),
                    int(sorted_ijk[idx][1]),
                    int(sorted_ijk[idx][2]),
                )
                r_c = float(sorted_r[idx])
                if r_c <= 0:
                    continue
                frac_c = ijks_to_fractional(ci, cj, ck)
                dist = pbc_distance_angstrom(frac_v, frac_c)
                if dist <= r_c:
                    center_ijk = (ci, cj, ck)
                    radius_val = r_c
                    break
            if center_ijk is None:
                continue
            key = (center_ijk[0], center_ijk[1], center_ijk[2], radius_val)
            if key in seen:
                continue
            seen.add(key)
            sx, sy, sz = ijks_to_fractional(
                center_ijk[0], center_ijk[1], center_ijk[2]
            )
            r_ang = max(radius_val * radius_scale, 0.5)
            out_spheres.append([float(sx), float(sy), float(sz), r_ang])

    return out_spheres


def percolation_region_labels_union_winding(axis_clusters: dict) -> set[int]:
    """
    Union of segmentation region IDs appearing in any winding-percolation cluster.

    Parameters
    ----------
    axis_clusters : dict
        ``percolation_all_axes_winding`` output ``axis_clusters``: keys
        ``'x'``, ``'y'``, ``'z'`` mapping to lists of label groups.
    """
    labels: set[int] = set()
    for axis in ("x", "y", "z"):
        for group in axis_clusters[axis]:
            labels.update(int(x) for x in group)
    return labels


def blocking_spheres_greedy(
    results,
    ase_atoms,
    threshold,
    spacing=None,
    block_accessible=False,
    radius_scale=1.0,
):
    """
    Greedy blocking spheres (org-note style): place spheres at pocket maxima for
    each target region, then repeatedly add the sphere with the largest
    remaining ``dgrid`` value until all masked voxels are covered.

    Coverage uses **Euclidean** distances in Å (**no** minimum-image / PBC).
    Use :func:`blocking_spheres_largest_including` if you need PBC-consistent
    including spheres on the same mask.

    Parameters
    ----------
    results : dict
        Segmentation output: ``regions``, ``dgrid``, ``all_maxima``,
        ``maxima_info``, plus fields required by winding percolation
        (:func:`mofography.winding_percolation.percolation_all_axes_winding`).
    ase_atoms : ase.Atoms
    threshold : float
        Probe radius (Å) for winding percolation / ``dgrid`` mask when
        ``block_accessible`` is True.
    spacing : float, optional
        Grid spacing (Å); must match ``results['regions']``. Inferred from
        ``results`` / cell if omitted.
    block_accessible : bool, default False
        * **False** — Mask all voxels in **non-percolating** regions at
          ``threshold``. For standard RASPA ``BlockPockets`` (block side pockets).
        * **True** — Mask voxels in **percolating** regions with
          ``dgrid >= threshold``. For ``InvertBlockPockets``.
    radius_scale : float, default 1.0
        Multiply every sphere radius by this factor; each radius is at least
        0.5 Å (same floor as :func:`blocking_spheres_largest_including`).

    Returns
    -------
    spheres : list of [sx, sy, sz, r_Å]
        Fractional centers, radii in Å (RASPA convention).
    meta : dict
        ``block_accessible``, ``percolation_labels``, ``inaccessible_labels``,
        ``target_region_labels``, ``n_initial``, ``n_greedy``, ``grid_shape``,
        ``n_masked_voxels``.
    """
    import numpy as np

    from .compute import dask_grid_over_atoms
    from .helper import find_coord_from_indices, get_fractional_coordinates
    from .winding_percolation import (
        build_base_region_labels_from_results,
        make_directed_graph_from_segmentation_results,
        percolation_all_axes_winding,
    )

    dgrid = np.asarray(results["dgrid"], dtype=float)
    regions = np.asarray(results["regions"])
    if dgrid.shape != regions.shape:
        raise ValueError("dgrid and regions shapes differ")

    if spacing is None:
        spacing = _infer_spacing_from_results(ase_atoms, regions.shape, results=results)

    graph = make_directed_graph_from_segmentation_results(results)
    base_labels = build_base_region_labels_from_results(results)
    axis_clusters, _, _ = percolation_all_axes_winding(
        results,
        threshold=threshold,
        graph=graph,
        base_region_labels=base_labels,
        include_dangling=True,
    )
    percolation_labels = percolation_region_labels_union_winding(axis_clusters)
    all_labels = {int(x) for x in np.unique(regions) if int(x) != 0}
    inaccessible_labels = all_labels - percolation_labels

    dgrid_1d = dgrid.ravel()
    regions_1d = regions.ravel()
    mask = np.zeros(len(regions_1d), dtype=bool)

    if block_accessible:
        target_labels = sorted(percolation_labels)
        for i, n in enumerate(regions_1d):
            nn = int(n)
            if nn in percolation_labels and dgrid_1d[i] >= threshold:
                mask[i] = True
    else:
        target_labels = sorted(inaccessible_labels)
        for i, n in enumerate(regions_1d):
            if int(n) in inaccessible_labels:
                mask[i] = True

    if not target_labels or not mask.any():
        meta_empty = {
            "block_accessible": block_accessible,
            "percolation_labels": percolation_labels,
            "inaccessible_labels": inaccessible_labels,
            "target_region_labels": target_labels,
            "n_initial": 0,
            "n_greedy": 0,
            "grid_shape": tuple(int(x) for x in dgrid.shape),
            "n_masked_voxels": int(mask.sum()),
        }
        return [], meta_empty

    grid_da, grid_shape = dask_grid_over_atoms(ase_atoms, spacing=spacing)
    if tuple(grid_shape) != tuple(dgrid.shape):
        raise ValueError(
            f"dask_grid shape {grid_shape} != dgrid.shape {dgrid.shape} "
            f"(spacing={spacing} must match segmentation)"
        )
    grid_coordinate = grid_da.compute()

    all_maxima = np.asarray(results["all_maxima"])
    maxima_info = np.asarray(results["maxima_info"])
    n_reg = int(regions.max())
    maxima = all_maxima[:n_reg]
    centroid_coordinate = find_coord_from_indices(maxima, regions.shape, ase_atoms)
    all_maxima_radii = dgrid[tuple(maxima_info[:, 1:4].T)]
    maxima_radii = all_maxima_radii[:n_reg]

    target_arr = np.array(target_labels, dtype=int)
    initial_centers = centroid_coordinate[target_arr - 1]
    initial_scaled = get_fractional_coordinates(initial_centers, ase_atoms)
    initial_radii = maxima_radii[target_arr - 1].reshape(-1)

    grid_masked = grid_coordinate[mask]
    grid_d_masked = dgrid_1d[mask]

    m = np.ones(len(grid_masked), dtype=bool)
    for c, r in zip(initial_centers, initial_radii):
        dist2 = np.sum((grid_masked - c) ** 2, axis=1)
        m &= dist2 > r**2
    outside_coord = grid_masked[m]
    outside_dist = grid_d_masked[m]

    greedy_scaled = []
    greedy_r = []
    while len(outside_coord) > 0:
        idx = int(np.argmax(outside_dist))
        center = outside_coord[idx]
        r = float(outside_dist[idx])
        greedy_scaled.append(
            get_fractional_coordinates(center.reshape(1, 3), ase_atoms)[0]
        )
        greedy_r.append(r)
        dist2 = np.sum((outside_coord - center) ** 2, axis=1)
        keep = dist2 > r**2
        outside_coord = outside_coord[keep]
        outside_dist = outside_dist[keep]

    rs = float(radius_scale)
    spheres = []
    for i, s in enumerate(initial_scaled):
        r_ang = max(float(initial_radii[i]) * rs, 0.5)
        spheres.append([float(s[0]), float(s[1]), float(s[2]), r_ang])
    for s, r in zip(greedy_scaled, greedy_r):
        r_ang = max(float(r) * rs, 0.5)
        spheres.append([float(s[0]), float(s[1]), float(s[2]), r_ang])

    meta = {
        "block_accessible": block_accessible,
        "percolation_labels": percolation_labels,
        "inaccessible_labels": inaccessible_labels,
        "target_region_labels": target_labels,
        "n_initial": len(initial_scaled),
        "n_greedy": len(greedy_r),
        "grid_shape": tuple(int(x) for x in dgrid.shape),
        "n_masked_voxels": int(mask.sum()),
    }
    return spheres, meta


# Disabled 2026-03: one sphere per pocket at watershed maximum (dgrid at maxima) can
# **fail to cover** irregular / fragmented inaccessible regions. Prefer
# :func:`blocking_spheres_largest_including` or :func:`blocking_spheres_greedy`.
#
# def blocking_pockets_from_percolation(
#     results,
#     ase_atoms,
#     percolated_region_labels,
#     radius_scale=1.5,
# ):
#     """
#     One blocking sphere per *inaccessible* region: center at pocket maximum,
#     radius = pore radius at maximum × ``radius_scale``.
#
#     Faster than :func:`blocking_spheres_largest_including` but may not cover
#     irregular pockets fully.
#     """
#     import numpy as np
#
#     regions = np.asarray(results["regions"])
#     all_maxima = results["all_maxima"]
#     dgrid = np.asarray(results["dgrid"])
#     shape = regions.shape
#     n_regions = int(regions.max())
#     maxima_radii = np.zeros(n_regions)
#     for L in range(1, n_regions + 1):
#         i, j, k = all_maxima[L - 1]
#         maxima_radii[L - 1] = float(dgrid[i, j, k])
#     ijks_to_fractional, _ = _grid_frac_cart_helpers(shape, ase_atoms.get_cell())
#
#     percolated = set(int(x) for x in percolated_region_labels)
#     out = []
#     for region_label in range(1, n_regions + 1):
#         if region_label in percolated:
#             continue
#         i, j, k = all_maxima[region_label - 1]
#         sx, sy, sz = ijks_to_fractional(i, j, k)
#         r_ang = maxima_radii[region_label - 1] * radius_scale
#         r_ang = max(r_ang, 0.5)
#         out.append([float(sx), float(sy), float(sz), r_ang])
#     return out


def blocking_spheres_for_probe_radius(
    results,
    ase_atoms,
    probe_radius,
    radius_scale=1.0,
    method="largest_including",
    greedy_spacing=None,
    **percolation_kwargs,
):
    """
    Convenience: inaccessible regions at ``probe_radius``, then blocking spheres.

    Parameters
    ----------
    results, ase_atoms
        Segmentation results and framework.
    probe_radius : float
        Percolation threshold (Å).
    radius_scale : float
        Passed to the blocking builder.
    method : str
        * ``"largest_including"`` — :func:`blocking_spheres_largest_including`.
        * ``"greedy"`` — :func:`blocking_spheres_greedy` with
          ``block_accessible=False`` (standard ``BlockPockets``).
    greedy_spacing : float, optional
        Grid spacing for ``method=\"greedy\"``; inferred from ``results`` if omitted.
    **percolation_kwargs
        Forwarded to :func:`percolated_region_labels_winding`.

    Returns
    -------
    spheres : list of [sx, sy, sz, r]
    percolated : set of int
    inaccessible : set of int
    """
    import numpy as np

    regions = np.asarray(results["regions"])
    all_labels = {int(x) for x in np.unique(regions) if int(x) != 0}
    percolated = percolated_region_labels_winding(
        results, probe_radius, **percolation_kwargs
    )
    inaccessible = all_labels - percolated
    if method == "largest_including":
        spheres = blocking_spheres_largest_including(
            results, ase_atoms, inaccessible, radius_scale=radius_scale
        )
    elif method == "greedy":
        spheres, _meta = blocking_spheres_greedy(
            results,
            ase_atoms,
            probe_radius,
            spacing=greedy_spacing,
            block_accessible=False,
            radius_scale=radius_scale,
        )
    # elif method == "pocket_maximum":
    #     # Was :func:`blocking_pockets_from_percolation` (commented out; see above).
    #     spheres = blocking_pockets_from_percolation(
    #         results, ase_atoms, percolated, radius_scale=max(radius_scale, 1.0)
    #     )
    else:
        raise ValueError(
            f"method must be 'largest_including' or 'greedy', got {method!r}"
        )
    return spheres, percolated, inaccessible


# gRASPA LAMMPS snapshot I/O (Movies/System_0/result_*.data)
GRASPA_FRAMEWORK_MOL_ID = 1
GRASPA_XYLENE_TYPE_IDS = frozenset({12, 13, 14})


def read_graspa_lammps_snapshot(path_to_file):
    """Parse a gRASPA ``result_*.data`` snapshot (see :func:`mofography.data.read_graspa_lammps_data`)."""
    from .data import read_graspa_lammps_data

    return read_graspa_lammps_data(path_to_file)


def graspa_lammps_data_to_ase(
    data,
    *,
    framework_mol_id=GRASPA_FRAMEWORK_MOL_ID,
    guest_type_ids=GRASPA_XYLENE_TYPE_IDS,
):
    """Build ASE ``Atoms`` objects from :func:`read_graspa_lammps_data` output.

    Returns
    -------
  framework : ase.Atoms
        Framework atoms (default ``mol_id == 1``).
    guests : ase.Atoms
        All adsorbate atoms in one ``Atoms`` object.
    molecules : dict[int, ase.Atoms]
        One ``Atoms`` per guest ``mol_id`` (keys sorted).
    """
    from ase import Atoms
    import numpy as np

    cell = np.asarray(data["cell"], dtype=float)
    fw_idx, guest_idx, mol_groups = [], [], {}
    for i, rec in enumerate(data["atoms"]):
        is_guest = rec["mol_id"] != framework_mol_id or rec["type_id"] in guest_type_ids
        if is_guest:
            guest_idx.append(i)
            mol_groups.setdefault(rec["mol_id"], []).append(i)
        else:
            fw_idx.append(i)

    def _build(indices):
        if not indices:
            return Atoms(cell=cell, pbc=True)
        recs = [data["atoms"][i] for i in indices]
        atoms = Atoms(
            symbols=[r["symbol"] for r in recs],
            positions=np.array([r["position"] for r in recs], dtype=float),
            cell=cell,
            pbc=True,
        )
        atoms.arrays["charge"] = np.array([r["charge"] for r in recs], dtype=float)
        atoms.arrays["mol_id"] = np.array([r["mol_id"] for r in recs], dtype=int)
        atoms.arrays["type_id"] = np.array([r["type_id"] for r in recs], dtype=int)
        atoms.arrays["atom_id"] = np.array([r["id"] for r in recs], dtype=int)
        atoms.arrays["mass_label"] = np.array(
            [r["mass_label"] for r in recs], dtype=object
        )
        if any(r.get("molecule_name") for r in recs):
            atoms.arrays["molecule_name"] = np.array(
                [r.get("molecule_name") or "" for r in recs], dtype=object
            )
        return atoms

    framework = _build(fw_idx)
    guests = _build(guest_idx)
    molecules = {mid: _build(idxs) for mid, idxs in sorted(mol_groups.items())}
    for mid, at in molecules.items():
        names = at.arrays.get("molecule_name")
        if names is not None and len(names) and names[0]:
            at.info["molecule_name"] = str(names[0])
        at.info["mol_id"] = int(mid)
    return framework, guests, molecules


def get_atoms_from_graspa_lammps_snapshot(
    path_to_file,
    *,
    framework_mol_id=GRASPA_FRAMEWORK_MOL_ID,
    guest_type_ids=GRASPA_XYLENE_TYPE_IDS,
    split_guests=True,
):
    """Load one gRASPA ``.data`` snapshot into ASE objects.

    Parameters
    ----------
    path_to_file : str or Path
        ``result_*.data`` from gRASPA ``Movies/System_0/``.
    split_guests : bool
        If True (default), return ``(framework, guests, molecules)``.
        If False, return a single combined ``Atoms`` with ``arrays['is_guest']``.

    Returns
    -------
    tuple or ase.Atoms
    """
    from .data import read_graspa_lammps_data

    data = read_graspa_lammps_data(path_to_file)
    framework, guests, molecules = graspa_lammps_data_to_ase(
        data,
        framework_mol_id=framework_mol_id,
        guest_type_ids=guest_type_ids,
    )
    if split_guests:
        return framework, guests, molecules

    from ase import Atoms
    import numpy as np

    combined = framework + guests
    is_guest = np.zeros(len(combined), dtype=bool)
    is_guest[len(framework) :] = True
    combined.arrays["is_guest"] = is_guest
    return combined


def write_blocking_spheres_file(path, spheres):
    """
    Write RASPA-style blocking list: first line = count, then ``sx sy sz r``
    per line (fractional coordinates, radius in Å).

    Parameters
    ----------
    path : str or pathlib.Path
    spheres : sequence of (sx, sy, sz, r)
    """
    from pathlib import Path

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sph = list(spheres)
    with path.open("w") as f:
        f.write(f"{len(sph)}\n")
        for row in sph:
            sx, sy, sz, r = (float(row[0]), float(row[1]), float(row[2]), float(row[3]))
            f.write(f"{sx: .8f} {sy: .8f} {sz: .8f} {r: .8f}\n")

#
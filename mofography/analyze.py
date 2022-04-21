"""
 Functions for analyzing the chemical structure in mofography
  █████╗ ███╗   ██╗ █████╗ ██╗  ██╗   ██╗███████╗███████╗
 ██╔══██╗████╗  ██║██╔══██╗██║  ╚██╗ ██╔╝╚══███╔╝██╔════╝
 ███████║██╔██╗ ██║███████║██║   ╚████╔╝   ███╔╝ █████╗  
 ██╔══██║██║╚██╗██║██╔══██║██║    ╚██╔╝   ███╔╝  ██╔══╝  
 ██║  ██║██║ ╚████║██║  ██║███████╗██║   ███████╗███████╗
 ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝╚═╝   ╚══════╝╚══════╝
                                                         
 
"""
def detect_smiles_in_mol(smiles_molecule, smiles_pattern):
    """
    It takes a SMILES string for a molecule and a SMILES string for a pattern, and returns a list of
    dictionaries that map the nodes of the pattern to the nodes of the molecule
    
    :param smiles_molecule: The SMILES string of the molecule you want to search in
    :param smiles_pattern: The SMILES pattern you want to search for
    :return: A list of dictionaries. Each dictionary contains the mapping of the nodes in the pattern to
    the nodes in the molecule.
    """
    import pysmiles as ps
    from networkx.algorithms.isomorphism import ISMAGS, numerical_edge_match, categorical_node_match

    pattern = ps.read_smiles(smiles_pattern)
    molecule = ps.read_smiles(smiles_molecule)
    isomo = ISMAGS(molecule, pattern, node_match =categorical_node_match('element', default='none'), edge_match=numerical_edge_match('order', default='100') )
    return list(isomo.subgraph_isomorphisms_iter())

def find_patterns_in_mofid(mofid, pattern_smiles):
    """
    It takes a MOFid string and a list of SMILES strings, and returns a list of the number of times each
    SMILES string is found in the MOFid string
    
    :param mofid: the MOFid string
    :param pattern_smiles: a list of SMILES strings that you want to search for in the MOFid
    :return: A list of lists of the patterns found in the MOFid.
    """
    import numpy as np
    mofid_parts = mofid.split('.')
    
    return [detect_smiles_in_mol(mofid_parts[i].split(' ')[0], pattern_smiles) if 'MOFid' in mofid_parts[i] else detect_smiles_in_mol(mofid_parts[i], pattern_smiles) for i in range(1, len(mofid_parts)-2)]


def find_covalent_bonds_nl(data, nlist):
    """
    It takes the data and the neighbor list and returns the bonds and bond orders. 
    
    The function is a bit long, but it's not too complicated. 
    
    The first thing we do is import the necessary modules. 
    
    Then we create a molecule object from the data. 
    
    We then create a list of bonds and bond orders. 
    
    We then create a list of covalent radii. 
    
    We then loop over the neighbor list. 
    
    We then check if the bond is covalent. 
    
    If it is, we add it to the list of bonds and bond orders. 
    
    We then return the bonds and bond orders. 
    
    Let's test it out.
    
    :param data: the atoms object
    :param nlist: a list of nearest neighbors, which is a list of tuples. Each tuple is a pair of
    indices of atoms
    :return: the bonds and bond orders.
    """
    import itertools
    from pymatgen.core.bonds import CovalentBond
    from pymatgen.io.ase import AseAtomsAdaptor 
    from ase.data import covalent_radii
    import numpy as np

    molecule = AseAtomsAdaptor.get_molecule(data)
    number_of_sites = len(molecule.sites)
    bonds = []
    bond_order=[]
    radii = covalent_radii[data.get_atomic_numbers()]
    frac = molecule.frac_coords
    for (a, b) in nlist:
        # First we check if the nearest image b from a is in the unit cell
        # if not we will add that site, similar to vesta
        
        
        
        cov1 = radii[a]
        cov2 = radii[b]
        bond_test = CovalentBond(molecule.sites[a], molecule.sites[b])
        
        
        if CovalentBond.is_bonded(molecule.sites[a], molecule.sites[b], default_bl=cov1+cov2, tol=0.1):
            bonds.append([a, b])
            bond_order.append(bond_test.get_bond_order(default_bl=cov1+cov2, tol=0.1))
            # print(bonds, bond_order)
    return np.array(bonds), np.array(bond_order)
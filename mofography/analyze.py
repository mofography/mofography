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
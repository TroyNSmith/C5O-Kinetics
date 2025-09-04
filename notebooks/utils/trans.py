from protomol import rd


def intra_proton_transfer(smiles: str, idx1: int, idx2: int):
    mol = rd.mol.from_smiles(smiles)  # Initialize mol from SMILES string
    mol = rd.mol.with_coordinates(mol)  # Assign initial coordinates to mol
    bmat = rd.mol.dg_bounds_set_dist(
        mol, idx1, idx2
    )  # Calculate a new coordinates matrix with the updated idx1:idx2 distance
    return rd.mol.with_coordinates(mol, bmat=bmat)  # Assign new coordinates to output

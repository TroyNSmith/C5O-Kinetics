from rdkit import Chem
from rdkit.Chem import AllChem


def map_atomic_symbols_to_indices(smiles: str, return_tuple: bool = True) -> int:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    indices = [atom.GetIdx() for atom in mol.GetAtoms()]
    if return_tuple:
        identifiers = [f"{atom.GetSymbol()}:{atom.GetIdx()}" for atom in mol.GetAtoms()]
        tuples = [(identifiers[i], indices[i]) for i in range(len(indices))]
        return tuples

    return indices

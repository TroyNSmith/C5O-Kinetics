from rdkit import Chem
from rdkit.Chem import Descriptors


def b_cleavage(mol: Chem.Mol):
    """
    Identifies beta cleavage mechanisms in a molecule represented by an RDKit Mol object.
    Returns a list of [transition_smiles, stable_smiles, radical_smiles] or None if no cleavage is possible.
    """
    results = []
    structs = [mol]
    if Descriptors.NumRadicalElectrons(mol) == 0 or mol.GetNumAtoms() < 3:
        return None

    for initial_mol in structs:
        for atom in initial_mol.GetAtoms():
            n_rads = atom.GetNumRadicalElectrons()
            if n_rads > 0:
                r_idx = atom.GetIdx()
                alphas = [a for a in atom.GetNeighbors() if a.GetAtomicNum() > 1]
                for alpha in alphas:
                    a_idx = alpha.GetIdx()
                    betas = [
                        b
                        for b in alpha.GetNeighbors()
                        if b.GetIdx() != r_idx and b.GetAtomicNum() > 1
                    ]
                    for beta in betas:
                        b_idx = beta.GetIdx()
                        # --- Create editable mol and identify critical bonds
                        rw_mol = Chem.RWMol(initial_mol)
                        r_a_bond = rw_mol.GetBondBetweenAtoms(r_idx, a_idx)
                        a_b_bond = initial_mol.GetBondBetweenAtoms(a_idx, b_idx)
                        # --- Check for resonance
                        if a_b_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            results.append(["Check for resonance", None, None])
                            continue

                        # --- Break the bond between the alpha and beta atoms
                        trans = initial_mol.__copy__()
                        a_b_bond = trans.GetBondBetweenAtoms(a_idx, b_idx)
                        a_b_bond.SetBondType(Chem.rdchem.BondType.ZERO)
                        # --- Export the transition state
                        trans_smiles = Chem.MolToSmiles(trans)
                        # --- Simulate alpha-beta cleavage
                        rw_mol.RemoveBond(a_idx, b_idx)
                        # Remove one radical electron from the radical atom
                        rw_mol.GetAtomWithIdx(r_idx).SetNumRadicalElectrons(0)
                        # Add a radical electron to the beta atom
                        rw_mol.GetAtomWithIdx(b_idx).SetNumRadicalElectrons(
                            rw_mol.GetAtomWithIdx(b_idx).GetNumRadicalElectrons() + 1
                        )
                        # Change the rad-alpha bond from single to double (if possible)
                        if (
                            r_a_bond is not None
                            and r_a_bond.GetBondType() == Chem.rdchem.BondType.SINGLE
                        ):
                            rw_mol.RemoveBond(r_idx, a_idx)
                            rw_mol.AddBond(r_idx, a_idx, Chem.rdchem.BondType.DOUBLE)
                        elif (
                            r_a_bond is not None
                            and r_a_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                        ):
                            rw_mol.RemoveBond(r_idx, a_idx)
                            rw_mol.AddBond(r_idx, a_idx, Chem.rdchem.BondType.TRIPLE)

                        # Sanitize and get products
                        try:
                            Chem.SanitizeMol(rw_mol)
                            prod_smiles = Chem.MolToSmiles(rw_mol).split(".")
                        except Exception:
                            continue

                        # Assign stable/radical products
                        stable_smiles, radical_smiles = None, None
                        for s in prod_smiles:
                            m = Chem.MolFromSmiles(s)
                            if m and Descriptors.NumRadicalElectrons(m) > 0:
                                radical_smiles = s
                            else:
                                stable_smiles = s

                        results.append([trans_smiles, stable_smiles, radical_smiles])

    return results if results else None


def iterate(smiles, mechanism):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string provided.")

    cleavage = b_cleavage(mol)
    if cleavage is None:
        # End product reached, return a copy of the mechanism path
        return [mechanism]

    results = []
    for trans, stable, radical in cleavage:
        if radical is not None:
            # Continue the iteration with the new radical
            results.extend(iterate(radical, mechanism + [trans, stable, radical]))
        else:
            # No further radical, just append the current path
            results.append(mechanism + [trans, stable, radical])
    return results

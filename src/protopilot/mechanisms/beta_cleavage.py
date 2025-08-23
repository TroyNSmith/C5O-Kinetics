from IPython.display import display
from rdkit import Chem
from rdkit.Chem import Draw


def b_cleavage(initial_mol: Chem.Mol) -> Chem.Mol:
    """
    Identifies beta cleavage mechanisms in a molecule represented by a SMILES string.
    """
    for atom in initial_mol.GetAtoms():
        n_rads = atom.GetNumRadicalElectrons()  # - Count radical electrons
        if n_rads > 0:
            r_idx = atom.GetIdx()
            rad = initial_mol.GetAtomWithIdx(r_idx)  # - Fetch the radical atom
            alphas = rad.GetNeighbors()  # ----- Fetch neighboring (heavy) atoms
            for alpha in alphas:
                smiles = [Chem.MolToSmiles(initial_mol)]
                mols = [initial_mol]
                a_idx = alpha.GetIdx()  # - Fetch the neighboring (beta) atoms
                betas = alpha.GetNeighbors()
                for beta in betas:
                    b_idx = beta.GetIdx()
                    if b_idx != r_idx and b_idx != a_idx:
                        a_b_bond = initial_mol.GetBondBetweenAtoms(a_idx, b_idx)
                        if (
                            a_b_bond is not None
                            and a_b_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                        ):  # --- Check if the alpha-beta bond is a double bond; if so, skip
                            print("Implement check for resonance here.")
                            continue
                trans = Chem.RWMol(initial_mol.__copy__())
                mols.append(trans)
                # --- Change the bond type between the radical and alpha atom
                r_a_bond = trans.GetBondBetweenAtoms(r_idx, a_idx)
                if r_a_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    r_a_bond.SetBondType(Chem.rdchem.BondType.TRIPLE)
                else:
                    r_a_bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
                # --- Break the bond between the alpha and beta atoms
                a_b_bond = trans.GetBondBetweenAtoms(a_idx, b_idx)
                a_b_bond.SetBondType(Chem.rdchem.BondType.ZERO)
                # --- Give the beta atom a radical
                beta_atom = trans.GetAtomWithIdx(b_idx)
                beta_atom.SetNumRadicalElectrons(beta_atom.GetNumRadicalElectrons() + 1)
                # --- Export the transition state
                trans_smiles = Chem.MolToSmiles(trans)
                smiles.append(trans_smiles)
                mols.append(Chem.MolFromSmiles(trans_smiles))
                # --- Remove the transition bond and export the products
                trans.RemoveBond(a_idx, b_idx)
                prods_smiles = Chem.MolToSmiles(trans).split(".")
                for smile in prods_smiles:
                    smiles.append(smile)
                    mols.append(Chem.MolFromSmiles(smile))

                # Combine reactant and products for drawing
                img = Draw.MolsToGridImage(
                    mols,
                    molsPerRow=len(mols),
                    subImgSize=(200, 200),
                    legends=[f"{smile}" for smile in smiles],
                    highlightAtomLists=[
                        [
                            atom.GetIdx()
                            for atom in m.GetAtoms()
                            if atom.GetNumRadicalElectrons() > 0
                        ]
                        for m in mols
                    ],
                    highlightAtomColors=[
                        {
                            atom.GetIdx(): (1.0, 0.0, 0.0)
                            for atom in m.GetAtoms()
                            if atom.GetNumRadicalElectrons() > 0
                        }
                        for m in mols
                    ],
                    highlightAtomRadii=[
                        {
                            atom.GetIdx(): 0.6
                            for atom in m.GetAtoms()
                            if atom.GetNumRadicalElectrons() > 0
                        }
                        for m in mols
                    ],
                )
                label = f"{smiles[0]} → {smiles[1]} → {' + '.join(prods_smiles)}"
                print(label)
                display(img)

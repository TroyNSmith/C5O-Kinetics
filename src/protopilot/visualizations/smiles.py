from IPython.display import display
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


def visualize_smiles(smiles: str):
    """
    Visualizes a given SMILES string.

    Parameters:
    smiles (str): The SMILES representation of the molecule to visualize.

    Returns:
    None: Displays the molecule image.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    xyz_block = Chem.MolToXYZBlock(mol)

    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(xyz_block, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    natoms = mol.GetNumAtoms()
    for i in range(natoms):
        viewer.addLabel(
            i,
            {
                "backgroundOpacity": 0,
                "fontColor": "blue",
                "alignment": "center",
                "inFront": True,
            },
            {"index": i},
        )
    viewer.zoomTo()
    viewer.show()


def visualize_b_cleavage(smiles: str):
    """
    Visualizes the beta cleavage mechanism for a given SMILES string.

    Parameters:
    smiles (str): The SMILES representation of the molecule to visualize.

    Returns:
    None: Displays the reaction images and labels.
    """
    mol = Chem.MolFromSmiles(smiles)  # - Convert SMILES to RDKit molecule

    mechanisms = []
    for atom in mol.GetAtoms():
        n_rads = atom.GetNumRadicalElectrons()  # - Count radical electrons
        if n_rads > 0:
            r_idx = atom.GetIdx()
            rad = mol.GetAtomWithIdx(r_idx)  # - Fetch the radical atom
            alphas = rad.GetNeighbors()  # ----- Fetch neighboring (heavy) atoms
            for alpha in alphas:
                a_idx = alpha.GetIdx()  # - Fetch the neighboring (beta) atoms
                betas = alpha.GetNeighbors()
                for beta in betas:
                    b_idx = beta.GetIdx()
                    if b_idx != r_idx and b_idx != a_idx:
                        r_a_bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
                        if (
                            r_a_bond is not None
                            and r_a_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                        ):  # --- Check if the alpha-beta bond is a double bond; if so, skip
                            continue
                        mechanisms.append((r_idx, a_idx, b_idx))

    rxn_images = []
    rxn_labels = []

    for mechanism in mechanisms:
        r_idx, a_idx, b_idx = mechanism
        prod = Chem.RWMol(mol.__copy__())
        # --- Change the bond type between the radical and alpha atom
        r_a_bond = prod.GetBondBetweenAtoms(r_idx, a_idx)
        if r_a_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            r_a_bond.SetBondType(Chem.rdchem.BondType.TRIPLE)
        else:
            r_a_bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
        # --- Break the bond between the alpha and beta atoms
        a_b_bond = prod.GetBondBetweenAtoms(a_idx, b_idx)
        a_b_bond.SetBondType(Chem.rdchem.BondType.ZERO)

        # --- Give the beta atom a radical
        beta_atom = prod.GetAtomWithIdx(b_idx)
        beta_atom.SetNumRadicalElectrons(beta_atom.GetNumRadicalElectrons() + 1)

        prod.CommitBatchEdit()
        Chem.SanitizeMol(prod)
        smi = Chem.MolToSmiles(prod)
        split = smi.split(".")

        # --- Prepare reactants and products
        reactant = Chem.MolFromSmiles(smiles)
        products = [Chem.MolFromSmiles(p) for p in split if p]

        # --- Make label
        label = f"{smiles} → B {mechanism[1]} {mechanism[2]} → {' + '.join(split)}"
        rxn_labels.append(label)

        # --- Highlight radical atoms in red and increase their size
        highlight_atoms = []
        highlight_atom_radii = {}
        highlight_atom_colors = {}

        for m in [reactant] + products:
            for atom in m.GetAtoms():
                if atom.GetNumRadicalElectrons() > 0:
                    idx = atom.GetIdx()
                    highlight_atoms.append(idx)
                    highlight_atom_radii[idx] = 0.3
                    highlight_atom_colors[idx] = (1.0, 0.0, 0.0)  # Red

        # --- Drawing
        mols = [reactant] + products
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=len(mols),
            subImgSize=(200, 200),
            legends=["Reactant"] + [f"Product {i + 1}" for i in range(len(products))],
        )
        rxn_images.append(img)

    for img, label in zip(rxn_images, rxn_labels):
        print(label)
        display(img)

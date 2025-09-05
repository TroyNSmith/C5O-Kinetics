import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem


def _get_radical_indices(
    smiles: str = None, xyz_block: str = None, mol: Chem.Mol = None
):
    if smiles is not None:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol)
    elif xyz_block is not None:
        mol = Chem.MolFromXYZBlock(xyz_block)

    radicals = [a.GetIdx() for a in mol.GetAtoms() if a.GetNumRadicalElectrons() != 0]
    return radicals


def set_viewer(
    xyz_block: str,
    add_labels: bool = False,
    radicals: list = [],
) -> py3Dmol.view:
    """
    Return a py3Dmol viewer from an xyz_block.
    """
    viewer = py3Dmol.view()
    viewer.addModel(xyz_block, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    if add_labels:
        for i in range(int(xyz_block.split("\n")[0])):
            color = "red" if i in radicals else "blue"
            viewer.addLabel(
                str(i),
                {
                    "backgroundOpacity": 0,
                    "fontColor": color,
                    "alignment": "center",
                    "inFront": True,
                },
                {"index": i},
            )
    return viewer.zoomTo()


def from_smiles(smiles: str, add_labels: bool = False):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    radicals = _get_radical_indices(smiles=smiles)
    viewer = set_viewer(Chem.MolToXYZBlock(mol), add_labels, radicals=radicals)
    viewer.zoomTo()
    viewer.show()


def from_mol(mol: Chem.Mol, add_labels: bool = False):
    radicals = _get_radical_indices(mol=mol)
    viewer = set_viewer(Chem.MolToXYZBlock(mol), add_labels, radicals=radicals)
    viewer.zoomTo()
    viewer.show()


def from_xyz_block(xyz_block: str, add_labels: bool = False):
    radicals = _get_radical_indices(xyz_block=xyz_block)
    viewer = set_viewer(xyz_block, add_labels, radicals=radicals)
    viewer.zoomTo()
    viewer.show()


def from_allxyz_block(allxyz_block: str):
    viewer = py3Dmol.view()
    allxyz_block = "\n".join(
        line for line in allxyz_block.splitlines() if ">" not in line
    )
    viewer.addModelsAsFrames(allxyz_block, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    viewer.zoomTo()
    viewer.animate({"loop": "backAndForth"})
    viewer.show()


def visualize_traj(xyz_path, width: int = 600, height: int = 400):
    with open(xyz_path, "r") as f:
        xyz_text = "".join(line for line in f if ">" not in line)
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModelsAsFrames(xyz_text, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    viewer.zoomTo()
    viewer.animate({"loop": "backAndForth"})
    viewer.show()

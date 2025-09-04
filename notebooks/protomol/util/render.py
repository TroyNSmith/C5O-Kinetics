import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem


def set_viewer(
    xyz_block: str, add_labels: bool = False, width: int = 350, height: int = 350
) -> py3Dmol.view:
    """
    Return a py3Dmol viewer from an xyz_block.
    """
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(xyz_block, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    if add_labels:
        for i in range(int(xyz_block.split("\n")[0])):
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
    return viewer.zoomTo()


def from_smiles(
    smiles: str, add_labels: bool = False, width: int = 600, height: int = 400
):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    viewer = set_viewer(Chem.MolToXYZBlock(mol), add_labels)
    viewer.show()


def from_mol(
    mol: Chem.Mol, add_labels: bool = False, width: int = 600, height: int = 400
):
    viewer = set_viewer(Chem.MolToXYZBlock(mol), add_labels)
    viewer.show()


def from_xyz_block(
    xyz_block: str, add_labels: bool = False, width: int = 600, height: int = 400
):
    viewer = set_viewer(xyz_block, add_labels)
    viewer.show()


def from_allxyz_block(allxyz_block):
    viewer = py3Dmol.view()
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

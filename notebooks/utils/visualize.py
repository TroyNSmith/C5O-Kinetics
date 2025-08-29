import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

from pathlib import Path


def visualize_smiles(smiles: str, width: int = 600, height: int = 350):
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

    viewer = py3Dmol.view(width=width, height=height)
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


def visualize_xyz(xyz_file: Path, width: int = 600, height: int = 350):
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(Path(xyz_file).read_text(), "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    natoms = int(Path(xyz_file).read_text().split()[0])
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


def visualize_traj(xyz_path, width: int = 600, height: int = 350):
    with open(xyz_path, "r") as f:
        xyz_text = "".join(line for line in f if ">" not in line)
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModelsAsFrames(xyz_text, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    viewer.zoomTo()
    viewer.animate({"loop": "backAndForth"})
    viewer.show()

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem


def _set_viewer(
    xyz_block: str,
    add_labels: bool = False,
) -> py3Dmol.view:
    """
    Return a py3Dmol viewer from an xyz_block.
    """
    viewer = py3Dmol.view()
    viewer.addModel(xyz_block, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    if add_labels:
        for i in range(int(xyz_block.split("\n")[0])):
            color = "blue"
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
    viewer = _set_viewer(Chem.MolToXYZBlock(mol), add_labels)
    viewer.show()


def from_mol(mol: Chem.Mol, add_labels: bool = False):
    viewer = _set_viewer(Chem.MolToXYZBlock(mol), add_labels)
    viewer.show()


def from_xyz_block(xyz_block: str):
    viewer = _set_viewer(xyz_block)
    natoms = int(xyz_block.splitlines()[0])
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
    viewer.show()


def visualize_traj(xyz_block: str):
    xyz_text = xyz_block.replace(">\n", "")
    viewer = py3Dmol.view()
    viewer.addModelsAsFrames(xyz_text, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    viewer.zoomTo()
    viewer.animate({"loop": "backAndForth"})
    viewer.show()

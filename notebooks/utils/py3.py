import py3Dmol


def render_xyz(xyz_block: str, width: int = 600, height: int = 400):
    """
    Visualizes a given SMILES string.

    Parameters:
    smiles (str): The SMILES representation of the molecule to visualize.

    Returns:
    None: Displays the molecule image.
    """
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(xyz_block, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    natoms = xyz_block.split("\n")[0]
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

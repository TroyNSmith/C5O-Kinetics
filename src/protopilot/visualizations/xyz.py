import py3Dmol

from pathlib import Path


def labelled_3d_structure(xyz_file: str | Path):
    natoms = int(Path(xyz_file).read_text().split()[0])

    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(Path(xyz_file).read_text(), "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})

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


def animate_xyz(xyz_path):
    with open(xyz_path, "r") as f:
        xyz_text = "".join(line for line in f if ">" not in line)
    natoms = int(xyz_text.split()[0])
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModelsAsFrames(xyz_text, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
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
    viewer.animate({"loop": "backAndForth"})
    viewer.show()

import altair as alt
import pandas as pd
import py3Dmol
import pyparsing as pp
from pyparsing import pyparsing_common as ppc

from pathlib import Path
import re


def labelled_3d_structure(xyz_file: str):
    viewer = py3Dmol.view(width=400, height=400)
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


def find_xyz(smiles_str: str):
    identifier = (
        smiles_str.replace("[", "_")
        .replace("]", "_")
        .replace("(", "ch_")
        .replace(")", "_hc")
        .replace("=", "dbl")
        .replace("#", "trpl")
        .replace("/", "up")
        .replace("\\", "dwn")
    )
    revdsd_dir = Path.home() / f"C5O-Kinetics/calc/{identifier}/Optimization/run/Guess"

    for file in revdsd_dir.rglob("guess.xyz"):
        return file


def find_dirs(smiles_str: str):
    identifier = (
        smiles_str.replace("[", "_")
        .replace("]", "_")
        .replace("(", "ch_")
        .replace(")", "_hc")
        .replace("=", "dbl")
        .replace("#", "trpl")
        .replace("/", "up")
        .replace("\\", "dwn")
    )
    par_dir = Path.home() / f"C5O-Kinetics/calc/{identifier}"
    output = []
    for dir in par_dir.iterdir():
        output.append((dir.name, dir))
    return output


def allxyz_energies(scan_dir: str | Path):
    allxyz_paths = [Path(x) for x in scan_dir.rglob("REVDSD.allxyz")]

    comment = pp.Group(
        pp.Keyword("Scan Step")
        + ppc.integer("index")
        + pp.Keyword("E")
        + ppc.fnumber("energy")
    )
    expr = pp.OneOrMore(pp.SkipTo(comment, include=True))
    text = re.sub(">\n", "", allxyz_paths[0].read_text())

    indices = []
    energies = []
    results = expr.parse_string(text)
    for result in results[1::2]:
        index, energy = result[1:4:2]
        indices.append(index)
        energies.append(energy)

    df = pd.DataFrame(
        {
            "index": indices,
            "energy": energies,
        }
    )
    alt.Chart(df).mark_point().encode(
        x="index",
        y=alt.Y("energy", scale=alt.Scale(zero=False)),
    ).show()

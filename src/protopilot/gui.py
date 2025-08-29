"""
Backend functions for the ipywidgets-based GUI in gui.ipynb.
This module provides functions for:
- Parsing SMILES and generating 3D structures
- Extracting ZPV and SPC energies
- Rendering xyz, trajectories, and vibrational modes
- Managing method options and data
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from pathlib import Path

# Example: Generate 3D structure from SMILES


def smiles_to_xyz_block(smiles: str, method: str = "Optimization") -> str:
    dir_name = (
        smiles.replace("[", "_")
        .replace("]", "_")
        .replace("(", "ch_")
        .replace(")", "_hc")
        .replace("=", "dbl")
        .replace("#", "trpl")
        .replace("/", "up")
        .replace("\\", "dwn")
    )
    calc_dir = Path.home() / f"C5O-Kinetics/calc/{dir_name}/{method}/run/REVDSD"
    for dir in calc_dir.iterdir():
        print(dir)
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    return Chem.MolToXYZBlock(mol)


# Example: Extract ZPV and SPC energies from log files (stub)
def extract_energies(calc_dir: Path) -> tuple[str, str]:
    # TODO: Implement real file parsing
    zpv = "-123.45"  # Placeholder
    spc = "-456.78"  # Placeholder
    return zpv, spc


# Example: Get available methods (stub)
def get_methods(calc_dir: Path) -> list[str]:
    # TODO: Implement real method discovery
    return ["Method1", "Method2"]


# Example: Get trajectory data (stub)
def get_trajectory(calc_dir: Path) -> pd.DataFrame:
    # TODO: Implement real trajectory extraction
    return pd.DataFrame({"step": [0, 1, 2], "energy": [-123.4, -122.8, -123.1]})


# Example: Get vibrational mode data (stub)
def get_vibrational_mode(calc_dir: Path):
    # TODO: Implement real vibrational mode extraction
    return None

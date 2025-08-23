from pathlib import Path


def main(smiles_string: str):
    # --- Clean the name of the smiles string
    name = smiles_string.replace("[", "_").replace("]", "_").replace("=", "dbl")

    # --- Define methods and method directory labels
    methods = [
        "XTB",
        "REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c",
        "CCSD(T)-F12 cc-pVTZ-F12 cc-pVTZ-F12-CABS",
    ]
    labels = ["XTB", "REVDSD", "CCSDT"]

    # --- Create directories for calculations
    par_dir = Path(f"calc/{name}").resolve()
    for label in labels:
        method_dir = par_dir / label
        method_dir.mkdir(parents=True, exist_ok=True)

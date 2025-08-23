from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from pathlib import Path
import shutil
import sys

sys.path.append(str(Path.home() / "C5O-Kinetics"))

from protopilot.io.utils import orca_input_and_sh


def write(
    smiles_string: str = None,
    initial_xyz: Path = None,
    transition_state: bool = False,
    multiplicity: int = None,
    charge: int = 0,
):
    """
    If smiles_string is provided, optimize with XTB, then run REVDSD and CCSDT.
    If initial_xyz is provided, skip XTB and run REVDSD and CCSDT.
    """
    labels = ["XTB", "REVDSD", "CCSDT"]

    if smiles_string and initial_xyz:
        raise ValueError("Provide either smiles_string or initial_xyz, not both.")
    elif smiles_string:
        name = (
            smiles_string.replace("[", "_")
            .replace("]", "_")
            .replace("(", "ch_")
            .replace(")", "_hc")
            .replace("=", "dbl")
            .replace("#", "trpl")
        )
        par_dir = Path.home() / f"C5O-Kinetics/calc/{name}/Optimization/run"
        if par_dir.exists():
            raise FileExistsError(f"Directory {par_dir} already exists.")
        initial_dir = par_dir / "Guess"
        initial_dir.mkdir(parents=True, exist_ok=True)
    elif initial_xyz:
        name = initial_xyz.stem
        labels.remove("XTB")
        par_dir = initial_xyz.parent / "run"
        if multiplicity is None:
            raise ValueError("Multiplicity must be provided when using initial_xyz.")
    else:
        raise ValueError("Either smiles_string or initial_xyz must be provided.")

    for label in labels:
        method_dir = par_dir / label
        method_dir.mkdir(parents=True, exist_ok=True)

    # --- Prepare molecule and properties
    if smiles_string:
        mol_format = Chem.AddHs(Chem.MolFromSmiles(smiles_string))
        AllChem.EmbedMolecule(mol_format)
        charge = Chem.GetFormalCharge(mol_format)
        multiplicity = Descriptors.NumRadicalElectrons(mol_format) + 1
        # Write guess xyz
        init_out = initial_dir / "guess.xyz"
        xyz = Chem.MolToXYZBlock(mol_format)
        init_out.write_text(xyz.lstrip())
        # Copy xyz to XTB
        shutil.copyfile(init_out, par_dir / labels[0] / "guess.xyz")
        # Write XTB input and submission script
        xtb_job_opts = (
            "# --- Run XTB and copy output\n"
            f"$(which orca) {labels[0]}.inp > {labels[0]}.log\n"
            f"cp {labels[0]}.globalminimum.xyz ../{labels[1]}/{labels[0]}.xyz\n\n"
        )
        orca_input_and_sh(
            method="XTB",
            job_opts=xtb_job_opts,
            smiles=name,
            charge=charge,
            multiplicity=multiplicity,
            output_dir=par_dir / "XTB",
            xyz_name="guess",
            job_type="goat",
        )
        # Write REVDSD input (uses XTB output)
        revdsd_job_opts = (
            "# --- Run REVDSD and copy output\n"
            f"$(which orca) {labels[1]}.inp > {labels[1]}.log\n"
            f"cp {labels[1]}.xyz ../{labels[2]}/\n\n"
        )
        orca_input_and_sh(
            method="REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c",
            job_opts=revdsd_job_opts,
            smiles=name,
            charge=charge,
            multiplicity=multiplicity,
            output_dir=par_dir / "REVDSD",
            xyz_name="XTB",
            job_type="opt_freq",
        )

        ccsdt_job_opts = (
            f"# --- Run CCSDT\n$(which orca) {labels[2]}.inp > {labels[2]}.log\n\n"
        )
        # Write CCSDT input (uses REVDSD output)
        orca_input_and_sh(
            method="CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/c",
            job_opts=ccsdt_job_opts,
            smiles=name,
            charge=charge,
            multiplicity=multiplicity,
            output_dir=par_dir / "CCSDT",
            xyz_name="REVDSD",
            job_type="elec",
        )
        # Write parent submission script
        sh_script = (
            "#!/bin/bash\n\n"
            f"cd {par_dir}/XTB\n"
            "xtb_id=$(sbatch submit_XTB.sh | awk '{print $NF}')\n\n"
            f"cd {par_dir}/REVDSD\n"
            "revdsd_id=$(sbatch --dependency=afterok:$xtb_id submit_REVDSD.sh | awk '{print $NF}')\n\n"
            f"cd {par_dir}/CCSDT\n"
            "sbatch --dependency=afterok:$revdsd_id submit_CCSDT.sh"
        )
        path_out = Path(par_dir) / "submit.sh"
        path_out.write_text(sh_script)

        print(f"bash {path_out}")

    elif initial_xyz:
        raise NotImplementedError("Initial XYZ input not yet implemented.")
    else:
        raise ValueError("Either smiles_string or initial_xyz must be provided.")

    return par_dir

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from pathlib import Path

from ..io import intermediate_writing as writing

methods = [
    "XTB",
    "REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c",
    "CCSD(T)-F12 cc-pVTZ-F12 cc-pVTZ-F12-CABS",
]


def clean_method(method):
    return (
        method.replace(" ", "_")
        .replace("/", "")
        .replace(")", "")
        .replace("(", "")
        .replace("-", "_")
    )


def main(smiles_string: str):
    name = smiles_string.replace("[", "_").replace("]", "_").replace("=", "dbl")
    valid_names = [clean_method(m) for m in methods]
    par_dir = Path(f"calc/{name}").resolve()
    work_dirs = {m: par_dir / Path(f"{clean_method(m)}/run") for m in methods}
    for m in methods:
        Path(work_dirs[m]).mkdir(parents=True, exist_ok=True)

    # Prepare molecule and properties
    mol_format = Chem.AddHs(Chem.MolFromSmiles(smiles_string))
    AllChem.EmbedMolecule(mol_format)
    charge = Chem.GetFormalCharge(mol_format)
    multiplicity = Descriptors.NumRadicalElectrons(mol_format) + 1

    # Write initial guess xyz for the first method if not already present
    guess_xyz_path = work_dirs[methods[0]] / f"{valid_names[0]}.xyz"
    if not guess_xyz_path.exists():
        writing.guess_xyz(
            mol=mol_format,
            output_dir=work_dirs[methods[0]],
            xyz_name=valid_names[0],
        )

    # Write input files for all methods, skipping if output already exists
    methods_to_run = []
    for i, method in enumerate(methods):
        method_dir = work_dirs[method]
        log_file = method_dir / f"{valid_names[i]}.log"
        is_complete = False
        if log_file.exists() and log_file.stat().st_size > 0:
            try:
                with open(log_file, "r") as f:
                    lines = f.readlines()
                    for line in lines[-10:]:
                        if "****ORCA TERMINATED NORMALLY****" in line:
                            is_complete = True
                            break
            except Exception as e:
                print(f"Warning: Could not read {log_file}: {e}")
        if is_complete:
            print(f"Skipping {method}: log file exists and terminated normally.")
            continue
        methods_to_run.append(method)
        if i == 0:
            job_type = "goat"
            xyz_name = valid_names[0]
        elif i == len(methods) - 2:
            job_type = "opt_freq"
            xyz_name = valid_names[i - 1]
        elif i == len(methods) - 1:
            job_type = "elec"
            xyz_name = valid_names[i - 1]
        else:
            job_type = "opt"
            xyz_name = valid_names[i - 1]
        writing.write_orca_input(
            method=method,
            charge=charge,
            multiplicity=multiplicity,
            output_dir=method_dir,
            xyz_name=xyz_name,
            job_type=job_type,
        )

    # Write a parent .sh to orchestrate only the calculations that need to be run
    if methods_to_run:
        writing.sub_bash(
            name=name,
            methods=[clean_method(m) for m in methods_to_run],
            work_dirs={m: work_dirs[m] for m in methods_to_run},
            par_dir=par_dir,
        )
    else:
        print("All steps are already complete. No bash script written.")

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from pathlib import Path
import shutil

from .writing import orca_input_and_sh


def generate_guess(smiles: str, job_type: str, label: str):
    identifier = (
        smiles.replace("[", "_")
        .replace("]", "_")
        .replace("(", "ch_")
        .replace(")", "_hc")
        .replace("=", "dbl")
        .replace("#", "trpl")
        .replace("/", "up")
        .replace("\\", "dwn")
    )
    output_dir = (
        Path.home()
        / f"C5O-Kinetics/calc/{identifier}/Optimization/run/{job_type}_{label}/init.xyz"
    )
    mol_format = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol_format)
    charge = Chem.GetFormalCharge(mol_format)
    multiplicity = Descriptors.NumRadicalElectrons(mol_format) + 1
    # Write guess xyz
    xyz = Chem.MolToXYZBlock(mol_format)
    output_dir.write_text(xyz.lstrip())


def get_xyz(smiles: str = "[CH]1CO1", method: str = "REVDSD"):
    identifier = (
        smiles.replace("[", "_")
        .replace("]", "_")
        .replace("(", "ch_")
        .replace(")", "_hc")
        .replace("=", "dbl")
        .replace("#", "trpl")
        .replace("/", "up")
        .replace("\\", "dwn")
    )
    xyz_dir = Path.home() / f"C5O-Kinetics/calc/{identifier}/Optimization/run/{method}"
    unique_files = {f.name: f for f in xyz_dir.rglob("*xyz")}
    return sorted(unique_files.items(), key=lambda x: x[0])


def get_transitions(smiles: str = "[CH]1CO1"):
    identifier = (
        smiles.replace("[", "_")
        .replace("]", "_")
        .replace("(", "ch_")
        .replace(")", "_hc")
        .replace("=", "dbl")
        .replace("#", "trpl")
        .replace("/", "up")
        .replace("\\", "dwn")
    )
    trans_dir = Path.home() / f"C5O-Kinetics/calc/{identifier}"
    return sorted(
        [(f.name, f) for f in trans_dir.iterdir() if f.is_dir()], key=lambda x: x[0]
    )


def get_scan_xyzs(trans_dir: Path):
    scan_dir = trans_dir / "run/REVDSD_Scan"
    files = [(str(f).split("REVDSD_Scan")[1], f) for f in scan_dir.rglob("*xyz")]
    return sorted(files, key=lambda x: x[0])


def get_energies(xyz_path: Path):
    run_dir = Path(str(xyz_path).split("run")[0]) / "run"
    # ZPV energies
    revdsd_logs = list((run_dir / "REVDSD").rglob("REVDSD.log"))
    zpv_energies = []
    for log in revdsd_logs:
        with open(log, "r") as f:
            zpv_energies += [
                line.split("Eh")[1].strip() for line in f if "Zero point energy" in line
            ]
    zpv_text = (
        f"ZPV: {','.join(zpv_energies)}"
        if zpv_energies
        else "Zero Point Vibrational Energy: --"
    )
    # SPC energies
    ccsdt_logs = list((run_dir / "CCSDT").rglob("CCSDT.log"))
    spc_energies = []
    for log in ccsdt_logs:
        with open(log, "r") as f:
            spc_energies += [
                line.split("ENERGY")[1].strip()
                for line in f
                if "FINAL SINGLE POINT ENERGY" in line
            ]
    spc_text = (
        f"SPC (hartrees): {','.join(spc_energies)}"
        if spc_energies
        else "Single Point Calculation Energy: --"
    )
    return zpv_text, spc_text


def orca_optimization(
    method: str,
    smiles: str,
    charge: int = 0,
    multiplicity: int = 2,
    job_type: str = "opt",
    guess: bool = False,
) -> None:
    """
    Write an ORCA input file for various calculation types.

    job_type options:
        - "opt": Optimization only (no frequency calculation) (default)
        - "opt_freq": Optimization and frequency calculation
        - "goat": Geometry optimization and frequency calculation (GOAT)
        - "elec": Single-point electronic structure calculation
    """
    identifier = (
        smiles.replace("[", "_")
        .replace("]", "_")
        .replace("(", "ch_")
        .replace(")", "_hc")
        .replace("=", "dbl")
        .replace("#", "trpl")
        .replace("/", "up")
        .replace("\\", "dwn")
    )
    if "XTB" in method:
        label = "XTB"
        job_name = f"{smiles}_xtb"
        num_cpus = 8
        mem_per_cpu = 1
        lscratch_size = 10
        time = "00:30:00"
    elif "REVDSD" in method:
        label = "REVDSD"
        job_name = f"{smiles}_revdsd"
        num_cpus = 16
        mem_per_cpu = 1
        lscratch_size = 20
        time = "8:00:00"
    elif "CCSD(T)" in method:
        label = "CCSDT"
        job_name = f"{smiles}_ccsdt"
        if smiles.count("C") + smiles.count("O") >= 8:
            num_cpus = 20
            mem_per_cpu = 12
            lscratch_size = 200
            time = "3:00:00"
        elif smiles.count("C") + smiles.count("O") == 7:
            num_cpus = 20
            mem_per_cpu = 6
            lscratch_size = 200
            time = "1:00:00"
        elif smiles.count("C") + smiles.count("O") == 6:
            num_cpus = 20
            mem_per_cpu = 4
            lscratch_size = 200
            time = "1:00:00"
        elif smiles.count("C") + smiles.count("O") == 5:
            num_cpus = 16
            mem_per_cpu = 4
            lscratch_size = 200
            time = "1:00:00"
        elif smiles.count("C") + smiles.count("O") == 4:
            num_cpus = 12
            mem_per_cpu = 4
            lscratch_size = 100
            time = "00:30:00"
        elif smiles.count("C") + smiles.count("O") == 3:
            num_cpus = 12
            mem_per_cpu = 2
            lscratch_size = 50
            time = "00:15:00"
        elif smiles.count("C") + smiles.count("O") == 2:
            num_cpus = 4
            mem_per_cpu = 1
            lscratch_size = 50
            time = "00:15:00"
        elif smiles.count("C") + smiles.count("O") == 1:
            num_cpus = 4
            mem_per_cpu = 1
            lscratch_size = 50
            time = "00:15:00"
    else:
        raise ValueError(f"Unknown method: {method}")

    # --- ORCA Input File
    par_line = "# --- ORCA Parameters\n"
    par_line += f"%PAL NPROCS {num_cpus} END\n"  # Here
    par_line += f"%maxcore {mem_per_cpu * 750}\n\n"  # Here
    job_line = "# --- Job Parameters\n"
    if job_type == "opt":
        job_line += f"! {method} OPT\n\n"
    elif job_type == "opt_freq":
        job_line += f"! {method} OPT NumFreq\n\n"
    elif job_type == "opt_freq_ts":
        job_line += f"! {method} OPTTS NumFreq\n\n%geom\n    Calc_Hess true\n    NumHess true\nend\n\n"
    elif job_type == "goat":
        job_line += f"! {method} GOAT\n\n"
    elif job_type == "elec":
        job_line += f"! {method}\n\n"
    else:
        raise ValueError(f"Unknown job_type: {job_type}")

    xyz_line = "# --- XYZ Parameters\n"
    xyz_line += f"* xyzfile {charge} {multiplicity} init.xyz\n"

    txt = f"{par_line}{job_line}{xyz_line}"
    output_dir = (
        Path.home()
        / f"C5O-Kinetics/calc/{identifier}/Optimization/run/{job_type}_{label}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    path_out = output_dir / f"{label}.inp"
    if not path_out.exists():
        path_out.write_text(txt)
    # ----- Submission Shell Script
    # --- SLURM Parameters
    slurm_opts = (
        "#!/bin/bash\n"
        f"#SBATCH --partition=batch\n"
        f"#SBATCH --gres=lscratch:{lscratch_size}\n"
        f"#SBATCH --job-name={job_name}\n"
        "#SBATCH --nodes=1\n"
        f"#SBATCH --ntasks={num_cpus}\n"
        f"#SBATCH --ntasks-per-node={num_cpus}\n"  # Here
        "#SBATCH --cpus-per-task=1\n"
        f"#SBATCH --time={time}\n"
        f"#SBATCH --mem-per-cpu={mem_per_cpu}G\n\n"  # Here
    )
    # --- ORCA Module Initialization
    module_opts = "# --- Load ORCA 6.1\n"
    module_opts += "module load ORCA/6.1\n"
    module_opts += "cd ${SLURM_SUBMIT_DIR}\n\n"
    # --- Lscratch Directory Initialization
    lscratch_opts = "# --- Initialize lscratch directory\n"
    lscratch_opts += "mkdir -p /lscratch/${USER}/${SLURM_JOB_ID}\n"
    lscratch_opts += 'cp -r "${SLURM_SUBMIT_DIR}/." /lscratch/${USER}/${SLURM_JOB_ID}\n'
    lscratch_opts += "cd /lscratch/${USER}/${SLURM_JOB_ID}\n\n"
    # --- Append Cleanup and Copy Commands
    job_opts = (
        "# --- Run orca\n"
        f"$(which orca) {label}.inp > {label}.log\n\n"
        "# --- Extract important files\n"
        "cp -r /lscratch/${USER}/${SLURM_JOB_ID} ${SLURM_SUBMIT_DIR}\n"
        "rm -rf /lscratch/${USER}/${SLURM_JOB_ID}\n"
    )
    # --- Combine All Parts and Write to File
    script = slurm_opts + module_opts + lscratch_opts + job_opts
    path_out = Path(output_dir) / f"submit_{label}.sh"
    if not path_out.exists():
        path_out.write_text(script)

    if guess is True:
        generate_guess(smiles=smiles, job_type=job_type, label=label)

    return output_dir, label, identifier


def optimization(
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
            .replace("/", "up")
            .replace("\\", "dwn")
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

        return f"bash {path_out}"

    elif initial_xyz:
        raise NotImplementedError("Initial XYZ input not yet implemented.")
    else:
        raise ValueError("Either smiles_string or initial_xyz must be provided.")

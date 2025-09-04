from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from pathlib import Path

from .resources import optimization_pars


def methods():
    methods = [
        ("NONE", None),
        ("XTB", "XTB"),
        ("REVDSD", "REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c"),
        (
            "CCSDT",
            "CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/c",
        ),
    ]
    return methods


def types():
    types = [
        ("NONE", None),
        ("GOAT", "goat"),
        ("OPT", "opt"),
        ("FREQUENCY", "zpv"),
        ("SPC", "spc"),
    ]
    return types


def smiles_identifier(smiles: str):
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
    return identifier


def generate_guess(smiles: str, job_type: str, label: str):
    identifier = smiles_identifier(smiles)
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
    identifier = smiles_identifier(smiles)
    trans_dir = Path.home() / f"C5O-Kinetics/calc/{identifier}"
    return sorted(
        [(f.name, f) for f in trans_dir.iterdir() if f.is_dir()], key=lambda x: x[0]
    )


def get_scan_xyzs(trans_dir: Path):
    scan_dir = trans_dir / "run/REVDSD_Scan"
    files = [(str(f).split("REVDSD_Scan")[1], f) for f in scan_dir.rglob("*xyz")]
    return sorted(files, key=lambda x: x[0])


def get_energies(xyz_path: Path):
    zpv_energies = []
    spc_energies = []
    for inp in list((xyz_path.parent.parent.parent).rglob("*inp")):
        with open(inp, "r") as f:
            content = f.read()
            if "GOAT" in content:
                pass
            elif "ZPV" in content:
                for log in list((inp.parent).rglob("*log")):
                    with open(log, "r") as g:
                        zpv_energies += [
                            line.split("Eh")[1].strip()
                            for line in g
                            if "Zero point energy" in line
                        ]
            elif "SPC" in content:
                for log in list((inp.parent).rglob("*log")):
                    with open(log, "r") as g:
                        spc_energies += [
                            line.split("ENERGY")[1].strip()
                            for line in g
                            if "FINAL SINGLE POINT ENERGY" in line
                        ]
    zpv_text = (
        f"ZPV: {','.join(set(zpv_energies))}"
        if zpv_energies
        else "Zero Point Vibrational Energy: --"
    )
    spc_text = (
        f"SPC (hartrees): {','.join(set(spc_energies))}"
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
    cp_cmd: str = "",
    log_cmd: str = "",
) -> None:
    """
    Write an ORCA input file for various calculation types.

    job_type options:
        - "opt": Optimization only (no frequency calculation) (default)
        - "opt_freq": Optimization and frequency calculation
        - "goat": Geometry optimization and frequency calculation (GOAT)
        - "elec": Single-point electronic structure calculation
    """
    identifier = smiles_identifier(smiles)
    label = method[:3]
    job_name = f"{smiles}_{label}"
    num_heavy_ats = smiles.count("C") + smiles.count("O")
    pars = optimization_pars(num_heavy_ats, label)
    # ----- ORCA Input File
    par_line = (
        "# --- ORCA Parameters\n"
        f"%PAL NPROCS {pars['num_cpus']} END\n"
        f"%maxcore {pars['mem_per_cpu'] * 750}\n\n"
    )
    job_line = "# --- Job Parameters\n"
    if job_type == "opt":
        job_line += f"! {method} OPT\n\n"
    elif job_type == "zpv":
        job_line += f"# --- ZPV\n! {method} OPT NumFreq\n\n"
    elif job_type == "goat":
        job_line += f"! {method} GOAT\n\n"
    elif job_type == "spc":
        job_line += f"# --- SPC\n! {method}\n\n"
    else:
        raise ValueError(f"Unknown job_type: {job_type}")
    xyz_line = f"# --- XYZ Parameters\n* xyzfile {charge} {multiplicity} init.xyz\n"
    output_dir = (
        Path.home()
        / f"C5O-Kinetics/calc/{identifier}/Optimization/run/{job_type}_{label}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    path_out = output_dir / f"{label}.inp"
    path_out.write_text(f"{par_line}{job_line}{xyz_line}")
    # ----- Submission Shell Script
    # --- SLURM Parameters
    slurm_opts = (
        "#!/bin/bash\n"
        f"#SBATCH --partition=batch\n"
        f"#SBATCH --gres=lscratch:{pars['lscratch_size']}\n"
        f"#SBATCH --job-name={job_name}\n"
        "#SBATCH --nodes=1\n"
        f"#SBATCH --ntasks={pars['num_cpus']}\n"
        f"#SBATCH --ntasks-per-node={pars['num_cpus']}\n"  # Here
        "#SBATCH --cpus-per-task=1\n"
        f"#SBATCH --time={pars['time']}\n"
        f"#SBATCH --mem-per-cpu={pars['mem_per_cpu']}G\n\n"  # Here
    )
    # --- ORCA Module Initialization
    module_opts = (
        "# --- Load ORCA 6.1\nmodule load ORCA/6.1\ncd ${SLURM_SUBMIT_DIR}\n\n"
    )
    # --- lscratch Directory Initialization
    lscratch_opts = (
        "# --- Initialize lscratch directory\n"
        "mkdir -p /lscratch/${USER}/${SLURM_JOB_ID}\n"
        'cp -r "${SLURM_SUBMIT_DIR}/." /lscratch/${USER}/${SLURM_JOB_ID}\n'
        "cd /lscratch/${USER}/${SLURM_JOB_ID}\n\n"
    )
    # --- Append Cleanup and Copy Commands
    job_opts = (
        "# --- Run orca\n"
        f"$(which orca) {label}.inp > {label}.log\n\n"
        "# --- Extract important files\n"
        f"{log_cmd}"
        f"{cp_cmd}"
        "cp -r /lscratch/${USER}/${SLURM_JOB_ID} ${SLURM_SUBMIT_DIR}\n"
        "rm -rf /lscratch/${USER}/${SLURM_JOB_ID}\n"
    )
    # --- Combine All Parts and Write to File
    script = slurm_opts + module_opts + lscratch_opts + job_opts
    path_out = Path(output_dir) / f"submit_{label}.sh"
    path_out.write_text(script)
    if guess is True:
        generate_guess(smiles=smiles, job_type=job_type, label=label)

from pathlib import Path


def orca_input_and_sh(
    method: str,
    smiles: str,
    job_opts: str,
    charge: int,
    multiplicity: int,
    output_dir: Path,
    xyz_name: str,
    job_type: str = "opt",
    rxn_coord: str = None,
) -> None:
    """
    Write an ORCA input file for various calculation types.

    job_type options:
        - "opt": Optimization only (no frequency calculation) (default)
        - "opt_freq": Optimization and frequency calculation
        - "opt_freq_ts": Optimization and frequency calculation for transition states
        - "goat": Geometry optimization and frequency calculation (GOAT)
        - "elec": Single-point electronic structure calculation
        - "rxn_coord": Reaction coordinate to scan
    """
    if "XTB" in method:
        label = "XTB"
        job_name = f"{smiles}_xtb"
        energy_cmd = "ENERGY=0\n"
        partition = "batch"
        num_cpus = 8
        mem_per_cpu = 1
        lscratch_size = 10
        time = "00:30:00"
        if job_type == "goat":
            cp_command = f"cp XTB.globalminimum.xyz {output_dir.parent}/REVDSD/XTB.xyz"
        else:
            cp_command = ""
    elif "REVDSD" in method:
        label = "REVDSD"
        job_name = f"{smiles}_revdsd"
        energy_cmd = f'ENERGY=$(grep -R --include "REVDSD.log" "Zero point energy")\n'
        partition = "batch"
        num_cpus = 16
        mem_per_cpu = 1
        lscratch_size = 20
        time = "8:00:00"
        cp_command = f"cp REVDSD.xyz {output_dir.parent}/CCSDT/REVDSD.xyz"
    elif "CCSD(T)" in method:
        label = "CCSDT"
        job_name = f"{smiles}_ccsdt"
        partition = "batch"
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
        cp_command = ""
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

    if rxn_coord:
        job_line += f"# --- Define the scanning path\n% geom\n   scan\n       {rxn_coord}\n   end\nend\n\n"

    xyz_line = "# --- XYZ Parameters\n"
    xyz_line += f"* xyzfile {charge} {multiplicity} {xyz_name}.xyz\n"

    txt = f"{par_line}{job_line}{xyz_line}"
    path_out = output_dir / f"{label}.inp"
    path_out.write_text(txt)
    # ----- Submission Shell Script
    # --- SLURM Parameters
    slurm_opts = (
        "#!/bin/bash\n"
        f"#SBATCH --partition={partition}\n"
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
    job_opts += (
        "# --- Extract important files\n"
        "cp -r /lscratch/${USER}/${SLURM_JOB_ID} ${SLURM_SUBMIT_DIR}\n"
        f"{cp_command}\n"
        "rm -rf /lscratch/${USER}/${SLURM_JOB_ID}\n\n"
    )
    log_opts = (
        "# --- Update the JSON log\n"
        'if [! -f $HOME/C5O-Kinetics/results.json]; then\n'
        '   echo "{}" > $HOME/C5O-Kinetics/results.json\n'
        "fi\n\n"
    )
    if 
    (
        f"jq --arg s {smiles} --arg m {method} --arg t {job_type}"

    )
    # --- Combine All Parts and Write to File
    script = slurm_opts + module_opts + lscratch_opts + job_opts
    path_out = Path(output_dir) / f"submit_{label}.sh"
    path_out.write_text(script)

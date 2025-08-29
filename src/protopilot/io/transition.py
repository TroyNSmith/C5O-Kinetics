from pathlib import Path
import shutil
import sys

sys.path.append(str(Path.home() / "C5O-Kinetics"))
from protopilot.io.utils import orca_input_and_sh


def write_scan(
    intermediate_xyz: Path,
    trans_smiles: str,
    scan_coordinates: str,
    rxn_coordinates: str,
    charge: int = 0,
    multiplicity: int = 2,
):
    par_dir = Path(intermediate_xyz.split("Optimization")[0])
    work_dir = par_dir / f"{scan_coordinates.replace(' ', '_')}/run"
    if work_dir.exists():
        raise FileExistsError(f"{work_dir} already exists.")
    work_dir.mkdir(parents=True, exist_ok=True)
    scan_dir = work_dir / "REVDSD_Scan"
    scan_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy(intermediate_xyz, scan_dir / "OPT.xyz")

    scan_job_opts = "# --- Run XTB\n$(which orca) REVDSD.inp > REVDSD.log\n"

    orca_input_and_sh(
        method="REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c",
        smiles=trans_smiles,
        job_opts=scan_job_opts,
        charge=charge,
        multiplicity=multiplicity,
        output_dir=scan_dir,
        xyz_name="OPT",
        job_type="opt",
        rxn_coord=rxn_coordinates,
    )


def write_opt(
    transition_xyz: Path,
    trans_smiles: str,
    charge: int = 0,
    multiplicity: int = 2,
):
    work_dir = Path(transition_xyz.split("REVDSD")[0])

    revdsd_dir = work_dir / "REVDSD"
    revdsd_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(Path(transition_xyz), revdsd_dir / "Scan.xyz")
    revdsd_job_opts = (
        "# --- Run REVDSD and copy output\n"
        "$(which orca) REVDSD.inp > REVDSD.log\n"
        "cp REVDSD.xyz ../CCSDT/\n\n"
        "# --- Extract imaginary node from REVDSD log\n"
        'imaginary_mode=$(grep -R --include "REVDSD.log" "imaginary mode")\n'
        'trimmed_str="${imaginary_mode##*([[:space:]])}"\n'
        'node="$(trimmed_str:0:1)"\n'
        "orca_pltvib REVDSD.log $node\n\n"
    )
    orca_input_and_sh(
        method="REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c",
        smiles=trans_smiles,
        job_opts=revdsd_job_opts,
        charge=charge,
        multiplicity=multiplicity,
        output_dir=revdsd_dir,
        xyz_name="Scan",
        job_type="opt_freq_ts",
    )
    # --- CCSDT
    ccsdt_dir = work_dir / "CCSDT"
    ccsdt_dir.mkdir(parents=True, exist_ok=True)
    ccsdt_job_opts = "# --- Run CCSDT\n$(which orca) CCSDT.inp > CCSDT.log\n\n"
    # Write CCSDT input (uses REVDSD output)
    orca_input_and_sh(
        method="CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/c",
        job_opts=ccsdt_job_opts,
        smiles=trans_smiles,
        charge=charge,
        multiplicity=multiplicity,
        output_dir=ccsdt_dir,
        xyz_name="REVDSD",
        job_type="elec",
    )

    # Write parent submission script
    sh_script = (
        "#!/bin/bash\n\n"
        f"cd {work_dir}/REVDSD\n"
        "revdsd_id=$(sbatch submit_REVDSD.sh | awk '{print $NF}')\n\n"
        f"cd {work_dir}/CCSDT\n"
        "sbatch --dependency=afterok:$revdsd_id submit_CCSDT.sh"
    )
    #     f"cd {work_dir}/XTB\n"
    # "xtb_id=$(sbatch submit_XTB.sh | awk '{print $NF}')\n\n"
    path_out = Path(work_dir) / "submit_opt.sh"
    path_out.write_text(sh_script)
    print(f"bash {path_out}")

from pathlib import Path

from . import io


def from_smiles(smiles: str, multiplicity: int):
    sh_script = "#!/bin/bash\n\n"
    identifier = io.smiles_identifier(smiles)
    workdir = Path.home() / f"C5O-Kinetics/calc/{identifier}/Optimization/run"
    db_sh = Path.home() / "C5O-Kinetics/db/log_calc.sh"
    # --- XTB
    cp_cmd = f"cp XTB.globalminimum.xyz {workdir}/zpv_REV/init.xyz\n"
    log_cmd = f'module load SQLite/3.45.3\nxyz_block=$(<"XTB.globalminimum.xyz")\nbash {db_sh} "{smiles}" "XTB" "OPT" 0 "${{xyz_block//\'\\n\'/}}"\n\n'
    io.orca_optimization(
        method="XTB",
        smiles=smiles,
        job_type="goat",
        guess=True,
        multiplicity=multiplicity,
        cp_cmd=cp_cmd,
        log_cmd=log_cmd,
    )
    sh_script += (
        f"cd {workdir}/goat_XTB\n"
        "xtb_id=$(sbatch submit_XTB.sh | awk '{{print $NF}}')\n\n"
    )
    # --- REVDSD
    cp_cmd = f"cp REV.xyz {workdir}/spc_CCS/init.xyz\n"
    log_cmd = (
        "module load SQLite/3.45.3\n"
        "energy=$(tac REV.log | grep 'Zero point energy' | head -n 1 | awk '{print $(NF-1)}')\n"
        'xyz_block=$(<"REV.xyz")\n'
        f'bash {db_sh} "{smiles}" "REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c" "OPT" "$energy" "${{xyz_block//\'\\n\'/}}"\n\n'
    )
    io.orca_optimization(
        method="REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c",
        smiles=smiles,
        job_type="zpv",
        multiplicity=multiplicity,
        cp_cmd=cp_cmd,
        log_cmd=log_cmd,
    )
    sh_script += (
        f"cd {workdir}/zpv_REV\n"
        f"revdsd_id=$(sbatch --dependency=afterok:$xtb_id submit_REV.sh | awk '{{print $NF}}')\n\n"
    )
    # --- CCSD(T)
    log_cmd = (
        "module load SQLite/3.45.3\n"
        "energy_hartree=$(tac CCS.log | grep 'FINAL SINGLE POINT ENERGY' | head -n 1 | awk '{print $(NF)}')\n"
        'energy=$(echo "$energy_hartree * 627.509" | bc -l)\n'
        'xyz_block=$(<"init.xyz")\n'
        f'bash {db_sh} "{smiles}" "CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/c" "OPT" "$energy" "${{xyz_block//\'\\n\'/}}"\n\n'
    )
    io.orca_optimization(
        method="CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/c",
        smiles=smiles,
        job_type="spc",
        multiplicity=multiplicity,
        log_cmd=log_cmd,
    )
    sh_script += (
        f"cd {workdir}/spc_CCS\nsbatch --dependency=afterok:$revdsd_id submit_CCS.sh"
    )

    path_out = workdir / "submit.sh"
    path_out.write_text(sh_script)

    return f"bash {path_out}"

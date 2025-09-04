from protomol.rd.mol import from_smiles, intra_proton_transfer
from rdkit import Chem
from rdkit.Chem import AllChem

from pathlib import Path

RSRCS_DCT = {
    "XTB": {
        0: [0, 0, 0, "00:00:00"],
        1: [8, 1, 10, "00:30:00"],
        2: [8, 1, 10, "00:30:00"],
        3: [8, 1, 10, "00:30:00"],
        4: [8, 1, 10, "00:30:00"],
        5: [8, 1, 10, "00:30:00"],
        6: [8, 1, 10, "00:30:00"],
        7: [8, 1, 10, "00:30:00"],
        8: [8, 1, 10, "00:30:00"],
    },  # Heavy atoms: [cpus, mem, lscratch, time]
    "REV": {
        0: [0, 0, 0, "00:00:00"],
        1: [16, 1, 20, "02:00:00"],
        2: [16, 1, 20, "02:00:00"],
        3: [16, 1, 20, "02:00:00"],
        4: [16, 1, 20, "02:00:00"],
        5: [16, 1, 20, "02:00:00"],
        6: [16, 1, 20, "02:00:00"],
        7: [16, 1, 20, "02:00:00"],
        8: [16, 1, 20, "02:00:00"],
    },
    "CCS": {
        0: [0, 0, 0, "00:00:00"],
        1: [4, 1, 50, "00:30:00"],
        2: [4, 1, 50, "00:30:00"],
        3: [12, 2, 100, "00:30:00"],
        4: [12, 4, 200, "00:60:00"],
        5: [16, 4, 200, "02:00:00"],
        6: [20, 4, 200, "02:00:00"],
        7: [20, 6, 200, "02:00:00"],
        8: [20, 12, 200, "06:00:00"],
    },
}


def _smiles_dir_identifier(smiles: str):
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


def _workdir(
    smiles,
    method_abbr,
    procedure: str = "Optimization",
):
    identifier = _smiles_dir_identifier(smiles)
    return Path.home() / f"C5O-Kinetics/calc/{identifier}/{procedure}/run/{method_abbr}"


def _files_from_template(
    smiles: str,
    method: list,
    log_command: str,
    cp_command: str,
    charge: int,
    multiplicity: int,
    procedure: str = "Optimization",
    scan_command: str = "",
    gen_guess: bool = False,
):
    identifier = _smiles_dir_identifier(smiles)
    num_heavy_ats = smiles.count("C") + smiles.count("O")
    rsrcs_lst = RSRCS_DCT[method[0][:3]][num_heavy_ats]
    replacements = {
        "[smiles]": smiles,
        "[num_cpus]": rsrcs_lst[0],
        "[size_mem]": rsrcs_lst[1] * 750,
        "[size_lscratch]": rsrcs_lst[2],
        "[time]": rsrcs_lst[3],
        "[method]": method[0],
        "[method_abbr]": str(method[0])[:3],
        "[procedure]": method[1],
        "[log_command]": log_command,
        "[cp_command]": cp_command,
        "[scan_command]": scan_command,
        "[charge]": charge,
        "[multiplicity]": multiplicity,
    }

    workdir = (
        Path.home()
        / f"C5O-Kinetics/calc/{identifier}/{procedure}/run/{replacements['[method_abbr]']}"
    )
    workdir.mkdir(parents=True, exist_ok=True)
    with open(Path(__file__).parent / "orca.inp.template", "r") as file:
        template_content = file.read()
    for placeholder, value in replacements.items():
        template_content = template_content.replace(placeholder, str(value))
    with open(workdir / f"{replacements['[method_abbr]']}.inp", "w") as file:
        file.write(template_content)
    replacements["[size_mem]"] = rsrcs_lst[1]
    with open(Path(__file__).parent / "slurm.sh.template", "r") as file:
        template_content = file.read()
    for placeholder, value in replacements.items():
        template_content = template_content.replace(placeholder, str(value))
    with open(workdir / "submit.sh", "w") as file:
        file.write(template_content)

    if gen_guess:
        _gen_guess_xyz(smiles, workdir)


def _gen_guess_xyz(smiles: str, workdir: Path):
    mol = from_smiles(smiles)
    AllChem.EmbedMolecule(mol)
    xyz = Chem.MolToXYZBlock(mol)
    (workdir / "init.xyz").write_text(xyz.lstrip())


def standard_optimization_from_smiles(
    smiles: str,
    charge: int = 0,
    multiplicity: int = None,
):
    identifier = _smiles_dir_identifier(smiles)
    if multiplicity is None and (smiles.count("[") + smiles.count("]")) == 0:
        multiplicity = 1
    else:
        multiplicity = 2

    # XTB GOAT
    method = ("XTB", "GOAT")
    log_command = (
        "module load SQLite/3.45.3\n"
        f'bash /home/tns97255/C5O-Kinetics/db/log_calc.sh "{smiles}" "XTB" "GOAT" {multiplicity} 0 "$(<"XTB.globalminimum.xyz")" ""'
    )
    cp_command = f"cp {method[0][:3]}.globalminimum.xyz /home/tns97255/C5O-Kinetics/calc/{identifier}/Optimization/run/REV/init.xyz"
    _files_from_template(
        smiles,
        method,
        log_command,
        cp_command,
        charge,
        multiplicity,
        gen_guess=True,
    )

    # REVDSD OPT NumFreq
    method = ("REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c", "OPT NumFreq")
    log_command = (
        "module load SQLite/3.45.3\n"
        "energy=$(tac REV.log | grep 'Zero point energy' | head -n 1 | awk '{print $(NF-1)}')\n"
        f'bash /home/tns97255/C5O-Kinetics/db/log_calc.sh "{smiles}" "REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c" "OPT NumFreq" {multiplicity} "$energy" "$(<"REV.xyz")" "$(<"REV.hess")"'
    )
    cp_command = f"cp {method[0][:3]}.xyz /home/tns97255/C5O-Kinetics/calc/{identifier}/Optimization/run/CCS/init.xyz"
    _files_from_template(smiles, method, log_command, cp_command, charge, multiplicity)

    # CCSDT
    method = ("CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/c", "")
    log_command = (
        "module load SQLite/3.45.3\n"
        "energy_hartree=$(tac CCS.log | grep 'FINAL SINGLE POINT ENERGY' | head -n 1 | awk '{print $(NF)}')\n"
        'energy=$(echo "$energy_hartree * 627.509" | bc -l)\n'
        f'bash /home/tns97255/C5O-Kinetics/db/log_calc.sh "{smiles}" "CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/c" "SPC" {multiplicity} "$energy" "$(<"init.xyz")" ""'
    )
    cp_command = ""
    _files_from_template(smiles, method, log_command, cp_command, charge, multiplicity)

    workdir = _workdir(smiles, "")
    with open(Path(__file__).parent / "batch.sh.template", "r") as file:
        template_content = file.read()
    template_content = template_content.replace(
        "[identifier]", _smiles_dir_identifier(smiles)
    )
    with open(workdir / "batch.sh", "w") as file:
        file.write(template_content)

    return f"bash {workdir / 'batch.sh'}"


def standard_scan_from_smiles(
    smiles: str,
    idx1: int,
    idx2: int,
    charge: int = 0,
    multiplicity: int = None,
):
    if multiplicity is None and (smiles.count("[") + smiles.count("]")) == 0:
        multiplicity = 1
    else:
        multiplicity = 2

    # XTB OPT SCAN
    method = ("XTB", "OPT")
    log_command = (
        "module load SQLite/3.45.3\n"
        "xyz_block=\n"
        f'bash /home/tns97255/C5O-Kinetics/db/log_calc.sh "{smiles}" "XTB SCAN {idx1} {idx2}" "OPT" {multiplicity} 0 "$(<"XTB.allxyz")" ""'
    )
    cp_command = ""
    xyz, scan_command = intra_proton_transfer(smiles, idx1, idx2)
    _files_from_template(
        smiles=smiles,
        method=method,
        log_command=log_command,
        cp_command=cp_command,
        charge=charge,
        multiplicity=multiplicity,
        procedure=f"Scan_{idx1}_{idx2}",
        scan_command=scan_command,
    )
    workdir = _workdir(smiles, "", f"Scan_{idx1}_{idx2}")
    (workdir / "XTB/init.xyz").write_text(xyz)
    return f"sbatch {workdir / 'XTB/submit.sh'}"


def standard_optimization_from_scan(
    smiles: str,
    xyz_block: str,
    idx1: int,
    idx2: int,
    multiplicity: int = 2,
    charge: int = 0,
):
    # REVDSD OPTTS FREQ
    method = (
        "REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c",
        "OPTTS NumFreq",
    )
    log_command = (
        "module load SQLite/3.45.3\n"
        "xyz_block=\n"
        f'bash /home/tns97255/C5O-Kinetics/db/log_calc.sh "{smiles}" "REV OPT {idx1} {idx2}" "OPTTS NumFreq" 0 "$(<"REV.xyz")" "$(<"REV.hess")"'
    )
    cp_command = ""
    scan_command = "%geom\n   Calc_Hess true\n   NumHess true\n end"
    _files_from_template(
        smiles=smiles,
        method=method,
        log_command=log_command,
        cp_command=cp_command,
        charge=charge,
        multiplicity=multiplicity,
        procedure=f"Scan_{idx1}_{idx2}",
        scan_command=scan_command,
    )
    workdir = _workdir(smiles, "", f"Scan_{idx1}_{idx2}")
    (workdir / "REV/init.xyz").write_text(xyz_block)
    return f"sbatch {workdir / 'XTB/submit.sh'}"

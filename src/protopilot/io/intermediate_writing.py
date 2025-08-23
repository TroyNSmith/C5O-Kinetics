from rdkit import Chem

from pathlib import Path


def sub_bash(
    name: str,
    methods: list,
    work_dirs: dict,
    par_dir: Path,
) -> None:
    """Generate a SLURM submission script for a sequence of ORCA jobs."""
    # Remove 'guess' if present
    work_dirs = work_dirs.copy()
    methods = methods.copy()
    work_dirs.pop("guess", None)
    if "guess" in methods:
        methods.remove("guess")

    slurm_opts = (
        "#!/bin/bash\n\n"
        "#SBATCH --partition=batch\n"
        f"#SBATCH --job-name={name}\n"
        "#SBATCH --nodes=1\n"
        "#SBATCH --ntasks=8\n"
        "#SBATCH --ntasks-per-node=8\n"
        "#SBATCH --cpus-per-task=1\n"
        "#SBATCH --time=12:00:00\n"
        "#SBATCH --mem-per-cpu=4G\n\n"
    )
    module_opts = "module load ORCA/6.1\n\n"

    keys = list(work_dirs.keys())
    job_steps = []
    for i, method in enumerate(methods):
        job = f"$(which orca) {work_dirs[keys[i]]}/{method}.inp > {work_dirs[keys[i]]}/{method}.log"
        print(i, method)
        if i < len(methods) - 1:
            # For the first method, check for *.globalminimum.xyz and copy/rename if present
            if i == 0:
                job += (
                    f"\nif compgen -G '{work_dirs[keys[i]]}/*.globalminimum.xyz' > /dev/null; then\n"
                    f"  cp {work_dirs[keys[i]]}/*.globalminimum.xyz {work_dirs[keys[i + 1]]}/{method}.xyz;\n"
                    "else\n"
                    f"  cp {work_dirs[keys[i]]}/{method}.xyz {work_dirs[keys[i + 1]]};\n"
                    "fi"
                )
            else:
                job += (
                    f"\ncp {work_dirs[keys[i]]}/{method}.xyz {work_dirs[keys[i + 1]]}"
                )
        job_steps.append(job + "\n")

    script = slurm_opts + module_opts + "".join(job_steps)
    path_out = par_dir / "submit.sh"
    path_out.write_text(script)


def guess_xyz(
    mol: Chem.rdchem.Mol,
    output_dir: Path,
    xyz_name: str,
) -> None:
    """Write an RDKit molecule to an XYZ file in the specified directory."""
    path_out = output_dir / f"{xyz_name}.xyz"
    xyz = Chem.MolToXYZBlock(mol)
    path_out.write_text(xyz.lstrip())


def write_orca_input(
    method: str,
    charge: int,
    multiplicity: int,
    output_dir: Path,
    xyz_name: str,
    job_type: str = "opt",
) -> None:
    """
    Write an ORCA input file for various calculation types.

    job_type options:
        - "opt": Optimization only (no frequency calculation) (default)
        - "opt_freq": Optimization and frequency calculation
        - "goat": Geometry optimization and frequency calculation (GOAT)
        - "elec": Single-point electronic structure calculation
    """
    pal = "! PAL8\n"
    if job_type == "opt_freq":
        job_line = f"! {method}\n! OPT NumFreq\n* xyzfile {charge} {multiplicity} {output_dir}/{xyz_name}.xyz\n"
    elif job_type == "opt":
        job_line = f"! {method}\n! OPT\n* xyzfile {charge} {multiplicity} {output_dir}/{xyz_name}.xyz\n"
    elif job_type == "goat":
        job_line = f"! {method}\n! GOAT\n* xyzfile {charge} {multiplicity} {output_dir}/{xyz_name}.xyz\n"
    elif job_type == "elec":
        job_line = f"! {method}\n* xyzfile {charge} {multiplicity} {output_dir}/{xyz_name}.xyz\n"
    else:
        raise ValueError(f"Unknown job_type: {job_type}")
    txt = f"{pal}{job_line}"
    inp_name = (
        method.replace(" ", "_")
        .replace("/", "")
        .replace(")", "")
        .replace("(", "")
        .replace("-", "_")
    )
    path_out = output_dir / f"{inp_name}.inp"
    path_out.write_text(txt)

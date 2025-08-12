from rdkit import Chem

from pathlib import Path


class Optimization:
    def xyz(mol: Chem.rdchem.Mol, output_dir: Path):
        path_out = output_dir / Path("guess.xyz")
        xyz = Chem.MolToXYZBlock(mol)
        path_out.write_text(" " + xyz)

    def inp(method: str, charge: int, multiplicity: int, output_dir: Path, n_cpus: int):
        txt = f"! PAL{n_cpus}\n! {method}\n! OPT\n* xyzfile {charge} {multiplicity} guess.xyz\n"
        path_out = output_dir / Path("opt.inp")
        path_out.write_text(txt)

    def bash(name: str, output_dir: Path, n_cpus: int):
        txt = (
            "#!/bin/bash\n\n#SBATCH"
            "--partition=batch\n"
            f"#SBATCH --job-name=opt_{name}\n"
            f"#SBATCH --ntasks={n_cpus}\n"
            "#SBATCH --time=4:00:00\n"
            "#SBATCH --mem-per-cpu=10G\n\n"
            "module load ORCA/6.1\n"
            "$(which orca) opt.inp > opt.log"
        )
        path_out = output_dir / Path("run_opt.sh")
        path_out.write_text(txt)

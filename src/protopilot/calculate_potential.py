import click
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from pathlib import Path

from .util import writing


@click.group()
def cli():
    """Main CLI entry point."""
    pass


@cli.command
@click.option(
    "-smi",
    "--smiles",
    "smiles_string",
    help="Structure denoted by SMILES string format.",
    type=str,
    default=None,
)
@click.option(
    "-m",
    "--method",
    "method",
    help="Method for performing quantum calculations.",
    type=str,
    default="wB97X-3c",
)
@click.option(
    "-ncpu",
    "--number_cpus",
    "n_cpus",
    help="Number of cpus to utilize for parallel processing.",
    type=int,
    default=8,
)
def optimize(smiles_string: str, method: str, n_cpus: int):
    name = smiles_string.replace("[", "_").replace("]", "_").replace("=", "dbl")
    output_dir = Path(f"calc/{name}/run_opt/")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if not Path.exists(output_dir / Path("guess.xyz")) or True is True:
        mol_format = Chem.AddHs(Chem.MolFromSmiles(smiles_string))
        AllChem.EmbedMolecule(mol_format)

        charge = Chem.GetFormalCharge(mol_format)
        multiplicity = Descriptors.NumRadicalElectrons(mol_format) + 1
        writing.Optimization.xyz(mol_format, output_dir)
        writing.Optimization.inp(method, charge, multiplicity, output_dir, n_cpus)
        writing.Optimization.bash(name, output_dir, n_cpus)

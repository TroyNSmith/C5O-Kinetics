"""
Main CLI entry point for ORCA workflows.
Dispatches between intermediate and transition state (TS) workflows.
"""

import click


@click.group()
def cli():
    """Main CLI for protopilot workflows."""
    pass


@cli.command()
@click.option(
    "-smi",
    "--smiles",
    "smiles_string",
    required=True,
    help="Structure denoted by SMILES string format.",
    type=str,
)
def intermediate(smiles_string):
    """Run the intermediate structure workflow."""
    from protopilot.workflow import intermediate_bu

    intermediate_bu.main(smiles_string)


@cli.command()
def ts():
    """Run the transition state workflow (to be implemented)."""
    from protopilot.workflow import transition

    transition.main()


if __name__ == "__main__":
    cli()

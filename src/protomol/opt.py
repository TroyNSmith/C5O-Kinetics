import click
from orca.writing import write_orca
import pandas as pd
from rd.mol import smiles_to_xyz_block
from util.sql import execute_append, refresh

from pathlib import Path
import sqlite3


@click.group()
def cli():
    refresh()
    pass


def display_methods_table(method_dict: dict[int, str], title: str):
    print(f"\n{title}")
    print(f"  {'Method ID':<10} | {'Method':<15} | Functional/Basis")
    print(f"  {'=' * 10:<10} | {'=' * 15:<15} | {'=' * 20}")
    for _, m in method_dict.items():
        print(m)


@cli.command()
@click.option("-smi", "--smiles", "smiles", help="Input corresponding SMILES string")
def optimize(smiles: str):
    """Generate an ORCA command for a given SMILES string."""
    db = Path.home() / "C5O-Kinetics/db/data.db"
    conn = sqlite3.connect(db)

    # Get all SMILES and methods
    df_smiles = pd.read_sql_query("SELECT smiles_id, smiles_text FROM smiles", conn)
    smiles_dict = {row.smiles_id: row.smiles_text for row in df_smiles.itertuples()}

    df_methods = pd.read_sql_query("SELECT * FROM methods", conn)
    methods_dict = {
        row.method_id: f"  {row.method_id:<10} | {row.method:<15} | {row.functional} {row.basis}"
        for row in df_methods.itertuples()
    }

    if smiles in smiles_dict.values():
        # SMILES already exists — retrieve existing calculations
        smiles_id = next(k for k, v in smiles_dict.items() if v == smiles)
        df_calculations = pd.read_sql_query(
            f"SELECT * FROM calculations WHERE smiles_id = {smiles_id}", conn
        )

        if df_calculations.empty:
            print("No previous calculations found for this SMILES.")
            return

        # Build calc_id → method summary mapping
        existing_methods = {
            row.calc_id: methods_dict.get(row.method_id, "Unknown method")
            for row in df_calculations.itertuples()
        }

        # Display existing methods for this SMILES
        print("\n  Existing Calculation Methods")
        print("  " + "=" * 55)
        for val in existing_methods.values():
            print(val)

        # Prompt user to select method for initial geometry
        init_method_id = int(input("\nSelect a method id for initial run: "))
        init_calc_id = next(
            (
                row.calc_id
                for row in df_calculations.itertuples()
                if row.method_id == init_method_id
            ),
            None,
        )
        assert init_calc_id, (
            f"Method ID {init_method_id} not found in existing calculations."
        )
        df_xyzs = pd.read_sql_query(
            "SELECT * FROM xyz WHERE calc_id = ?", conn, params=(init_calc_id,)
        )
        if df_xyzs.shape[0] != 1:
            raise ValueError(
                f"Expected 1 match in SMILES table, found {df_smiles.shape[0]} for {smiles}"
            )
        initial_xyz = df_xyzs.iloc[0]["xyz_text"]

        # Show all methods
        print("\n  Available Methods")
        print("  " + "=" * 55)
        for val in methods_dict.values():
            print(val)

    else:
        # New SMILES — prompt user to add
        multiplicity = int(input("\nIndicate multiplicity of SMILES molecule: "))
        initial_xyz = smiles_to_xyz_block(smiles)
        execute_append(
            "INSERT INTO smiles (smiles_text, multiplicity, initial) VALUES (?, ?, ?)",
            (smiles, multiplicity, initial_xyz),
            db,
        )

        # Show only "GOAT" methods (substring check)
        goat_methods = {
            row.method_id: methods_dict[row.method_id]
            for row in df_methods.itertuples()
            if "GOAT" in {row.method}
        }

        if not goat_methods:
            print("No GOAT methods available.")
            return

        print("\n  Available GOAT Methods")
        print("  " + "=" * 55)
        for val in goat_methods.values():
            print(val)

    new_method_id = int(input("\nSelect a method id for new run: "))
    shell_cmd = write_orca(smiles, new_method_id, initial_xyz)
    print("\nCopy and paste the following command into your terminal:\n")
    print(shell_cmd, "\n")


if __name__ == "__main__":
    cli()

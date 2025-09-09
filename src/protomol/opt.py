import click
from orca.writing import write_orca
import pandas as pd
from rd.mol import beta_cleavage_indices, proton_transfer_indices, smiles_to_xyz_block
from util.ref import METHOD_MAP
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

    idx1 = -1
    idx2 = -1

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
        existing_methods = {}
        seen_method_ids = set()

        for row in df_calculations.itertuples():
            if row.method_id not in seen_method_ids:
                existing_methods[row.method_id] = methods_dict.get(
                    row.method_id, "Unknown method"
                )
                seen_method_ids.add(row.method_id)

        # Display existing methods for this SMILES
        print("\n  Existing Calculation Methods")
        print("  " + "=" * 55)
        for val in existing_methods.values():
            print(val)

        init_method_id = int(input("\nSelect a method id for initial run: "))
        matching_calcs = [
            row
            for row in df_calculations.itertuples()
            if row.method_id == init_method_id
        ]

        if not matching_calcs:
            raise ValueError(f"No calculations found for method ID {init_method_id}.")

        if len(matching_calcs) == 1:
            selected_calc = matching_calcs[0]
        else:
            print(f"\nMultiple calculations found for method ID {init_method_id}.")
            print(
                "Select which calculation you'd like to use as the starting geometry:"
            )
            for i, row in enumerate(matching_calcs):
                scan_str = (
                    f"Scan {row.scan_idx1}:{row.scan_idx2}"
                    if row.scan_idx1 != -1 and row.scan_idx2 != -1
                    else "Original Geometry"
                )
                print(f"  {i + 1}. {scan_str} (calc_id = {row.calc_id})")

            selected_idx = int(input("\nEnter the number of the desired calculation: "))
            if not (1 <= selected_idx <= len(matching_calcs)):
                raise ValueError("Invalid selection.")

            selected_calc = matching_calcs[selected_idx - 1]

        init_calc_id = selected_calc.calc_id
        idx1 = selected_calc.scan_idx1
        idx2 = selected_calc.scan_idx2

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

        # Get method name for the selected init method
        init_method_row = df_methods[df_methods.method_id == init_method_id].iloc[0]
        init_method_name = init_method_row["method"]

        # Determine allowed next methods
        allowed_methods = METHOD_MAP.get(init_method_name, [])

        # Filter methods to only allowed ones
        filtered_methods = {
            row.method_id: methods_dict[row.method_id]
            for row in df_methods.itertuples()
            if row.method in allowed_methods
        }

        if not filtered_methods:
            print(
                f"\nNo valid next-step methods found for initial method '{init_method_name}'."
            )
            return

        # Show only allowed new methods
        print(f"\n  Available Methods after '{init_method_name}'")
        print("  " + "=" * 55)
        for val in filtered_methods.values():
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

        # Get method name for GOAT methods
        df_goat = df_methods[df_methods.method == "GOAT"]

        filtered_methods = {
            row.method_id: methods_dict[row.method_id] for row in df_goat.itertuples()
        }

        # Show only allowed new methods
        print("\n  Available GOAT Methods'")
        print("  " + "=" * 55)
        for val in filtered_methods.values():
            print(val)
    allowed_ids = set(filtered_methods.keys())
    new_method_id = int(input("\nSelect a method id for new run: "))
    if new_method_id not in allowed_ids:
        raise ValueError(
            f"Method ID {new_method_id} is not valid after '{init_method_name}'."
        )

    method = methods_dict[new_method_id]
    mechanism = None
    if "SCAN" in method:
        print("\n  Available Mechanisms")
        print("  " + "=" * 55)
        print("  1 | Beta Cleavage\n  2 | Proton Transfer")
        mechanism_inp = int(input("\nSelect a mechanism: "))

        if mechanism_inp == 1:
            mechanism = "beta cleavage"
            alpha_list, beta_list = beta_cleavage_indices(smiles)

            idx1 = int(
                input(f"\nSelect an alpha index for beta cleavage {alpha_list}: ")
            )
            if idx1 not in alpha_list:
                raise ValueError(
                    f"Index {idx1} is not a valid alpha index for beta cleavage."
                )

            idx2 = int(input(f"\nSelect a beta index for beta cleavage {beta_list}: "))
            if idx2 not in beta_list:
                raise ValueError(
                    f"Index {idx2} is not a valid beta index for beta cleavage."
                )

        elif mechanism_inp == 2:
            mechanism = "proton transfer"
            radical_list, proton_list = proton_transfer_indices(smiles)

            idx1 = int(input(f"\nSelect a proton index {proton_list}: "))
            if idx1 not in proton_list:
                raise ValueError(f"Index {idx1} is not a valid proton index.")

            idx2 = int(input(f"\nSelect a radical index {radical_list}: "))
            if idx2 not in radical_list:
                raise ValueError(f"Index {idx2} is not a valid radical index.")

        else:
            raise ValueError("Select the integer corresponding to a valid mechanism.")

    shell_cmd = write_orca(smiles, new_method_id, initial_xyz, idx1, idx2, mechanism)
    print("\nCopy and paste the following command into your terminal:\n")
    print(shell_cmd, "\n")


if __name__ == "__main__":
    cli()

from protomol.util.backend import Query_SQL, write_orca, Append_SQL
from protomol.util.update_db import refresh
import click


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
    smiles_list = Query_SQL.smiles()
    methods = sorted(set(Query_SQL.methods()))



    if smiles in smiles_list:
        # SMILES already exists — retrieve existing calculations
        calculations = Query_SQL.calculations(smiles)

        # Build mapping of calc_id → method
        existing_methods = {
            calc[0]: m
            for calc in calculations
            for m in methods
            if str(calc[2]) in m[:5]
        }

        display_methods_table(
            {int(m[:5]): m for m in existing_methods.values()},
            "Existing Calculation Methods",
        )

        # Prompt for initial method ID (matches a calc_id)
        init_method_id = int(input("\nSelect a method id for initial run: "))
        init_calc_id = next(
            (
                calc_id
                for calc_id, m in existing_methods.items()
                if str(init_method_id) in m[:5]
            ),
            None,
        )

        if init_calc_id is None:
            raise ValueError(
                f"Method ID {init_method_id} not found in existing calculations."
            )

        # Show all methods for new calculation
        all_method_dict = {int(m[:5]): m for m in methods}
        display_methods_table(all_method_dict, "Available Methods")

        # Prompt for new method selection
        new_method_id = int(input("\nSelect a method id for new run: "))

        # Generate initial files and shell command
        shell_cmd = write_orca(smiles, new_method_id, init_calc_id)
        print("\nCopy and paste the following command into your terminal:\n")
        print(shell_cmd, "\n")

    else:
        # SMILES is new — prompt for multiplicity and method
        multiplicity = int(input("\nIndicate multiplicity of SMILES molecule: "))
        Append_SQL.smiles(smiles, multiplicity)

        # Filter methods that contain "GOAT" between positions 14–31
        goat_methods = [m for m in methods if "GOAT" in m[14:31]]
        goat_method_dict = {int(m[:5]): m for m in goat_methods}

        display_methods_table(goat_method_dict, "Available GOAT Methods")

        new_method_id = int(input("\nSelect a method id for new run: "))

        shell_cmd = write_orca(smiles, new_method_id, multiplicity)
        print("\nCopy and paste the following command into your terminal:\n")
        print(shell_cmd, "\n")


@cli.command()
@click.option("-smi", "--smiles", "smiles", help="input corresponding SMILES string")
def visualize(smiles: str):
    smi_lst = Query_SQL.smiles()
    mthds = sorted(set(Query_SQL.methods()))
    if smiles in smi_lst:
        print(f"\n  {'Method ID':<10} | {'Method':<15} | Functional/Basis")
        print(f"  {'=' * 10:<10} | {'=' * 15:<15} | {'=' * 20}")
        calcs = Query_SQL.calculations(smiles)
        init_mthds = {
            calc[0]: m for calc in calcs for m in mthds if str(calc[2]) in m[:5]
        }
        for calc_id, m in init_mthds.items():
            print(f"{m}")
        init_mthd_id = int(input("\nSelect a method id for initial run: "))
        init_calc_id = next(
            (
                calc_id
                for calc_id, m in init_mthds.items()
                if str(init_mthd_id) in m[:5]
            ),
            None,
        )
        xyz = ""


if __name__ == "__main__":
    cli()

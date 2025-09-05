from ipywidgets import Layout
import matplotlib.pyplot as plt
from rdkit import Chem

from pathlib import Path
import re
import sqlite3

from ..rd.mol import from_smiles
from .render import from_allxyz_block, from_xyz_block

calc_db = Path.home() / "C5O-Kinetics/db/results.db"
freq_db = Path.home() / "C5O-Kinetics/db/frequencies.db"

new_db = Path.home() / "C5O-Kinetics/db/data.db"


class Query_SQL:
    def _execute_query(
        execute: str,
        db: str | Path,
        target_value: str = None,
    ):
        conn = sqlite3.connect(db)
        cursor = conn.cursor()
        if target_value:
            cursor.execute(execute, (target_value,))
        else:
            cursor.execute(execute)
        rows = cursor.fetchall()
        cursor.close()
        return rows

    def smiles() -> list:
        execute = "SELECT DISTINCT smiles_text FROM smiles"
        rows = Query_SQL._execute_query(execute, new_db)
        return [rows[i][0] for i in range(len(rows))]

    def methods() -> list:
        execute = "SELECT DISTINCT method, method_id FROM methods"
        rows = Query_SQL._execute_query(execute, new_db)
        return [(rows[i][0], rows[i][1]) for i in range(len(rows))]

    def functionals(method_id: int) -> list:
        execute = "SELECT DISTINCT functional, basis FROM methods WHERE method_id = ?"
        rows = Query_SQL._execute_query(execute, new_db, method_id)
        return [f"{rows[i][0]} {rows[i][1]}" for i in range(len(rows))]

    def xyzs(smiles: str, method_id: int) -> str:
        execute = "SELECT DISTINCT smiles_id FROM smiles WHERE smiles_text = ?"
        smiles_id = Query_SQL._execute_query(execute, new_db, smiles)[0][0]
        execute = f"SELECT DISTINCT xyz_text FROM XYZS WHERE smiles_id = {smiles_id} AND method_id = {method_id}"
        rows = Query_SQL._execute_query(execute, new_db)
        if len(rows) > 1:
            return [f"{rows[i][0]} {rows[i][1]}" for i in range(len(rows))]
        else:
            return []

    def calculations(smiles: str, method_id: int) -> list:
        return


class Append_SQL:
    def _execute_append(
        execute: str,
        val: list[tuple],
        db: str | Path,
    ):
        conn = sqlite3.connect(db)
        cursor = conn.cursor()
        cursor.execute(execute, val)
        conn.commit()
        conn.close()

    def smiles(smiles: str, mult: int):
        initial = Chem.MolToXYZBlock(from_smiles(smiles, with_coords=True))
        execute = (
            "INSERT INTO smiles (smiles_text, multiplicity, initial) VALUES (?, ?, ?)"
        )
        Append_SQL._execute_append(execute, (smiles, mult, initial), new_db)

    def _method(
        functional: str,
        basis: str,
        method: str,
        inp_template: str,
        submit_template: str,
    ):
        execute = "INSERT INTO methods (functional, basis, method, inp_template, submit_template) VALUES (?, ?, ?, ?, ?)"
        Append_SQL._execute_append(
            execute, (functional, basis, method, inp_template, submit_template), new_db
        )


class Styles:
    def scale_width(width: str = "50%"):
        return Layout(height="30", width=width)


def on_visualize_xyz(xyz_block: str):
    if not xyz_block == "N/A":
        if "Relaxed Surface Scan" in xyz_block.value:
            from_allxyz_block(xyz_block.replace(">\n", ""))
            matches = re.findall(r"Step (\d+) E ([-\d\.]+)", xyz_block[1, :2].value)
            if matches:
                steps = [int(step) for step, energy in matches]
                energies = [float(energy) for step, energy in matches]
                fig, ax = plt.subplots()
                ax.plot(steps, energies, marker="o")
                ax.set_xlabel("Step")
                ax.set_ylabel("Energy")
                plt.show()
        else:
            from_xyz_block(xyz_block, add_labels=True)

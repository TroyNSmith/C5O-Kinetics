from ipywidgets import Layout
import matplotlib.pyplot as plt
from rdkit import Chem

from pathlib import Path
import re
import sqlite3

from ..rd.mol import beta_cleavage, from_smiles, intra_proton_transfer
from ..util.ref import get_ccsdt_parameters
from .render import from_allxyz_block, from_xyz_block

calc_db = Path.home() / "C5O-Kinetics/db/results.db"
freq_db = Path.home() / "C5O-Kinetics/db/frequencies.db"

new_db = Path.home() / "C5O-Kinetics/db/data.db"


class Query_SQL:
    def _execute_query(
        execute: str,
        db: str | Path,
        target_value: str = None,
        many: bool = False,
    ):
        conn = sqlite3.connect(db)
        cursor = conn.cursor()

        if target_value:
            if many:  # Can only detect first match for each tuple
                rows = []
                for tup in target_value:
                    cursor.execute(execute, tup)
                    row = cursor.fetchall()
                    rows.append(row[0])
                return rows
            else:
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

    def smiles_id(smiles: str) -> list:
        execute = "SELECT smiles_id FROM smiles WHERE smiles_text = ?"
        rows = Query_SQL._execute_query(execute, new_db, smiles)
        if len(rows) > 1:
            raise ValueError(f"Multiple matching SMILES strings found in {new_db}")
        elif len(rows) < 1:
            raise ValueError(f"No matching SMILES strings found in {new_db}")
        else:
            return rows[0][0]

    def methods() -> list:
        execute = "SELECT DISTINCT method, functional, basis, method_id FROM methods"
        rows = Query_SQL._execute_query(execute, new_db)
        return [
            (f"{rows[i][0].strip():<20} | {rows[i][1]} {rows[i][2]}", rows[i][3])
            for i in range(len(rows))
        ]

    def calculations(
        smiles: str, method_id: int = None, idx1: int = None, idx2: int = None
    ):
        smiles_id = Query_SQL.smiles_id(smiles)
        if method_id:
            if idx1 or idx2:
                assert idx1 is not None and idx2 is not None, (
                    "Must provide two indices when doing an index search."
                )
                execute = f"SELECT * FROM calculations WHERE smiles_id = {smiles_id} AND method_id = {method_id} AND scan_idx1 = {idx1} AND scan_idx2 = {idx2}"
            else:
                execute = f"SELECT * FROM calculations WHERE smiles_id = {smiles_id} AND method_id = {method_id}"
        else:
            execute = f"SELECT * FROM calculations WHERE smiles_id = {smiles_id}"
        return Query_SQL._execute_query(execute, new_db)

    def xyzs(calc_id: int = None) -> str:
        execute = "SELECT xyz_text FROM xyz WHERE calc_id = ?"
        xyzs = Query_SQL._execute_query(execute, new_db, calc_id)
        if len(xyzs) != 1:
            raise ValueError(
                f"Error fetching xyz_text. Check for either multiple or no entries in {new_db}."
            )
        return xyzs[0][0]

    def traj(calc_id: int = None) -> str:
        execute = "SELECT traj_text FROM traj WHERE calc_id = ?"
        traj = Query_SQL._execute_query(execute, new_db, calc_id)
        if len(traj) != 1:
            raise ValueError(
                f"Error fetching traj_text. Check for either multiple or no entries in {new_db}."
            )
        return traj[0][0]


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

    def calculation(
        smiles_id: int,
        method_id: int,
        scan_idx1: int = -1,
        scan_idx2: int = -1,
        return_id: bool = True,
    ):
        execute = "INSERT INTO calculations (smiles_id, method_id, scan_idx1, scan_idx2) VALUES (?, ?, ?, ?)"
        Append_SQL._execute_append(
            execute, (smiles_id, method_id, scan_idx1, scan_idx2), new_db
        )
        if return_id:
            execute = f"SELECT calc_id FROM calculations WHERE smiles_id = {smiles_id} AND method_id = {method_id}"
            return Query_SQL._execute_query(execute, new_db)[0][0]

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


def write_orca(
    smiles: str,
    method_id: int,
    init_calc_id: int = 0,
    idx1: int = 0,
    idx2: int = 0,
    mechanism: str = None,
):
    execute = "SELECT smiles_id, multiplicity FROM smiles WHERE smiles_text = ?"
    smiles_id, multiplicity = Query_SQL._execute_query(execute, new_db, smiles)[0][:]
    execute = "SELECT functional, method, inp_template, submit_template FROM methods WHERE method_id = ?"
    rows = Query_SQL._execute_query(execute, new_db, method_id)
    functional, method, inp_template, submit_template = rows[0][:]
    rows = Query_SQL.calculations(smiles, method_id, idx1, idx2)
    # if len(rows) > 0:
    #     raise FileExistsError("Calculation matching indicated values already exists.")

    if idx1 == 0 and idx2 == 0:
        calc_id = Append_SQL.calculation(smiles_id, method_id, return_id=True)
        outdir = Path.home() / f"C5O-Kinetics/calc/{calc_id}"
        outdir.mkdir(exist_ok=True)
        if method == "GOAT":
            execute = "SELECT initial FROM smiles WHERE smiles_text = ?"
            initial_xyz = Query_SQL._execute_query(execute, new_db, smiles)[0][0]
            (outdir / "init.xyz").write_text(initial_xyz)
        else:
            assert init_calc_id > 0, (
                "Must provide an initial calculation id tying to the initial xyz structure for this calculation."
            )
            initial_xyz = Query_SQL.xyzs(init_calc_id)
            (outdir / "init.xyz").write_text(initial_xyz)
        substitutions = {
            "[SMILES]": re.sub(r"[^a-zA-Z0-9\s]", "", smiles),
            "[multiplicity]": multiplicity,
        }
        if functional == "CCSD(T)-F12/RI":
            substitutions = substitutions | get_ccsdt_parameters(
                smiles.count("C") + smiles.count("O")
            )
        for key, val in substitutions.items():
            inp_template = inp_template.replace(key, str(val))
            submit_template = submit_template.replace(key, str(val))
        (outdir / "calc.inp").write_text(inp_template)
        (outdir / "submit.sh").write_text(submit_template)
    else:
        assert idx1 + idx2 > 0, "Please provide indices for the path."
        calc_id = Append_SQL.calculation(
            smiles_id, method_id, idx1, idx2, return_id=True
        )
        outdir = Path.home() / f"C5O-Kinetics/calc/{calc_id}"
        outdir.mkdir(exist_ok=True)
        substitutions = {
            "[SMILES]": re.sub(r"[^a-zA-Z0-9\s]", "", smiles),
            "[multiplicity]": multiplicity,
            "[idx1]": idx1,
            "[idx2]": idx2,
        }
        if mechanism.lower() == "proton transfer":
            initial_xyz, min_dist, max_dist = intra_proton_transfer(smiles, idx1, idx2)
            substitutions["[start]"] = min_dist
            substitutions["[end]"] = max_dist
            raise ValueError(substitutions)
        elif mechanism.lower() == "beta cleavage":
            initial_xyz, min_dist, max_dist = beta_cleavage(smiles, idx1, idx2)
            substitutions["[start]"] = min_dist
            substitutions["[end]"] = max_dist
        (outdir / "init.xyz").write_text(initial_xyz)
        for key, val in substitutions.items():
            inp_template = inp_template.replace(key, str(val))
            submit_template = submit_template.replace(key, str(val))
        (outdir / "calc.inp").write_text(inp_template)
        (outdir / "submit.sh").write_text(submit_template)

    return f"cd {outdir}; sbatch submit.sh"


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

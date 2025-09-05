from pathlib import Path
import sqlite3

new_db = Path.home() / "C5O-Kinetics/db/data.db"
calc_dir = Path.home() / "C5O-Kinetics/calc"


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


# Map existing entries from result tables


def refresh():
    execute = "SELECT calc_id FROM xyz"
    rows = _execute_query(execute, new_db)
    xyz_map = [rows[i][0] for i in range(len(rows))]

    execute = "SELECT calc_id FROM traj"
    rows = _execute_query(execute, new_db)
    traj_map = [rows[i][0] for i in range(len(rows))]

    execute = "SELECT calc_id FROM energies"
    rows = _execute_query(execute, new_db)
    energy_map = [rows[i][0] for i in range(len(rows))]

    execute = "SELECT calc_id FROM calculations"
    rows = _execute_query(execute, new_db)
    calc_ids = [rows[i][0] for i in range(len(rows))]
    for calc_id in calc_ids:
        workdir = calc_dir / f"{calc_id}"
        if calc_id not in xyz_map:
            try:
                result = [f for f in workdir.rglob("calc.xyz")][0].read_text()
                execute = "INSERT INTO xyz (calc_id, xyz_text) VALUES (?, ?)"
                _execute_append(execute, (calc_id, result), new_db)
            except Exception:
                pass
        if calc_id not in traj_map:
            try:
                result = [f for f in workdir.rglob("calc_trj.xyz")][0].read_text()
                execute = "INSERT INTO traj (calc_id, traj_text) VALUES (?, ?)"
                _execute_append(execute, (calc_id, result), new_db)
            except Exception:
                pass
        if calc_id not in energy_map:
            try:
                log = [f for f in workdir.rglob("calc.log")][0].read_text()
                results = [
                    line for line in log.splitlines() if "Zero point energy" in line
                ]
                if len(results) == 1:
                    # Take the kcal/mol from string (second to last block)
                    result = float(results[0].split()[-2])
                else:
                    results = [
                        line
                        for line in log.splitlines()
                        if "FINAL SINGLE POINT ENERGY" in line
                    ]
                    if len(results) == 1:
                        # Take the hartrees from string (last block) & convert to kcal/mol
                        result = float(results[0].split()[-1]) * 627.509
                    else:
                        result = 0
                execute = "INSERT INTO energies (calc_id, energy_value) VALUES (?, ?)"
                _execute_append(execute, (calc_id, result), new_db)
            except Exception:
                pass

from pathlib import Path
import sqlite3
from typing import Optional, Union

# === Global Paths ===
db = Path.home() / "C5O-Kinetics/db/data.db"
calc_dir = Path.home() / "C5O-Kinetics/calc"


# === SQL Execution Helpers ===
def execute_query(
    execute: str,
    db: Union[str, Path],
    target_value: Optional[str] = None,
):
    """Execute a SELECT query and return all rows."""
    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()
        if target_value:
            cursor.execute(execute, (target_value,))
        else:
            cursor.execute(execute)
        return cursor.fetchall()


def execute_append(
    execute: str,
    val: tuple,
    db: Union[str, Path],
    return_id: bool = False,
):
    """Execute an INSERT/UPDATE/REPLACE command with values."""
    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()
        cursor.execute(execute, val)
        conn.commit()
        if return_id:
            return cursor.lastrowid


def get_existing_calc_ids(table: str) -> set[int]:
    rows = execute_query(f"SELECT calc_id FROM {table}", db)
    return {row[0] for row in rows}


def read_file_or_none(path: Path) -> Optional[str]:
    try:
        return path.read_text()
    except Exception:
        return None


def extract_energy_from_log(log_text: str) -> float:
    """Try to extract energy from log file text."""
    lines = log_text.splitlines()

    zpe_lines = [line for line in lines if "Zero point energy" in line]
    if len(zpe_lines) == 1:
        return float(zpe_lines[0].split()[-2])

    sp_lines = [line for line in lines if "FINAL SINGLE POINT ENERGY" in line]
    if len(sp_lines) == 1:
        return float(sp_lines[0].split()[-1]) * 627.509  # Convert Hartree to kcal/mol

    return 0.0


def refresh():
    xyz_ids = get_existing_calc_ids("xyz")
    traj_ids = get_existing_calc_ids("traj")
    energy_ids = get_existing_calc_ids("energies")
    calc_ids = get_existing_calc_ids("calculations")

    # xyz_ids = []
    # traj_ids = []
    # energy_ids = []

    for calc_id in calc_ids:
        workdir = calc_dir / str(calc_id)

        # === XYZ Handling ===
        if calc_id not in xyz_ids:
            xyz_file = next(workdir.rglob("calc.allxyz"), None)
            if xyz_file is None:
                xyz_file = next(workdir.rglob("calc.xyz"), None)
            if xyz_file:
                xyz_text = read_file_or_none(xyz_file)
                if xyz_text:
                    if "Scan Step" in xyz_text:
                        try:
                            lines = xyz_text.splitlines()

                            # Extract steps and energies
                            steps = [
                                int(line.split()[-3])
                                for line in lines
                                if "Scan Step" in line
                            ]
                            energies = [
                                float(line.split()[-1])
                                for line in lines
                                if "Scan Step" in line
                            ]

                            # Find the step with maximum energy
                            max_index = energies.index(max(energies))
                            max_step = steps[max_index]

                            # Determine number of digits from filenames, fallback to 3 digits
                            step_files = list(workdir.glob("calc.*.xyz"))
                            if step_files:
                                num_digits = max(
                                    len(f.stem.split(".")[-1])
                                    for f in step_files
                                    if f.stem.split(".")[-1].isdigit()
                                )
                            else:
                                num_digits = 3  # fallback

                            # Format the step number with leading zeros
                            step_str = str(max_step).zfill(num_digits)
                            step_xyz_file = workdir / f"calc.{step_str}.xyz"

                            # Read and use that specific step file if it exists
                            step_xyz_text = read_file_or_none(step_xyz_file)
                            if step_xyz_text:
                                xyz_text = step_xyz_text  # override with highest-energy structure
                        except Exception as e:
                            print(f"Failed to process scan: {e}")

                    execute_append(
                        "INSERT OR REPLACE INTO xyz (calc_id, xyz_text) VALUES (?, ?)",
                        (calc_id, xyz_text),
                        db,
                    )
        # === Energy from log file ===
        if calc_id not in energy_ids:
            log_file = next(workdir.rglob("calc.log"), None)
            if log_file:
                log_text = read_file_or_none(log_file)
                if log_text:
                    energy = extract_energy_from_log(log_text)
                    execute_append(
                        "INSERT OR REPLACE INTO energies (calc_id, energy_value) VALUES (?, ?)",
                        (calc_id, energy),
                        db,
                    )

        # === Trajectory and Scan Energy Handling ===
        if calc_id not in traj_ids:
            traj_file = next(workdir.rglob("calc.allxyz"), None)
            if traj_file is None:
                traj_file = next(workdir.rglob("calc_trj.xyz"), None)
            traj_text = read_file_or_none(traj_file) if traj_file else None

            if traj_text:
                lines = traj_text.splitlines()
                if "Scan Step" in traj_text:
                    try:
                        steps = [
                            int(line.split()[-3])
                            for line in lines
                            if "Scan Step" in line
                        ]
                        energies = [
                            float(line.split()[-1])
                            for line in lines
                            if "Scan Step" in line
                        ]
                        execute_append(
                            "INSERT OR REPLACE INTO energies (calc_id, energy_value, energy_array) VALUES (?, ?, ?)",
                            (calc_id, 0, f"[{steps}, {energies}]"),
                            db,
                        )
                    except Exception:
                        pass

                execute_append(
                    "INSERT OR REPLACE INTO traj (calc_id, traj_text) VALUES (?, ?)",
                    (calc_id, traj_text),
                    db,
                )

import numpy as np

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
    if len(zpe_lines) in (1, 2):
        return float(zpe_lines[-1].split()[-2])

    spc_lines = [line for line in lines if "FINAL SINGLE POINT ENERGY" in line]
    if len(spc_lines) == 1:
        return float(spc_lines[0].split()[-1]) * 627.509  # Convert Hartree to kcal/mol

    return 0.0


def refresh():
    xyz_ids = get_existing_calc_ids("xyz")
    traj_ids = get_existing_calc_ids("traj")
    imaginary_ids = get_existing_calc_ids("imaginaryfrequencies")

    # xyz_ids = []
    # traj_ids = []
    # energy_ids = []

    rows = execute_query("SELECT * FROM calculations", db)
    for row in rows:
        calc_id, smiles_id, method_id, idx1, idx2 = row
        workdir = calc_dir / str(calc_id)

        # === XYZ Handling ===
        if calc_id not in xyz_ids:
            xyz_file = next(workdir.rglob("calc.allxyz"), None)
            if xyz_file is None:
                xyz_file = next(workdir.rglob("calc.xyz"), None)
                selected_step = 0
            if xyz_file:
                xyz_text = read_file_or_none(xyz_file)
                if xyz_text is not None and "Scan Step" in xyz_text:
                    lines = xyz_text.splitlines()
                    # Extract steps and energies
                    steps = [
                        int(line.split()[-3]) for line in lines if "Scan Step" in line
                    ]
                    energies = [
                        float(line.split()[-1]) for line in lines if "Scan Step" in line
                    ]

                    def identify_ts_index(energies: list[float]) -> int:
                        first_derivative = np.gradient(energies)
                        second_derivative = np.gradient(first_derivative)

                        # Local maxima (sudden drop)
                        energy_diffs = np.diff(energies)
                        sudden_drop_indices = np.where(energy_diffs < -0.004)[0]

                        # Local maxima (concave down)
                        first_signs = np.sign(first_derivative)
                        cp_indices = np.where(np.diff(first_signs) < 0)[0]

                        # Inflection point
                        second_signs = np.sign(second_derivative)
                        flip_indices = np.where(np.diff(second_signs) < 0)[0]

                        for idx in cp_indices:
                            if second_derivative[idx] < 0:
                                if idx + 1 not in sudden_drop_indices:
                                    return idx

                        for idx in sudden_drop_indices:
                            # Prefer the point AFTER the drop, if concave down
                            if 0 < idx + 1 < len(second_derivative) - 1:
                                if (
                                    second_derivative[idx] < 0
                                    and second_derivative[idx + 1] < 0
                                ):
                                    return idx

                        if len(flip_indices):
                            idx = flip_indices[0]
                            return idx

                        # -------- Fallback: Max energy --------
                        return int(np.argmax(energies))

                    ts_idx = identify_ts_index(energies)

                    selected_step = steps[ts_idx]
                    step_str = str(selected_step).zfill(3)
                    step_xyz_file = next(workdir.rglob(f"calc.{step_str}.xyz"), None)

                    try:
                        xyz_text = read_file_or_none(step_xyz_file)
                    except Exception as e:
                        raise e

                execute_append(
                    "INSERT OR REPLACE INTO xyz (calc_id, selected_step, xyz_text) VALUES (?, ?, ?)",
                    (calc_id, selected_step, xyz_text),
                    db,
                )
        # === Energy from log file ===
        energy_ids = get_existing_calc_ids("energies")
        if calc_id not in energy_ids:
            log_file = next(workdir.rglob("calc.log"), None)
            if log_file:
                log_text = read_file_or_none(log_file)
                if log_text:
                    energy = extract_energy_from_log(log_text)
                    execute_append(
                        "INSERT INTO energies (calc_id, energy_value) VALUES (?, ?)",
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
        if calc_id not in imaginary_ids:
            imag_file = next(workdir.rglob("imaginary.xyz"), None)
            xyz_text = read_file_or_none(imag_file) if imag_file else None

            if xyz_text:
                lines = xyz_text.splitlines()
                frequency = [
                    float(line.split()[-1]) for line in lines if "Frequency" in line
                ]
                execute_append(
                    "INSERT OR REPLACE INTO imaginaryfrequencies (calc_id, frequency, xyz) VALUES (?, ?, ?)",
                    (calc_id, frequency[0], xyz_text),
                    db,
                )

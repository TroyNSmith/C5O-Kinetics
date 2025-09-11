import pandas as pd
from protomol.rd.mol import beta_cleavage, intra_proton_transfer
from protomol.util.ref import get_ccsdt_parameters
from protomol.util.sql import execute_append

from pathlib import Path
import re
import sqlite3
from typing import Optional

db = Path.home() / "C5O-Kinetics/db/data.db"
calc_dir = Path.home() / "C5O-Kinetics/calc"


def write_orca(
    smiles: str,
    method_id: int,
    initial_xyz: str,
    idx1: int = -1,
    idx2: int = -1,
    mechanism: Optional[str] = None,
):
    db = Path.home() / "C5O-Kinetics/db/data.db"
    conn = sqlite3.connect(db)

    # --- Retrieve SMILES entry ---
    df_smiles = pd.read_sql_query(
        "SELECT * FROM smiles WHERE smiles_text = ?",
        conn,
        params=(smiles,),
    )
    if df_smiles.shape[0] != 1:
        raise ValueError(
            f"Expected 1 match in SMILES table, found {df_smiles.shape[0]} for {smiles}"
        )

    smiles_id = int(df_smiles.iloc[0]["smiles_id"])
    multiplicity = df_smiles.iloc[0]["multiplicity"]

    # --- Retrieve method entry ---
    df_methods = pd.read_sql_query(
        "SELECT * FROM methods WHERE method_id = ?",
        conn,
        params=(method_id,),
    )
    if df_methods.shape[0] != 1:
        raise ValueError(
            f"Expected 1 match in METHODS table, found {df_methods.shape[0]} for method_id {method_id}"
        )

    functional = df_methods.iloc[0]["functional"]
    inp_template = df_methods.iloc[0]["inp_template"]
    submit_template = df_methods.iloc[0]["submit_template"]

    # --- Check for existing calculation ---
    df_calculations = pd.read_sql_query(
        "SELECT * FROM calculations WHERE method_id = ? AND smiles_id = ? AND scan_idx1 = ? AND scan_idx2 = ?",
        conn,
        params=(method_id, smiles_id, idx1, idx2),
    )

    if df_calculations.shape[0] > 1:
        raise ValueError(
            f"Expected 0 matches in CALCULATIONS table, found {df_calculations.shape[0]} "
            f"for method_id {method_id}, smiles_id {smiles_id}, scan_idx1 {idx1}, scan_idx2 {idx2}"
        )

    # --- Insert new calculation row ---
    calc_id = int(
        execute_append(
            "INSERT OR REPLACE INTO calculations (smiles_id, method_id, scan_idx1, scan_idx2) VALUES (?, ?, ?, ?)",
            (smiles_id, method_id, idx1, idx2),
            db,
            return_id=True,
        )
    )

    # --- Create calc directory ---
    outdir = Path.home() / f"C5O-Kinetics/calc/{calc_id}"
    outdir.mkdir(exist_ok=True)

    # --- Build substitutions dictionary ---
    substitutions = {
        "[SMILES]": re.sub(r"[^a-zA-Z0-9\s]", "", smiles),
        "[multiplicity]": multiplicity,
        "[idx1]": idx1,
        "[idx2]": idx2,
    }

    if functional == "CCSD(T)-F12/RI":
        substitutions |= get_ccsdt_parameters(smiles.count("C") + smiles.count("O"))

    # --- If a scan is involved, compute scan-specific distances ---
    if mechanism:
        assert idx1 >= 0, f"{idx1} is not a valid index."
        assert idx2 >= 0, f"{idx2} is not a valid index."

        if mechanism.lower() == "proton transfer":
            initial_xyz, max_dist, min_dist = intra_proton_transfer(smiles, idx1, idx2)
        elif mechanism.lower() == "beta cleavage":
            initial_xyz, min_dist, max_dist = beta_cleavage(smiles, idx1, idx2)
        else:
            raise ValueError(f"Unrecognized mechanism {mechanism}")

        substitutions["[start]"] = min_dist
        substitutions["[end]"] = max_dist
    
    # --- Perform substitutions ---
    for key, val in substitutions.items():
        inp_template = inp_template.replace(key, str(val))
        submit_template = submit_template.replace(key, str(val))

    # --- Write input and job files ---
    (outdir / "init.xyz").write_text(initial_xyz)
    (outdir / "calc.inp").write_text(inp_template)
    (outdir / "submit.sh").write_text(submit_template)

    return f"cd {outdir}; sbatch submit.sh"

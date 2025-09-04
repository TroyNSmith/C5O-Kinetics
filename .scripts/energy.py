import pandas as pd

from pathlib import Path
import sqlite3
import sys


def main():
    if len(sys.argv) < 2:
        print("Usage: pixi run find-energy '<SMILES>'")
        sys.exit(1)

    smiles_input = sys.argv[1]

    db_path = Path.home() / "C5O-Kinetics/db/results.db"

    try:
        conn = sqlite3.connect(db_path)
    except sqlite3.Error as e:
        print(f"Failed to connect to database: {e}")
        sys.exit(1)

    try:
        query = "SELECT energy, method, proced FROM calculations WHERE SMILES = ?"
        df = pd.read_sql_query(query, conn, params=(smiles_input,))

        if df.empty:
            print(f"No entries found for SMILES: {smiles_input}")
        else:
            print(f"Energies for SMILES = {smiles_input}:")
            methods = [method[:3] for method in df["method"].values]
            procedures = [proc for proc in df["proced"].values]
            energies = [e for e in df["energy"].values]
            for m, p, e in zip(methods, procedures, energies):
                print(f"{m} {p}: {e} kcal/mol")
    except Exception as e:
        print(f"Error querying database: {e}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()

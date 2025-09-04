from pathlib import Path
import sqlite3


def sqlite3_query(
    execute: str, db: str | Path = Path.home() / "C5O-Kinetics/db/results.db"
):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM calculations")
    rows = cursor.fetchall()
    cursor.close()
    return rows

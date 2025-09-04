#!/bin/bash

DB="$HOME/C5O-Kinetics/db/frequencies.db"

SMILES="$1"
METHOD="$2"
PROCED="$3"
MULT="$4"
MODE="$5"
FREQ="$6"
XYZ="$7"

sqlite3 "$DB" <<EOF
BEGIN;
INSERT OR REPLACE INTO calculations (
  smiles, method, proced, mult, mode, freq, xyz
) VALUES (
  '$SMILES',
  '$METHOD',
  '$PROCED',
  $MULT,
  $MODE,
  $FREQ,
  '$XYZ'
);
COMMIT;
EOF
#!/bin/bash

DB="$HOME/C5O-Kinetics/db/results.db"

SMILES="$1"
METHOD="$2"
PROCEDURE="$3"
MULTIPLICITY="$4"
ENERGY="$5"
XYZ="$6"
HESS="$7"

sqlite3 "$DB" <<EOF
BEGIN;
INSERT OR REPLACE INTO calculations (
  smiles, method, proced, mult, energy, xyz, hess
) VALUES (
  '$SMILES',
  '$METHOD',
  '$PROCEDURE',
  $MULTIPLICITY,
  $ENERGY,
  '$XYZ',
  '$HESS',
);
COMMIT;
EOF
#!/bin/bash

cd $HOME/C5O-Kinetics/db

cat data.dump | sqlite3 refresh.db
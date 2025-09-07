#!/bin/bash

cd $HOME/C5O-Kinetics/db

echo ".dump" | sqlite3 data.db > data.dump
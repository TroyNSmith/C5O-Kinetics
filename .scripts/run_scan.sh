#!/bin/bash

cd $1

if [ ! -f "submit_REVDSD.sh" ]; then
    echo "No submit.sh found, please check $1..."
    exit 1
fi

if [ -f "slurm*" ]; then
    echo "Slurm log file found, please check that this run is not already completed..."
    exit 1
fi

echo "Running from $1..."
sbatch submit_REVDSD.sh
#!/bin/bash

# Slurm Options
#SBATCH --partition=batch
#SBATCH --job-name=trans
#SBATCH --ntasks=8
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G

module load ORCA/6.1

$(which orca) trans.inp > trans.log

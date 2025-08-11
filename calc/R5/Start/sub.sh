#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=test
#SBATCH --ntasks=8
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G

# Import module
module load ORCA/6.1

# Run input file (need to provide absolute path to ORCA when parallelizing)
# Parallel creates overhead that may be more expensive than the initial run
# tail -f {out file} <- Track progress
$(which orca) opt.inp > opt.log
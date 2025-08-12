#!/bin/bash

#SBATCH--partition=batch
#SBATCH --job-name=opt_CdblC
#SBATCH --ntasks=8
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G

module load ORCA/6.1
$(which orca) opt.inp > opt.log
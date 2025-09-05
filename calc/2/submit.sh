#!/bin/bash
#SBATCH --partition=batch
#SBATCH --gres=lscratch:10
#SBATCH --job-name=OO_GOAT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1G

cd ${SLURM_SUBMIT_DIR}
mkdir -p /lscratch/${USER}/${SLURM_JOB_ID}
cp -r "${SLURM_SUBMIT_DIR}/." /lscratch/${USER}/${SLURM_JOB_ID}
cd /lscratch/${USER}/${SLURM_JOB_ID}

module load ORCA/6.1
$(which orca) calc.inp > calc.log
cp -r /lscratch/${USER}/${SLURM_JOB_ID} ${SLURM_SUBMIT_DIR}
rm -rf /lscratch/${USER}/${SLURM_JOB_ID}

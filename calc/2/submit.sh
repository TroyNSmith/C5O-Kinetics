#!/bin/bash
#SBATCH --partition=batch
#SBATCH --gres=lscratch:20
#SBATCH --job-name=CH2CCC1CO1_OPT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=1G

cd ${SLURM_SUBMIT_DIR}
mkdir -p /lscratch/${USER}/${SLURM_JOB_ID}
cp -r "${SLURM_SUBMIT_DIR}/." /lscratch/${USER}/${SLURM_JOB_ID}
cd /lscratch/${USER}/${SLURM_JOB_ID}

module load ORCA/6.1
$(which orca) calc.inp > calc.log
cp -r /lscratch/${USER}/${SLURM_JOB_ID} ${SLURM_SUBMIT_DIR}
rm -rf /lscratch/${USER}/${SLURM_JOB_ID}

PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE SMILES (
    smiles_id INTEGER PRIMARY KEY,
    smiles_text TEXT NOT NULL UNIQUE,
    multiplicity INTEGER NOT NULL,
    initial TEXT NOT NULL
);
CREATE TABLE METHODS (
    method_id INTEGER PRIMARY KEY,
    functional TEXT NOT NULL,
    basis TEXT NOT NULL,
    method TEXT NOT NULL,
    inp_template TEXT NOT NULL,
    submit_template TEXT NOT NULL,
    UNIQUE (functional, basis, method)
);
INSERT INTO METHODS VALUES(1,'REVDSD-PBEP86-D4/2021','def2-TZVPP def2-TZVPP/c','OPT NumFreq',unistr('%PAL NPROCS 16 END\u000a%maxcore 750\u000a! REVDSD-PBEP86-D4/2021 def2-TZVPP def2-TZVPP/c OPT NumFreq\u000a* xyzfile 0 [multiplicity] init.xyz\u000a'),unistr('#!/bin/bash\u000a#SBATCH --partition=batch\u000a#SBATCH --gres=lscratch:20\u000a#SBATCH --job-name=[SMILES]_OPT\u000a#SBATCH --nodes=1\u000a#SBATCH --ntasks=16\u000a#SBATCH --ntasks-per-node=16\u000a#SBATCH --cpus-per-task=1\u000a#SBATCH --time=02:00:00\u000a#SBATCH --mem-per-cpu=1G\u000a\u000acd ${SLURM_SUBMIT_DIR}\u000amkdir -p /lscratch/${USER}/${SLURM_JOB_ID}\u000acp -r "${SLURM_SUBMIT_DIR}/." /lscratch/${USER}/${SLURM_JOB_ID}\u000acd /lscratch/${USER}/${SLURM_JOB_ID}\u000a\u000amodule load ORCA/6.1\u000a$(which orca) calc.inp > calc.log\u000acp -r /lscratch/${USER}/${SLURM_JOB_ID} ${SLURM_SUBMIT_DIR}\u000arm -rf /lscratch/${USER}/${SLURM_JOB_ID}\u000a'));
INSERT INTO METHODS VALUES(2,'XTB','N/A','GOAT',unistr('%PAL NPROCS 8 END\u000a%maxcore 750\u000a! XTB GOAT\u000a* xyzfile 0 [multiplicity] init.xyz\u000a'),unistr('#!/bin/bash\u000a#SBATCH --partition=batch\u000a#SBATCH --gres=lscratch:10\u000a#SBATCH --job-name=[SMILES]_GOAT\u000a#SBATCH --nodes=1\u000a#SBATCH --ntasks=8\u000a#SBATCH --ntasks-per-node=8\u000a#SBATCH --cpus-per-task=1\u000a#SBATCH --time=00:30:00\u000a#SBATCH --mem-per-cpu=1G\u000a\u000acd ${SLURM_SUBMIT_DIR}\u000amkdir -p /lscratch/${USER}/${SLURM_JOB_ID}\u000acp -r "${SLURM_SUBMIT_DIR}/." /lscratch/${USER}/${SLURM_JOB_ID}\u000acd /lscratch/${USER}/${SLURM_JOB_ID}\u000a\u000amodule load ORCA/6.1\u000a$(which orca) calc.inp > calc.log\u000acp -r /lscratch/${USER}/${SLURM_JOB_ID} ${SLURM_SUBMIT_DIR}\u000arm -rf /lscratch/${USER}/${SLURM_JOB_ID}\u000a'));
INSERT INTO METHODS VALUES(3,'CCSD(T)-F12/RI','cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/c','CALC',unistr('%PAL NPROCS [num_cpus] END\u000a%maxcore [input_mem_per_cpu]\u000a! CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/c\u000a* xyzfile 0 [multiplicity] init.xyz\u000a'),unistr('#!/bin/bash\u000a#SBATCH --partition=batch\u000a#SBATCH --gres=lscratch:[lscratch_size]\u000a#SBATCH --job-name=[SMILES]_CALC\u000a#SBATCH --nodes=1\u000a#SBATCH --ntasks=[num_cpus]\u000a#SBATCH --ntasks-per-node=[num_cpus]\u000a#SBATCH --cpus-per-task=1\u000a#SBATCH --time=[time_requested]\u000a#SBATCH --mem-per-cpu=[mem_per_cpu]\u000a\u000acd ${SLURM_SUBMIT_DIR}\u000amkdir -p /lscratch/${USER}/${SLURM_JOB_ID}\u000acp -r "${SLURM_SUBMIT_DIR}/." /lscratch/${USER}/${SLURM_JOB_ID}\u000acd /lscratch/${USER}/${SLURM_JOB_ID}\u000a\u000amodule load ORCA/6.1\u000a$(which orca) calc.inp > calc.log\u000acp -r /lscratch/${USER}/${SLURM_JOB_ID} ${SLURM_SUBMIT_DIR}\u000arm -rf /lscratch/${USER}/${SLURM_JOB_ID}\u000a'));
INSERT INTO METHODS VALUES(4,'XTB','N/A','OPT SCAN',unistr('%PAL NPROCS 8 END\u000a%maxcore 750\u000a! XTB OPT\u000a% geom\u000a   scan\u000a       B [idx1] [idx2] = [start], [end], 18\u000a   end\u000aend\u000a* xyzfile 0 [multiplicity] init.xyz\u000a'),unistr('#!/bin/bash\u000a#SBATCH --partition=batch\u000a#SBATCH --gres=lscratch:10\u000a#SBATCH --job-name=[SMILES]_XTB\u000a#SBATCH --nodes=1\u000a#SBATCH --ntasks=8\u000a#SBATCH --ntasks-per-node=8\u000a#SBATCH --cpus-per-task=1\u000a#SBATCH --time=00:30:00\u000a#SBATCH --mem-per-cpu=1G\u000acd ${SLURM_SUBMIT_DIR}\u000amkdir -p /lscratch/${USER}/${SLURM_JOB_ID}\u000acp -r "${SLURM_SUBMIT_DIR}/." /lscratch/${USER}/${SLURM_JOB_ID}\u000acd /lscratch/${USER}/${SLURM_JOB_ID}\u000a\u000amodule load ORCA/6.1\u000a$(which orca) calc.inp > calc.log\u000acp -r /lscratch/${USER}/${SLURM_JOB_ID} ${SLURM_SUBMIT_DIR}\u000arm -rf /lscratch/${USER}/${SLURM_JOB_ID}\u000a'));
CREATE TABLE CALCULATIONS (
    calc_id INTEGER PRIMARY KEY,
    smiles_id INTEGER REFERENCES SMILES(smiles_id),
    method_id INTEGER REFERENCES Methods(method_id),
    scan_idx1 INTEGER DEFAULT -1,
    scan_idx2 INTEGER DEFAULT -1
);
CREATE TABLE XYZ (
    calc_id INTEGER REFERENCES CALCULATIONS(calc_id) UNIQUE,
    xyz_text TEXT NOT NULL
);
CREATE TABLE TRAJ (
    calc_id INTEGER REFERENCES CALCULATIONS(calc_id) UNIQUE,
    traj_text TEXT NOT NULL
);
CREATE TABLE HESSIANS (
    calc_id INTEGER REFERENCES CALCULATIONS(calc_id) UNIQUE,
    hessian_data TEXT
);
CREATE TABLE IMAGINARYFREQUENCIES (
    calc_id INTEGER REFERENCES CALCULATIONS(calc_id) UNIQUE,
    frequency DOUBLE PRECISION NOT NULL,
    xyz TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS "ENERGIES" (
    calc_id INTEGER REFERENCES CALCULATIONS(calc_id) UNIQUE,
    energy_value DOUBLE PRECISION NOT NULL,
    energy_array TEXT DEFAULT NULL,
    units TEXT DEFAULT 'kcal/mol'
);
COMMIT;
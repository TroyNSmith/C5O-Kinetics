from pathlib import Path


class Optimization:
    def write_inp():
        return

    def write_submission_sh(
        name: str,
        par_dir: str | Path,
        job_opts: str,
        slurm_opts: str = None,
        module_opts: str = None,
        method: str = None,
        next_method: str = None,
        time: str = "24:00:00",
        num_cpus: int = 8,
        mem_per_cpu: int = 16,
    ):
        """
        Writes a SLURM submission shell script to the specified directory.

        Parameters
        ----------
            name (str): The name of the job to be used in the SLURM script.
            par_dir (str | Path): The parent directory where the submission script will be written.
            job_opts (str): The job-specific commands or options to be included in the script.
            slurm_opts (str, optional): Custom SLURM options to override the default SLURM header. If None, a default header is used.
            module_opts (str, optional): Custom module loading commands. If None, a default ORCA 6.1 module load command is used.
            task (str, optional): The specific task being performed (e.g., "opt", "freq", "ts").

        Returns
        -------
            None
        """
        if slurm_opts is None:
            method = "_" + method if method is not None else ""
            partition = "batch" if num_cpus * mem_per_cpu < 480 else "highmem_p"

            slurm_opts = (
                "#!/bin/bash\n"
                f"#SBATCH --partition={partition}\n"
                "#SBATCH --constraint=AMD\n"
                f"#SBATCH --job-name={num_cpus}_cpus_{mem_per_cpu}G_per\n"  # HERE
                "#SBATCH --nodes=1\n"
                f"#SBATCH --ntasks={num_cpus}\n"
                f"#SBATCH --ntasks-per-node={num_cpus}\n"
                "#SBATCH --cpus-per-task=1\n"
                f"#SBATCH --time={time}\n"
                f"#SBATCH --mem-per-cpu={mem_per_cpu}G\n\n"
            )
        if module_opts is None:
            module_opts = "# --- Load ORCA 6.1\n"
            module_opts += "module load ORCA/6.1\n\n"
        if next_method is not None:
            job_opts += (
                f"sbatch --dependency=afterok:$SLURM_JOB_ID submit_{next_method}.sh\n"
            )
        script = slurm_opts + module_opts + job_opts
        path_out = Path(par_dir) / f"submit{method}.sh"
        path_out.write_text(script)

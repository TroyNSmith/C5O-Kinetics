cd .scripts/

job_id=$(sbatch echo.sh | awk '{print $NF}')

echo "Job ID: $job_id"

job_id_2=$(sbatch --dependency=afterok:$job_id echo_2.sh | awk '{print $NF}')

echo "Job ID 2: $job_id_2"
#!/bin/bash
#SBATCH --job-name=calmissed_v3
#SBATCH --partition=production
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err

set -euo pipefail

project_dir=/u/group/c-pionlt/USERS/jainsam/shmsPA
input_file=${project_dir}/out/pion5.422GeV100000.root
result_dir=${project_dir}/v3_results
output_prefix=${result_dir}/pion5p422_v3

mkdir -p "${result_dir}"
cd "${project_dir}"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "Python: $(command -v python3)"

python3 ./CalcMissed_v3.py \
    "${input_file}" \
    5.422 \
    --csv-prefix "${output_prefix}" \
    --json "${output_prefix}_summary.json"

echo "Finished: $(date --iso-8601=seconds)"
echo "Results: ${result_dir}"

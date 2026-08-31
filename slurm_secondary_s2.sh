#!/bin/bash
#SBATCH --job-name=secondary_s2
#SBATCH --partition=production
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err

set -euo pipefail

project_dir=/u/group/c-pionlt/USERS/jainsam/shmsPA
result_dir=${project_dir}/secondary_s2_results

mkdir -p "${result_dir}"
cd "${project_dir}"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "Python: $(command -v python3)"

python3 ./classify_secondary_s2.py \
    out/pion5.422GeV100000.root \
    --output-prefix "${result_dir}/pion5p422_secondary_s2"

echo "Finished: $(date --iso-8601=seconds)"

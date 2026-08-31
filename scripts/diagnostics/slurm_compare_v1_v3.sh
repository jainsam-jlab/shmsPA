#!/bin/bash
#SBATCH --job-name=compare_v1_v3
#SBATCH --partition=production
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err

set -euo pipefail

project_dir=/u/group/c-pionlt/USERS/jainsam/shmsPA
input_file=${project_dir}/out/pion5.422GeV100000.root
output_dir=${project_dir}/v1_v3_diagnostics
analysis_python=/u/group/c-pionlt/USERS/jainsam/replay_lt_env/bin/python3

mkdir -p "${output_dir}" "${project_dir}/slurm_logs"
cd "${project_dir}"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "Python: ${analysis_python}"

"${analysis_python}" scripts/diagnostics/compare_v1_v3_event_selection.py \
    "${input_file}" \
    5.422 \
    --output-dir "${output_dir}"

echo "Finished: $(date --iso-8601=seconds)"
echo "Results: ${output_dir}"

#!/bin/bash
#SBATCH --job-name=secondary_pi_truth
#SBATCH --partition=production
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err

set -euo pipefail
project_dir=/u/group/c-pionlt/USERS/jainsam/shmsPA
root_dir=/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/root/6.24.08-gcc11.4.0
input_file=${project_dir}/out/pion5.422GeV100000.root
output_dir=${project_dir}/secondary_pion_truth_results
output_file=${output_dir}/pion5.422GeV100000_secondary_pion_truth.root
mkdir -p "${output_dir}" "${project_dir}/slurm_logs"
set +u
source "${root_dir}/bin/thisroot.sh"
set -u
cd "${project_dir}"
test -x ./build/diagnostics/secondary_pion_truth_study
echo "Job ${SLURM_JOB_ID} on $(hostname), start $(date --iso-8601=seconds)"
/usr/bin/time -v ./build/diagnostics/secondary_pion_truth_study "${input_file}" "${output_file}" 5.422
/usr/bin/time -v python3 scripts/diagnostics/postprocess_secondary_pion_truth.py \
  "${output_file}" --output-dir "${output_dir}"
echo "Finished $(date --iso-8601=seconds); output ${output_file}"

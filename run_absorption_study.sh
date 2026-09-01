#!/bin/bash
set -euo pipefail

usage() {
  echo "Usage: ./run_absorption_study.sh PARTICLE MOMENTUM_GEV EVENTS"
  echo "Particles: pi+ pi- kaon+ kaon- proton anti_proton e+ e- mu+ mu-"
}
if [[ $# -ne 3 ]]; then usage >&2; exit 2; fi

particle=$1
momentum=$2
events=$3

case "${particle}" in
  pi+) pdg=211; label=pi_plus ;;
  pi-) pdg=-211; label=pi_minus ;;
  kaon+) pdg=321; label=kaon_plus ;;
  kaon-) pdg=-321; label=kaon_minus ;;
  proton) pdg=2212; label=proton ;;
  anti_proton) pdg=-2212; label=anti_proton ;;
  e+) pdg=-11; label=positron ;;
  e-) pdg=11; label=electron ;;
  mu+) pdg=-13; label=mu_plus ;;
  mu-) pdg=13; label=mu_minus ;;
  *) echo "Error: unsupported particle '${particle}'." >&2; usage >&2; exit 2 ;;
esac

if ! [[ "${events}" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: EVENTS must be a positive integer." >&2; exit 2
fi
if ! [[ "${momentum}" =~ ^[0-9]+([.][0-9]+)?$ ]] || [[ "${momentum}" == "0" ]]; then
  echo "Error: MOMENTUM_GEV must be a positive number." >&2; exit 2
fi
project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
kinematics="${project_dir}/run5055kine.txt"
if [[ ! -r "${kinematics}" ]]; then
  echo "Error: cannot read '${kinematics}'." >&2; exit 2
fi
kinematic_events=$(awk 'NF && $1 !~ /^#/' "${kinematics}" | wc -l)
if (( kinematic_events < events )); then
  echo "Error: requested ${events} events, but '${kinematics}' has only ${kinematic_events} data rows." >&2
  exit 2
fi

run_name="${label}_${momentum}GeV_${events}"
output_dir="${project_dir}/absorption_results/${run_name}"
macro_file="${output_dir}/run.mac"
geant_output="${output_dir}/geant4.root"
truth_output="${output_dir}/absorption_truth.root"
report_dir="${output_dir}/report"
mkdir -p "${output_dir}" "${report_dir}"
rm -f \
  "${report_dir}/summary.json" \
  "${report_dir}/ABSORPTION_STUDY.md" \
  "${report_dir}/SECONDARY_PION_VISIBILITY_STUDY.md" \
  "${report_dir}/AGC_vs_HGC_no_target_in_PID.png" \
  "${report_dir}/AGC_vs_HGC_no_pion_in_PID.png"

cat > "${macro_file}" <<EOF
## Macro for SHMS particle material-loss simulation

/run/initialize

## Use focal-plane kinematics
/PA/generator/useGenerated true
/PA/generator/setInFile run5055kine.txt

## Pion and kaon decay are already included in SIMC
/process/inactivate Decay pi+
/process/inactivate Decay kaon+
/run/physicsModified

## Central SHMS momentum
/PA/generator/momentum ${momentum} GeV

## Primary particle
/gun/particle ${particle}

## ROOT output
/analysis/setFileName ${geant_output}

## Number of generated events
/run/printProgress 100
/run/beamOn ${events}
EOF

if ! command -v geant4-config >/dev/null 2>&1; then
  echo "Error: geant4-config is unavailable. Source the farm Geant4 environment first." >&2
  exit 1
fi
if ! command -v root-config >/dev/null 2>&1; then
  echo "Error: root-config is unavailable. Source the farm ROOT environment first." >&2
  exit 1
fi
if [[ -z "${GEANT4_DATA_DIR:-}" || ! -d "${GEANT4_DATA_DIR}" ]]; then
  echo "Error: GEANT4_DATA_DIR is not set to a valid directory." >&2
  exit 1
fi
if ! python3 -c 'import ROOT' >/dev/null 2>&1; then
  echo "Error: Python cannot import ROOT. Source the farm ROOT Python environment first." >&2
  exit 1
fi
echo "Study ${run_name}: ${particle} (PDG ${pdg}), ${momentum} GeV, ${events} events"
echo "Output: ${output_dir}"
cd "${project_dir}"

if [[ ! -f "${project_dir}/build/CMakeCache.txt" ]]; then
  cmake -S "${project_dir}" -B "${project_dir}/build"
fi
cmake --build "${project_dir}/build" -j4
"${project_dir}/scripts/diagnostics/build_secondary_pion_truth_study.sh"

echo "Running Geant4..."
"${project_dir}/build/shmsPA" "${macro_file}"
test -s "${geant_output}"

echo "Running direct-daughter truth analysis..."
"${project_dir}/build/diagnostics/secondary_pion_truth_study" \
  "${geant_output}" "${truth_output}" "${momentum}" "${pdg}"
test -s "${truth_output}"

echo "Creating absorption summary and plots..."
python3 "${project_dir}/scripts/diagnostics/postprocess_secondary_pion_truth.py" \
  "${truth_output}" --output-dir "${report_dir}" --particle-label "${particle}"

echo "Study complete."
echo "Geant4 ROOT:  ${geant_output}"
echo "Truth ROOT:   ${truth_output}"
echo "Summary CSV:  ${report_dir}/absorption_summary.csv"
echo "Plots:        ${report_dir}"

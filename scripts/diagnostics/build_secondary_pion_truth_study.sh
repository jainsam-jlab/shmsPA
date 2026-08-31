#!/bin/bash
set -euo pipefail
project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
if ! command -v root-config >/dev/null 2>&1; then
  echo "Error: root-config is unavailable. Source the farm ROOT environment first." >&2
  exit 1
fi
mkdir -p "${project_dir}/build/diagnostics"
root-config --cxx --cflags --libs > /dev/null
"$(root-config --cxx)" -O3 -DNDEBUG -std=c++17 \
  $(root-config --cflags) \
  "${project_dir}/scripts/diagnostics/secondary_pion_truth_study.cpp" \
  -o "${project_dir}/build/diagnostics/secondary_pion_truth_study" \
  $(root-config --libs) -lROOTNTuple -lvdt
echo "Built ${project_dir}/build/diagnostics/secondary_pion_truth_study"

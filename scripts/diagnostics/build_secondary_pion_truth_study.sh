#!/bin/bash
set -euo pipefail
project_dir=/u/group/c-pionlt/USERS/jainsam/shmsPA
root_dir=/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/root/6.24.08-gcc11.4.0
mkdir -p "${project_dir}/build/diagnostics"
"${root_dir}/bin/root-config" --cxx --cflags --libs > /dev/null
"$("${root_dir}/bin/root-config" --cxx)" -O3 -DNDEBUG -std=c++17 \
  $("${root_dir}/bin/root-config" --cflags) \
  "${project_dir}/scripts/diagnostics/secondary_pion_truth_study.cpp" \
  -o "${project_dir}/build/diagnostics/secondary_pion_truth_study" \
  $("${root_dir}/bin/root-config" --libs) -lROOTNTuple -lvdt
echo "Built ${project_dir}/build/diagnostics/secondary_pion_truth_study"

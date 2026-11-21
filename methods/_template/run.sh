#!/usr/bin/env bash
# Usage: run.sh <dataset_dir> <out_dir> <params_json> <threads>
# Tool paths arrive as POLYBENCH_<TOOL> environment variables (e.g. POLYBENCH_ASTRAL).
# POLYBENCH_PYTHON is the python that runs polybench.
set -euo pipefail
dataset_dir=$1; out_dir=$2; params_json=$3; threads=$4

# Replace this line with the call to the real tool. The only requirement is that
# it writes $out_dir/network.nwk.
cp "$dataset_dir/species_tree.nwk" "$out_dir/network.nwk"

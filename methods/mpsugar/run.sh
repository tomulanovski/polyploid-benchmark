#!/usr/bin/env bash
# Usage: run.sh <dataset_dir> <out_dir> <params_json> <threads>
# PhyNetPy is a Python library: run_mpsugar.py runs with `python` of the method env (method.yaml env),
# the output is parsed with $POLYBENCH_PYTHON.
set -euo pipefail
dataset_dir=$1; out_dir=$2; params_json=$3; threads=$4
inputs="$out_dir/work/inputs"
run_dir="$out_dir/work/mpsugar"

mkdir -p "$run_dir"
python "$(dirname "${BASH_SOURCE[0]}")/run_mpsugar.py" \
    "$inputs/mpsugar_trees.nex" "$inputs/taxon_map.json" "$params_json" "$run_dir/mpsugar_results.txt"
"$POLYBENCH_PYTHON" -m polybench.prep.outputs mpsugar "$run_dir/mpsugar_results.txt" > "$out_dir/network.nwk"

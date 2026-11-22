#!/usr/bin/env bash
# Usage: run.sh <dataset_dir> <out_dir> <params_json> <threads>
# GRANDMA writes the final MUL-tree to <out>/final_multree.tre.
set -euo pipefail
dataset_dir=$1; out_dir=$2; params_json=$3; threads=$4
py="$POLYBENCH_PYTHON"
inputs="$out_dir/work/inputs"
run_dir="$out_dir/work/grandma"
flags_txt=$("$py" -m polybench.prep.params flags "$params_json" "$(dirname "${BASH_SOURCE[0]}")/method.yaml")
flags=()
if [[ -n "$flags_txt" ]]; then mapfile -t flags <<< "$flags_txt"; fi

mkdir -p "$run_dir"
"$py" "${POLYBENCH_GRANDMA:?tools.grandma is not set in the config}" \
    -g "$inputs/grampa_trees.tre" -s "$inputs/species_tree.tre" -o "$run_dir" \
    "${flags[@]}" -p "$threads" -v 3 --overwrite
"$py" -m polybench.prep.outputs grandma_split "$run_dir" "$inputs/taxa_map.tsv" > "$out_dir/network.nwk"

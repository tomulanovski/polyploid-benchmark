#!/usr/bin/env bash
# Usage: run.sh <dataset_dir> <out_dir> <params_json> <threads>
# Polyphest writes <output_dir>/polyphest_trees-polyphest.txt.
set -euo pipefail
dataset_dir=$1; out_dir=$2; params_json=$3; threads=$4
inputs="$out_dir/work/inputs"
run_dir="$out_dir/work/polyphest"
polyphest="${POLYBENCH_POLYPHEST:?tools.polyphest is not set in the config}"
flags_txt=$("$POLYBENCH_PYTHON" -m polybench.prep.params flags "$params_json" "$(dirname "${BASH_SOURCE[0]}")/method.yaml")
flags=()
if [[ -n "$flags_txt" ]]; then mapfile -t flags <<< "$flags_txt"; fi

mkdir -p "$run_dir"
python "$polyphest" \
    --gene_tree_file "$inputs/polyphest_trees.tre" \
    --consensus_multiset_file "$inputs/multi_set.txt" \
    --output_dir "$run_dir" \
    "${flags[@]}"
"$POLYBENCH_PYTHON" -m polybench.prep.outputs polyphest "$run_dir/polyphest_trees-polyphest.txt" > "$out_dir/network.nwk"

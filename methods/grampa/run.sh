#!/usr/bin/env bash
# Usage: run.sh <dataset_dir> <out_dir> <params_json> <threads>
# mode iter_first runs GRANDMA once (-m single) and keeps iteration (0, 0) of history.json,
# mode original runs standalone GRAMPA and keeps the best row of grampa-scores.txt.
set -euo pipefail
dataset_dir=$1; out_dir=$2; params_json=$3; threads=$4
py="$POLYBENCH_PYTHON"
inputs="$out_dir/work/inputs"
run_dir="$out_dir/work/grampa"
mode=$("$py" -m polybench.prep.params get "$params_json" mode)

mkdir -p "$run_dir"
case "$mode" in
    iter_first)
        cap=$("$py" -m polybench.prep.params get "$params_json" cap)
        "$py" "${POLYBENCH_GRANDMA:?tools.grandma is not set in the config}" \
            -g "$inputs/grampa_trees.tre" -s "$inputs/species_tree.tre" -o "$run_dir" \
            --mode single --cap "$cap" -p "$threads" -v 3 --overwrite
        "$py" -m polybench.prep.outputs grandma_first_iteration "$run_dir/history.json" "$inputs/taxa_map.tsv" \
            > "$out_dir/network.nwk"
        ;;
    original)
        "${POLYBENCH_GRAMPA:?tools.grampa is not set in the config}" \
            -g "$inputs/grampa_trees.tre" -s "$inputs/species_tree.tre" -o "$run_dir" --overwrite
        "$py" -m polybench.prep.outputs grampa_scores "$run_dir/grampa-scores.txt" "$inputs/taxa_map.tsv" \
            > "$out_dir/network.nwk"
        ;;
    *)
        echo "grampa: unknown mode '$mode', expected iter_first or original" >&2
        exit 1
        ;;
esac

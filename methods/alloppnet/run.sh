#!/usr/bin/env bash
# Usage: run.sh <dataset_dir> <out_dir> <params_json> <threads>
# Rscript (method env) writes the BEAST XML from the AlloppDT scripts, BEAST samples MUL-trees until the chain
# ends or beast_hours runs out, TreeAnnotator summarizes them as the MCC tree with 10% burn-in.
set -euo pipefail
dataset_dir=$1; out_dir=$2; params_json=$3; threads=$4
py="$POLYBENCH_PYTHON"
inputs="$out_dir/work/inputs"
run_dir="$out_dir/work/alloppnet"
alloppdt="${POLYBENCH_ALLOPPDT:?tools.alloppdt is not set in the config}"
beast="${POLYBENCH_BEAST:?tools.beast is not set in the config}"
treeannotator="${POLYBENCH_TREEANNOTATOR:?tools.treeannotator is not set in the config}"
param() { "$py" -m polybench.prep.params get "$params_json" "$1"; }

mkdir -p "$run_dir"
Rscript "$(dirname "${BASH_SOURCE[0]}")/generate_beast_xml.r" "$alloppdt" "$inputs" "$run_dir" \
    "$(param chain_length)" "$(param log_every)" "$(param screen_log_every)"

# timeout exits 124 when it stops BEAST (killed 5 minutes later if still alive), the trees sampled so far
# are kept, as in the paper runs
code=0
(cd "$run_dir" && timeout -k 5m "$(param beast_hours)h" "$beast" alloppnet.XML) || code=$?
if [[ $code -ne 0 && $code -ne 124 ]]; then
    echo "alloppnet: BEAST failed with exit code $code" >&2
    exit "$code"
fi

trees="$run_dir/sampledmultrees.txt"
[[ -s "$trees" ]] || { echo "alloppnet: BEAST wrote no $trees" >&2; exit 1; }
# a stopped BEAST can leave a half-written last tree
last=$(tail -n 1 "$trees" | tr -d '\r')
if [[ $last == *"tree STATE_"* && $last != *";" ]]; then
    sed -i '$d' "$trees"
    echo "alloppnet: dropped the incomplete last tree of $trees"
fi
total=$(grep -c "tree STATE_" "$trees" || true)
if [[ $code -eq 124 ]]; then
    echo "alloppnet: BEAST exit code 124 (stopped after beast_hours), $total sampled trees"
else
    echo "alloppnet: BEAST finished, $total sampled trees"
fi
if (( total < 2 )); then
    echo "alloppnet: $trees has $total sampled trees, need at least 2" >&2
    exit 1
fi
burnin=$(( (total - 1) / 10 ))  # 10% of the trees after STATE_0
"$treeannotator" -burninTrees "$burnin" -heights mean "$trees" "$run_dir/mcc.tre"
"$py" -m polybench.prep.outputs alloppnet_mcc "$run_dir/mcc.tre" "$dataset_dir/ploidy.tsv" > "$out_dir/network.nwk"

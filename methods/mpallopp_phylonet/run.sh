#!/usr/bin/env bash
# Usage: run.sh <dataset_dir> <out_dir> <params_json> <threads>
# PhyloNet runs the command block of the NEXUS file and prints the inferred network to stdout.
set -euo pipefail
dataset_dir=$1; out_dir=$2; params_json=$3; threads=$4
inputs="$out_dir/work/inputs"
jar="${POLYBENCH_PHYLONET_JAR:?tools.phylonet_jar is not set in the config}"

sed "s/{threads}/$threads/g" "$inputs/mpallopp.nex" > "$out_dir/work/mpallopp_run.nex"
java -jar "$jar" "$out_dir/work/mpallopp_run.nex" | tee "$out_dir/work/phylonet_out.txt"
"$POLYBENCH_PYTHON" -m polybench.prep.outputs phylonet "$out_dir/work/phylonet_out.txt" > "$out_dir/network.nwk"

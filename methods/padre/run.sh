#!/usr/bin/env bash
# Usage: run.sh <dataset_dir> <out_dir> <params_json> <threads>
# PADRE writes <input>-result.tre next to its input.
set -euo pipefail
dataset_dir=$1; out_dir=$2; params_json=$3; threads=$4
inputs="$out_dir/work/inputs"
jar="${POLYBENCH_PADRE_JAR:?tools.padre_jar is not set in the config}"
java_mem=$("$POLYBENCH_PYTHON" -m polybench.prep.params get "$params_json" java_mem)

cd "$inputs"
java -Djava.awt.headless=true -Xmx"$java_mem" -jar "$jar" -i "$inputs/padre_trees.tre"
"$POLYBENCH_PYTHON" -m polybench.prep.outputs padre "$inputs/padre_trees-result.tre" > "$out_dir/network.nwk"

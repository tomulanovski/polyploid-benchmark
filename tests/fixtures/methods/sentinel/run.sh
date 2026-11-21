#!/usr/bin/env bash
# fails when a previous run left its tool state in work/
set -euo pipefail
if [[ -e "$2/work/sentinel" ]]; then echo "stale work/sentinel" >&2; exit 1; fi
touch "$2/work/sentinel"
cp "$1/species_tree.nwk" "$2/network.nwk"

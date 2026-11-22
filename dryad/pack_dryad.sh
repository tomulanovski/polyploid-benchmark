#!/usr/bin/env bash
# Usage: dryad/pack_dryad.sh --gene2net-root <gene2net> --out <dryad_dir> [--dry-run] [--configs ...] [--no-archives]
# Runs with the polybench env's python ($POLYBENCH_PYTHON, else python). See dryad/README_template.md.
set -euo pipefail
exec "${POLYBENCH_PYTHON:-python}" "$(dirname "$0")/pack_dryad.py" "$@"

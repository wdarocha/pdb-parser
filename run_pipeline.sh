#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 5 ]; then
	cat >&2 <<'USAGE'
Usage:
  ./run_pipeline.sh <pdb_ids.dat> <parser_params.cfg> <instance_reorder.cfg> <pdb_data_dir> <out_dir>

Arguments:
  pdb_ids.dat            File with one PDB id per line
  parser_params.cfg      Parameter file for the parser stage
  instance_reorder.cfg   Parameter file for the instance reordering stage
  pdb_data_dir           Directory where downloaded PDB files are stored
  out_dir                Directory where parser and reordered outputs are written
USAGE
	exit 2
fi

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"

PDB_IDS_FILE="$1"
PARSER_PARAMS_FILE="$2"
REORDER_PARAMS_FILE="$3"
PDB_DATA_DIR="$4"
OUT_DIR="$5"
PYTHON_BIN="${PYTHON:-python3}"

export PYTHONPATH="$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

printf '======================================================\n'
printf 'Stage 1/2: parser\n'
printf '======================================================\n'
"$PYTHON_BIN" -m pdb_parser.pdb_parser \
	"$PDB_IDS_FILE" \
	"$PARSER_PARAMS_FILE" \
	"$PDB_DATA_DIR" \
	"$OUT_DIR"

printf '======================================================\n'
printf 'Stage 2/2: instance reordering\n'
printf '======================================================\n'
"$PYTHON_BIN" -m pdb_parser.instance_reorder \
	"$PDB_IDS_FILE" \
	"$REORDER_PARAMS_FILE" \
	"$OUT_DIR"

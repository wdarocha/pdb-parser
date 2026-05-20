#!/usr/bin/env bash
set -euo pipefail

usage() {
	cat >&2 <<'USAGE'
Usage:
  ./run_pipeline.sh [options] <pdb_ids.dat> <parser_params.cfg> <instance_reorder.cfg> <pdb_data_dir> <seed_dir> <out_dir>

Options:
  --mode <mode>       Execution mode: full, parser, reorder, or reorder-all (default: full)
  --orders <list>     Order ids for reorder-all, e.g. 1-10 or 1,3,9 (default: 1-10)
  --log-dir <dir>     Write one log file per stage to this directory
  -h, --help          Show this help message

Arguments:
  pdb_ids.dat            File with one PDB id per line
  parser_params.cfg      Parameter file for the parser stage
  instance_reorder.cfg   Base parameter file for the instance reordering stage
  pdb_data_dir           Directory where downloaded PDB files are stored
  seed_dir               Directory where per-PDB seed files are stored
  out_dir                Directory where parser and reordered outputs are written

Modes:
  full         Run parser and then instance reordering with the configured order_id
  parser       Run only the parser stage
  reorder      Run only instance reordering with the configured order_id
  reorder-all  Run only instance reordering for all order ids selected by --orders
USAGE
}

MODE="full"
ORDERS="1-10"
LOG_DIR=""

while [ "$#" -gt 0 ]; do
	case "$1" in
		--mode)
			MODE="${2:?missing value for --mode}"
			shift 2
			;;
		--orders)
			ORDERS="${2:?missing value for --orders}"
			shift 2
			;;
		--log-dir)
			LOG_DIR="${2:?missing value for --log-dir}"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		--)
			shift
			break
			;;
		-*)
			printf 'Unknown option: %s

' "$1" >&2
			usage
			exit 2
			;;
		*)
			break
			;;
	esac
done

case "$MODE" in
	full|parser|reorder|reorder-all) ;;
	*)
		printf 'Invalid mode: %s

' "$MODE" >&2
		usage
		exit 2
		;;
esac

if [ "$#" -ne 6 ]; then
	usage
	exit 2
fi

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"

PDB_IDS_FILE="$1"
PARSER_PARAMS_FILE="$2"
REORDER_PARAMS_FILE="$3"
PDB_DATA_DIR="$4"
SEED_DIR="$5"
OUT_DIR="$6"
PYTHON_BIN="${PYTHON:-python3}"

export PYTHONPATH="$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

if [ -n "$LOG_DIR" ]; then
	mkdir -p "$LOG_DIR"
fi

run_with_optional_log() {
	local log_name="$1"
	shift

	if [ -n "$LOG_DIR" ]; then
		"$@" > "$LOG_DIR/$log_name" 2>&1
	else
		"$@"
	fi
}

print_header() {
	printf '======================================================
'
	printf '%s
' "$1"
	printf '======================================================
'
}

run_parser() {
	print_header 'Stage 1/2: parser'
	run_with_optional_log parser.log 		"$PYTHON_BIN" -m pdb_parser.pdb_parser 		"$PDB_IDS_FILE" 		"$PARSER_PARAMS_FILE" 		"$PDB_DATA_DIR" 		"$SEED_DIR" 		"$OUT_DIR"
}

run_reorder() {
	local params_file="$1"
	local log_name="$2"

	run_with_optional_log "$log_name" 		"$PYTHON_BIN" -m pdb_parser.instance_reorder 		"$PDB_IDS_FILE" 		"$params_file" 		"$OUT_DIR"
}

expand_orders() {
	local spec="$1"
	local item start end order
	IFS=',' read -ra items <<< "$spec"
	for item in "${items[@]}"; do
		if [[ "$item" =~ ^[0-9]+-[0-9]+$ ]]; then
			start="${item%-*}"
			end="${item#*-}"
			if [ "$start" -gt "$end" ]; then
				printf 'Invalid order range: %s
' "$item" >&2
				exit 2
			fi
			for ((order=start; order<=end; order++)); do
				printf '%s
' "$order"
			done
		elif [[ "$item" =~ ^[0-9]+$ ]]; then
			printf '%s
' "$item"
		else
			printf 'Invalid order specification: %s
' "$item" >&2
			exit 2
		fi
	done
}

write_reorder_params_for_order() {
	local order_id="$1"
	local source_file="$2"
	local target_file="$3"

	awk -v order_id="$order_id" '
		BEGIN { replaced = 0 }
		/^[[:space:]]*order_id[[:space:]]*:/ {
			print "order_id: " order_id
			replaced = 1
			next
		}
		{ print }
		END {
			if (!replaced) {
				print "order_id: " order_id
			}
		}
	' "$source_file" > "$target_file"
}

case "$MODE" in
	full)
		run_parser
		print_header 'Stage 2/2: instance reordering'
		run_reorder "$REORDER_PARAMS_FILE" reorder.log
		;;
	parser)
		run_parser
		;;
	reorder)
		print_header 'Stage 2/2: instance reordering'
		run_reorder "$REORDER_PARAMS_FILE" reorder.log
		;;
	reorder-all)
		TMP_DIR="$(mktemp -d)"
		trap 'rm -rf "$TMP_DIR"' EXIT
		print_header "Stage 2/2: instance reordering for order_id=${ORDERS}"
		expand_orders "$ORDERS" | while read -r ORDER_ID; do
			TMP_REORDER_PARAMS="$TMP_DIR/instance_reorder_order${ORDER_ID}.cfg"
			write_reorder_params_for_order "$ORDER_ID" "$REORDER_PARAMS_FILE" "$TMP_REORDER_PARAMS"
			printf '%s
' '------------------------------------------------------'
			printf 'Running order_id=%s
' "$ORDER_ID"
			printf '%s
' '------------------------------------------------------'
			run_reorder "$TMP_REORDER_PARAMS" "reorder_order${ORDER_ID}.log"
		done
		;;
esac

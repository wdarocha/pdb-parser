"""CLI for parsing PDB structures and writing output files.

Usage:
	python -m pdb_parser.pdb_parser <pdb_ids.dat> <params.cfg> <pdb_data_dir> <out_dir>

Arguments
---------
pdb_ids.dat	: text file with one PDB id per line (comments with '#')
params.cfg	: simple key=value file (comments with '#')
pdb_data_dir	: pdb data directory to be created if missing
out_dir		: output directory to be created if missing
"""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

if __package__ in {None, ""}:
	sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pdb_parser.io import download_pdb
from pdb_parser.pipeline import parser
from pdb_parser.utils import ensure_dir, read_params, read_pdb_ids


# -----------------------------------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
	p = argparse.ArgumentParser(
		prog="pdb-parser",
		description="Download and parse PDB structures, writing the output files to disk.",
	)
	p.add_argument("pdb_ids_file", type=Path, help="Path to a file with one PDB id per line.")
	p.add_argument("params_file", type=Path, help="Path to a simple key=value params file.")
	p.add_argument("pdb_data_dir", type=Path, help="PDB data directory to create if missing.")
	p.add_argument("out_dir", type=Path, help="Output directory to create if missing.")
	return p


def print_job_banner(pdb_id: str) -> None:
	print("------------------------------------------------------")
	print(f"------------------------ {pdb_id} ------------------------")
	print("------------------------------------------------------")


# -----------------------------------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
	arg_parser = build_parser()
	args = arg_parser.parse_args(argv)

	try:
		pdb_ids = read_pdb_ids(args.pdb_ids_file)
		params = read_params(args.params_file)
		pdb_data_dir = ensure_dir(args.pdb_data_dir)
		out_dir = ensure_dir(args.out_dir)
	except Exception as exc:
		print(f"[ERROR] {exc}", file=sys.stderr)
		return 2

	for pdb_id in pdb_ids:
		print_job_banner(pdb_id)
		try:
			download_pdb(pdb_id, pdb_data_dir)
			parser(params, pdb_data_dir, out_dir, pdb_id, True)
		except Exception as exc:
			print(f"[{pdb_id}] Skipped: {exc}")
			continue

	return 0


# -----------------------------------------------------------------------------------------------------
if __name__ == "__main__":
	raise SystemExit(main())

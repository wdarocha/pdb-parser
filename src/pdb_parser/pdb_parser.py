"""CLI for parsing PDB structures and writing output files.

The positional ``params.cfg`` input is a UTF-8 text file containing one
``key: value`` pair per line. The command processes one selected chain per
PDB execution and continues after an ordinary per-structure error.

Usage:
	python -m pdb_parser.pdb_parser [--keep-tmp] <pdb_ids.dat> <params.cfg> <pdb_data_dir> <seed_dir> <out_dir>

Arguments
---------
pdb_ids.dat	: text file with one PDB id per line (comments with '#')
params.cfg	: simple key: value file (comments with '#')
pdb_data_dir	: pdb data directory to be created if missing
seed_dir	: seed directory to be created if missing
out_dir		: output directory to be created if missing
--keep-tmp	: keep the intermediate files in <out_dir>/<pdb_id>/tmp/ (removed by default)

Configuration files use one ``key: value`` pair per line; blank lines and
comments beginning with ``#`` are ignored.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import TypedDict

if __package__ in {None, ""}:
	source_root = str(Path(__file__).resolve().parents[1])
	if source_root not in sys.path:
		sys.path.insert(0, source_root)

from pdb_parser.api import run_parser_job
from pdb_parser.cli_contexts import (
	ensure_public_parser_dirs,
	load_public_parser_params,
	load_public_pdb_ids,
)
from pdb_parser.cli_support import run_cli_command
from pdb_parser.pipeline import parser


class ParserCliContext(TypedDict):
	pdb_ids: list[str]
	params: dict[str, str]
	pdb_data_dir: Path
	seed_dir: Path
	out_dir: Path
	keep_tmp: bool


def build_parser() -> argparse.ArgumentParser:
	"""Build the command-line argument parser for the parser stage."""
	p = argparse.ArgumentParser(
		prog="pdb-parser",
		description="Download and parse PDB structures, writing the output files to disk.",
	)
	p.add_argument("pdb_ids_file", type=Path, help="Path to a file with one PDB id per line.")
	p.add_argument("params_file", type=Path, help="Path to a simple key: value params file.")
	p.add_argument("pdb_data_dir", type=Path, help="PDB data directory to create if missing.")
	p.add_argument("seed_dir", type=Path, help="Seed directory to create if missing.")
	p.add_argument("out_dir", type=Path, help="Output directory to create if missing.")
	p.add_argument(
		"--keep-tmp",
		action="store_true",
		help=(
			"Keep the intermediate files in <out_dir>/<pdb_id>/tmp/ after each structure "
			"(removed by default; always removed when a structure fails)."
		),
	)
	return p


def load_parser_cli_context(args: argparse.Namespace) -> ParserCliContext:
	"""Load and validate the public parser CLI inputs before processing jobs."""
	params = load_public_parser_params(args.params_file)
	pdb_data_dir, seed_dir, out_dir = ensure_public_parser_dirs(
		args.pdb_data_dir,
		args.seed_dir,
		args.out_dir,
	)

	return {
		"pdb_ids": load_public_pdb_ids(args.pdb_ids_file),
		"params": params,
		"pdb_data_dir": pdb_data_dir,
		"seed_dir": seed_dir,
		"out_dir": out_dir,
		"keep_tmp": bool(getattr(args, "keep_tmp", False)),
	}


# -----------------------------------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
	"""Run the parser CLI and return a process-style exit code."""
	def build_job_runner(context: ParserCliContext):
		def run_job(pdb_id: str) -> None:
			run_parser_job(
				context["params"],
				context["pdb_data_dir"],
				context["seed_dir"],
				context["out_dir"],
				pdb_id,
				parser,
				remove_tmp_dir=not context["keep_tmp"],
			)

		return run_job

	return run_cli_command(argv, build_parser, load_parser_cli_context, build_job_runner)


# -----------------------------------------------------------------------------------------------------
if __name__ == "__main__":
	raise SystemExit(main())

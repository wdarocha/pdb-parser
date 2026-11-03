"""CLI for remapping previously written parser outputs.

The positional ``instance_remap.cfg`` input is a UTF-8 text file containing
one ``key: value`` pair per line. Remapping validates continuity independently
from the parser and accepts only the public method/atom-selection combinations
documented in ``docs/remap_stage.md``.

Usage:
	python -m pdb_parser.instance_remap <pdb_ids.dat> <instance_remap.cfg> <out_dir>

Arguments
---------
pdb_ids.dat		: text file with one PDB id per line (comments with '#')
instance_remap.cfg	: simple key: value file for the remapping stage
out_dir			: output directory containing files written by pdb_parser.py

The remapping configuration uses one ``key: value`` pair per line; blank
lines and comments beginning with ``#`` are ignored. ``order_id`` accepts
the fixed values ``1`` through ``10``, alone or as a list such as ``1-10`` or
``1,3,9``; ``order_type`` accepts ``ddgp``, ``redmdgp`` or both, comma-separated.
Every combination is written from one reading of the parser files. For ``backbone``
outputs, ``order_type: dmdgp`` writes the fixed DMDGP order ``N CA C`` without
``order_id``; it cannot be combined with the other types.
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

from pdb_parser.api import run_remap_jobs
from pdb_parser.cli_contexts import (
	ensure_public_remap_out_dir,
	load_public_pdb_ids,
	load_public_remap_order_choices,
	load_public_remap_order_types,
	load_public_remap_params,
)
from pdb_parser.cli_support import run_cli_command
from pdb_parser.remap_cli_text import (
	REMAP_COMMAND,
	REMAP_DESCRIPTION,
	REMAP_OUT_DIR_HELP,
	REMAP_PARAMS_ARGUMENT,
	REMAP_PARAMS_HELP,
)


class RemapCliContext(TypedDict):
	pdb_ids: list[str]
	remap_params: dict[str, str]
	order_choices: list[int]
	order_types: list[str]
	out_dir: Path


def build_parser() -> argparse.ArgumentParser:
	"""Build the command-line argument parser for the remapping stage."""
	parser = argparse.ArgumentParser(
		prog=REMAP_COMMAND,
		description=REMAP_DESCRIPTION,
	)
	parser.add_argument("pdb_ids_file", type=Path, help="Path to a file with one PDB id per line.")
	parser.add_argument(
		REMAP_PARAMS_ARGUMENT,
		type=Path,
		help=REMAP_PARAMS_HELP,
	)
	parser.add_argument("out_dir", type=Path, help=REMAP_OUT_DIR_HELP)
	return parser


def read_instance_remap_params(path: Path) -> dict[str, str]:
	"""Read and validate the public configuration for the remapping stage."""
	return load_public_remap_params(path)


def load_remap_cli_context(args: argparse.Namespace) -> RemapCliContext:
	"""Load and validate the public remapping CLI inputs before processing jobs."""
	out_dir = ensure_public_remap_out_dir(args.out_dir)
	remap_params = read_instance_remap_params(args.remap_params_file)
	return {
		"pdb_ids": load_public_pdb_ids(args.pdb_ids_file),
		"remap_params": remap_params,
		"order_choices": load_public_remap_order_choices(remap_params),
		"order_types": load_public_remap_order_types(remap_params),
		"out_dir": out_dir,
	}


def main(argv: list[str] | None = None) -> int:
	"""Run the remapping CLI and return a process-style exit code."""
	def build_job_runner(context: RemapCliContext):
		def run_job(pdb_id: str) -> None:
			run_remap_jobs(
				context["remap_params"],
				context["out_dir"],
				pdb_id,
				context["order_choices"],
				context["order_types"],
			)

		return run_job

	return run_cli_command(argv, build_parser, load_remap_cli_context, build_job_runner)


if __name__ == "__main__":
	raise SystemExit(main())

"""CLI for reordering previously written parser outputs.

Usage:
	python -m pdb_parser.instance_reorder <pdb_ids.dat> <instance_reorder.cfg> <out_dir>

Arguments
---------
pdb_ids.dat		: text file with one PDB id per line (comments with '#')
instance_reorder.cfg	: simple key=value file for the reordering stage
out_dir			: output directory containing files written by pdb_parser.py
"""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

import numpy as np

if __package__ in {None, ""}:
	sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pdb_parser.io import read_space_separated_file
from pdb_parser.pipeline import reorder_instance
from pdb_parser.utils import read_params, read_pdb_ids


SUPPORTED_ORDER_IDS = {1, 9}


# -----------------------------------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
	p = argparse.ArgumentParser(
		prog="instance-reorder",
		description="Reorder previously generated parser output files.",
	)
	p.add_argument("pdb_ids_file", type=Path, help="Path to a file with one PDB id per line.")
	p.add_argument(
		"reorder_params_file",
		type=Path,
		help="Path to a simple key=value params file for the reordering stage.",
	)
	p.add_argument("out_dir", type=Path, help="Output directory containing the parser files.")
	return p


def print_job_banner(pdb_id: str) -> None:
	print("------------------------------------------------------")
	print(f"------------------------ {pdb_id} ------------------------")
	print("------------------------------------------------------")


def parse_order_id(raw_value: str) -> int:
	try:
		order_id = int(raw_value.strip())
	except ValueError as exc:
		raise ValueError("Invalid order_id. Supported values are 1 and 9.") from exc

	if order_id not in SUPPORTED_ORDER_IDS:
		raise ValueError("Invalid order_id. Supported values are 1 and 9.")

	return order_id


def read_instance_reorder_params(path: Path) -> dict[str, str]:
	params = read_params(path)

	required_keys = ("model_number", "chain_id", "order_id")
	missing_keys = [key for key in required_keys if key not in params]
	if missing_keys:
		raise ValueError(
			f"Missing required reordering parameter(s) in {path}: {', '.join(missing_keys)}"
		)

	try:
		int(params["model_number"])
	except ValueError as exc:
		raise ValueError("model_number must be an integer in the reordering params file.") from exc

	chain_id = params["chain_id"].strip()
	if not chain_id:
		raise ValueError("chain_id cannot be empty in the reordering params file.")
	params["chain_id"] = chain_id

	parse_order_id(params["order_id"])

	return params


def build_ddgp_order_vector(
	reorder_params: dict[str, str],
	out_dir: str | Path,
	pdb_id: str,
	order_id: int,
) -> np.ndarray:
	chosen_model = int(reorder_params["model_number"])
	chosen_chain = str(reorder_params["chain_id"])
	xfile = Path(out_dir) / pdb_id / f"X_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"

	df_x = read_space_separated_file(xfile)
	res_ids = df_x.iloc[:, 2].unique().tolist()

	if not res_ids:
		raise ValueError(f"No residues were found in {xfile}")

	first_residue = int(res_ids[0])
	last_residue = int(res_ids[-1])
	nres = last_residue - first_residue + 1

	if nres <= 0:
		raise ValueError(f"Invalid residue interval found in {xfile}")

	return np.full(nres, order_id, dtype=int)


# -----------------------------------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
	arg_parser = build_parser()
	args = arg_parser.parse_args(argv)

	try:
		pdb_ids = read_pdb_ids(args.pdb_ids_file)
		reorder_params = read_instance_reorder_params(args.reorder_params_file)
		order_choice = parse_order_id(reorder_params["order_id"])
	except Exception as exc:
		print(f"[ERROR] {exc}", file=sys.stderr)
		return 2

	out_dir = Path(args.out_dir)
	if not out_dir.exists():
		print(f"[ERROR] Output directory does not exist: {out_dir}", file=sys.stderr)
		return 2

	for pdb_id in pdb_ids:
		print_job_banner(pdb_id)
		try:
			ddgp_order_vec = build_ddgp_order_vector(reorder_params, out_dir, pdb_id, order_choice)
			skip_flag = reorder_instance(reorder_params, out_dir, pdb_id, ddgp_order_vec)
			if skip_flag == 1:
				continue
		except Exception as exc:
			print(f"[{pdb_id}] Skipped: {exc}")
			continue

	return 0


# -----------------------------------------------------------------------------------------------------
if __name__ == "__main__":
	raise SystemExit(main())


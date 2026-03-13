"""Minimal runner for pdb-parser.

Usage:
	python -m pdb_parser.main <pdb_ids.dat> <params.cfg> <pdb_data_dir> <out_dir>

Arguments
---------
pdb_ids.dat	: text file with one PDB id per line (comments with '#')
params.cfg	: simple key=value file (comments with '#')
pdb_dat_dir	: pdb data directory to be created if missing
out_dir		: output directory to be created if missing

Notes
-----

"""

from __future__ import annotations
from pathlib import Path
import argparse
import sys
from typing import List, Dict
import numpy as np

from pdb_parser.utils import *
from pdb_parser.io import *
from pdb_parser.pipeline import *
from pdb_parser.reordering import *

from pdb_parser.geometry import *
# -----------------------------------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
	p = argparse.ArgumentParser(
		prog="pdb-parser",
		description="Read PDB id list and params, validate, and prepare output directory.",
	)
	p.add_argument("pdb_ids_file", type=Path, help="Path to a file with one PDB id per line.")
	p.add_argument("params_file", type=Path, help="Path to a simple key=value params file.")
	p.add_argument("pdb_data_dir", type=Path, help="PDB data directory to create if missing.")
	p.add_argument("out_dir", type=Path, help="Output directory to create if missing.")
	return p
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
	
	pdb_data_dir = Path(args.pdb_data_dir)
	out_dir = Path(args.out_dir)
	
	ensure_dir(pdb_data_dir)
	ensure_dir(out_dir)
		
	for pdb_id in pdb_ids:
		download_pdb(pdb_id, pdb_data_dir)
		try:
			parser(params, pdb_data_dir, out_dir, pdb_id, True)
						
			if params.get("atom_selection", "").strip().lower() == "backbone_plus_hydrogens":
				
				out_dir_i = Path(out_dir) / pdb_id
				chosen_model = int(params["model_number"])
				chosen_chain = str(params["chain_id"])
				Xfile = out_dir_i / f"X_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"
				df_X = read_space_separated_file(Xfile)
				n = df_X.shape[0]
				nres = int(df_X.iat[n - 1, 2])
				ddgp_order_vec = np.full(nres, 9, dtype=int)
					
				skip_flag = sort_instance(params, out_dir, pdb_id, ddgp_order_vec)
				
				if skip_flag == 1:
					continue
				
		except ValueError as e:
			print(f"[{pdb_id}] Skipped: {e}")
			continue
		
	return 0
# -----------------------------------------------------------------------------------------------------
if __name__ == "__main__":
	raise SystemExit(main())

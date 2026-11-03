"""Application-facing helpers for running parser and remapping jobs.

This module provides a small library surface that sits above the low-level
pipeline functions and below the CLI entry points. The goal is to make the
core workflows reusable without depending on command-line modules.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

import numpy as np

from pdb_parser.io import download_pdb, read_validated_positional_file
from pdb_parser.io.sequence import validate_residue_sequence
from pdb_parser.messages import print_warning
from pdb_parser.output_layout import (
	normalize_layout_components,
	parser_structure_file,
)
from pdb_parser.pipeline import prepare_remap_execution_context, remap_instance
from pdb_parser.public_config import DMDGP_ORDER_TYPE, parse_public_order_type
from pdb_parser.remap_cli_text import (
	PARSER_STRUCTURE_REQUIRED_MESSAGE,
)
from pdb_parser.remapping import detect_remap_atom_layout

OrderChoice = int | str
ParserConfigValue = str | int | float
ParserParams = dict[str, ParserConfigValue]
BACKBONE_REMAP_ORDER_IGNORED_WARNING = (
	"Remapping detected parser outputs with atom_selection='backbone'; "
	"the configured order_id and order_type are ignored for this layout "
	"(use order_type: dmdgp to request its fixed N CA C order explicitly)."
)
DMDGP_REQUIRES_BACKBONE_MESSAGE = (
	"order_type: dmdgp requires parser outputs with atom_selection: backbone (fixed N CA C "
	"order); these outputs use atom_selection: backbone_plus_hydrogens."
)


class ParserWorkflow(Protocol):
	"""Callable parser workflow used by the application helpers."""

	def __call__(
		self,
		params: ParserParams,
		pdb_data_dir: str | Path,
		seed_dir: str | Path,
		out_dir: str | Path,
		pdb_id: str,
		remove_tmp_dir: bool = False,
	) -> None: ...


def _require_existing_parser_stage_file(path: Path, description: str) -> Path:
	"""Require one parser-stage input file with a stage-specific error message."""
	if not path.exists():
		raise ValueError(f"{description} was not found: {path}")
	return path


def _validate_residue_interval(res_ids: list[int], xfile: Path) -> int:
	"""Validate the current residue-interval convention and return its size."""
	return len(validate_residue_sequence(res_ids, context=f"remapping input {xfile}"))


def _load_remap_structure_residue_ids(
	out_dir: str | Path,
	pdb_id: str,
	chosen_model: int,
	chosen_chain: str,
) -> list[int]:
	"""Load residue ids without requiring a complete remappable atom layout."""
	xfile = parser_structure_file(out_dir, pdb_id, chosen_model, chosen_chain)
	_require_existing_parser_stage_file(xfile, PARSER_STRUCTURE_REQUIRED_MESSAGE)
	df_x = read_validated_positional_file(xfile, "X")
	# X is the authoritative source for the residue sequence.
	res_ids = list(dict.fromkeys(df_x.iloc[:, 2].astype(int).tolist()))
	if not res_ids:
		raise ValueError(f"No residues were found in {xfile}")
	_validate_residue_interval([int(residue_id) for residue_id in res_ids], xfile)
	return [int(residue_id) for residue_id in res_ids]


def _detect_remap_input_layout(
	out_dir: str | Path,
	pdb_id: str,
	chosen_model: int,
	chosen_chain: str,
) -> tuple[list[int], str, Path]:
	"""Load the parser structure file and detect the public remap layout."""
	xfile = parser_structure_file(out_dir, pdb_id, chosen_model, chosen_chain)
	_require_existing_parser_stage_file(
		xfile,
		PARSER_STRUCTURE_REQUIRED_MESSAGE,
	)
	df_x = read_validated_positional_file(xfile, "X")
	res_ids = list(dict.fromkeys(df_x.iloc[:, 2].astype(int).tolist()))

	if not res_ids:
		raise ValueError(f"No residues were found in {xfile}")

	residue_ids = [int(residue_id) for residue_id in res_ids]
	_validate_residue_interval(residue_ids, xfile)
	return residue_ids, detect_remap_atom_layout(df_x, residue_ids), xfile


def build_ddgp_order_vector(
	remap_params: dict[str, str],
	out_dir: str | Path,
	pdb_id: str,
	order_id: OrderChoice,
	remap_context: object | None = None,
) -> np.ndarray:
	"""Build the residue-wise DDGP order vector for one parsed instance."""
	chosen_model = int(remap_params["model_number"])
	chosen_chain = str(remap_params["chain_id"])
	xfile = parser_structure_file(out_dir, pdb_id, chosen_model, chosen_chain)
	if remap_context is None:
		res_ids = _load_remap_structure_residue_ids(out_dir, pdb_id, chosen_model, chosen_chain)
	else:
		res_ids = [int(residue_id) for residue_id in remap_context["res_ids"]]
	nres = _validate_residue_interval(res_ids, xfile)

	return np.full(nres, order_id, dtype=int)


def run_parser_job(
	params: ParserParams,
	pdb_data_dir: str | Path,
	seed_dir: str | Path,
	out_dir: str | Path,
	pdb_id: str,
	parser_fn: ParserWorkflow,
	download_pdb_fn=download_pdb,
	remove_tmp_dir: bool = True,
) -> None:
	"""Run one parser job using injected workflow functions."""
	pdb_id, _, _, _ = normalize_layout_components(pdb_id, 1, "A")
	download_pdb_fn(pdb_id, pdb_data_dir)
	parser_fn(params, pdb_data_dir, seed_dir, out_dir, pdb_id, remove_tmp_dir)


def run_remap_job(
	remap_params: dict[str, str],
	out_dir: str | Path,
	pdb_id: str,
	order_choice: OrderChoice,
	remap_instance_fn=remap_instance,
	prepare_remap_execution_context_fn=prepare_remap_execution_context,
	build_ddgp_order_vector_fn=build_ddgp_order_vector,
	remap_context=None,
) -> int:
	"""Run one remapping job using injected workflow functions.

	``remap_context`` reuses inputs already read and validated by
	``prepare_remap_execution_context`` (see :func:`run_remap_jobs`).
	"""
	chosen_model = int(remap_params["model_number"])
	chosen_chain = str(remap_params["chain_id"])
	pdb_id, _, _, _ = normalize_layout_components(
		pdb_id, chosen_model, chosen_chain, order_choice
	)
	if remap_context is not None:
		nres = int(remap_context["nres"])
		atom_layout = str(remap_context["atom_layout"])
	elif remap_instance_fn is remap_instance or build_ddgp_order_vector_fn is build_ddgp_order_vector:
		remap_context = prepare_remap_execution_context_fn(remap_params, out_dir, pdb_id)
		nres = int(remap_context["nres"])
		atom_layout = str(remap_context["atom_layout"])
	else:
		res_ids, atom_layout, xfile = _detect_remap_input_layout(
			out_dir,
			pdb_id,
			chosen_model,
			chosen_chain,
		)
		nres = _validate_residue_interval(res_ids, xfile)

	order_type = parse_public_order_type(remap_params)
	if order_type == DMDGP_ORDER_TYPE and atom_layout != "backbone":
		raise ValueError(DMDGP_REQUIRES_BACKBONE_MESSAGE)
	if atom_layout == "backbone":
		if order_type != DMDGP_ORDER_TYPE:
			print_warning(BACKBONE_REMAP_ORDER_IGNORED_WARNING)
		ddgp_order_vec = np.zeros(nres, dtype=int)
	elif build_ddgp_order_vector_fn is build_ddgp_order_vector:
		ddgp_order_vec = build_ddgp_order_vector_fn(
			remap_params,
			out_dir,
			pdb_id,
			order_choice,
			remap_context=remap_context,
		)
	else:
		ddgp_order_vec = build_ddgp_order_vector_fn(remap_params, out_dir, pdb_id, order_choice)

	if remap_instance_fn is remap_instance:
		return remap_instance_fn(
			remap_params,
			out_dir,
			pdb_id,
			ddgp_order_vec,
			remap_context=remap_context,
		)
	return remap_instance_fn(remap_params, out_dir, pdb_id, ddgp_order_vec)


REMAP_SKIPPED_MESSAGE = "structure skipped by the remap stage; no remapped files were written."


def _single_order_params(remap_params: dict[str, str], order_type: str, order_id: int) -> dict[str, str]:
	"""Return the remap parameters of one ``(order_type, order_id)`` combination."""
	return {**remap_params, "order_type": order_type, "order_id": str(order_id)}


def run_remap_jobs(
	remap_params: dict[str, str],
	out_dir: str | Path,
	pdb_id: str,
	order_choices: list[int],
	order_types: list[str],
) -> None:
	"""Remap one structure for every ``(order_type, order_id)``, reading its inputs once.

	The X/A/I files are read and validated a single time and every combination is written
	from them. A combination that fails is reported with a warning and the others still run;
	the structure then fails with a summary. ``backbone`` outputs ignore the order, so they
	are written once. A single combination behaves exactly as :func:`run_remap_job`.
	"""
	combinations = [(order_type, order_id) for order_type in order_types for order_id in order_choices]
	if len(combinations) == 1:
		order_type, order_id = combinations[0]
		if run_remap_job(_single_order_params(remap_params, order_type, order_id), out_dir, pdb_id, order_id) == 1:
			raise ValueError(REMAP_SKIPPED_MESSAGE)
		return

	first_params = _single_order_params(remap_params, *combinations[0])
	pdb_id, _, _, _ = normalize_layout_components(
		pdb_id, int(remap_params["model_number"]), str(remap_params["chain_id"]), combinations[0][1]
	)
	shared_context = prepare_remap_execution_context(first_params, out_dir, pdb_id)
	# Order-independent results (residue checks, distance lookups) computed by the first order.
	shared_context["shared_cache"] = {}
	if shared_context["atom_layout"] == "backbone":
		combinations = combinations[:1]

	failures: list[str] = []
	for order_type, order_id in combinations:
		context = {**shared_context, "order_type": order_type, "chosen_order_id": order_id}
		params = _single_order_params(remap_params, order_type, order_id)
		try:
			skip_flag = run_remap_job(params, out_dir, pdb_id, order_id, remap_context=context)
		except Exception as exc:
			print_warning(f"{pdb_id}: order_type {order_type}, order_id {order_id} failed: {exc}")
			failures.append(f"{order_type} {order_id}")
			continue
		if skip_flag == 1:
			# Missing atoms belong to the structure, not to the order: every order would skip.
			raise ValueError(REMAP_SKIPPED_MESSAGE)
	if failures:
		raise ValueError(
			f"{len(failures)} of {len(combinations)} remap orders failed ({', '.join(failures)}); "
			"the other orders were written."
		)

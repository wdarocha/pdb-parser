"""Helpers for the public on-disk layout of generated files."""

from __future__ import annotations

import re
from pathlib import Path

from pdb_parser.public_config import DDGP_ORDER_TYPE, DMDGP_ORDER_TYPE, REDMDGP_ORDER_TYPE

_PDB_ID_RE = re.compile(r"^[A-Z0-9]{4}$")


def normalize_layout_components(
	pdb_id: str,
	model_number: int,
	chain_id: str,
	order_id: int | str | None = None,
) -> tuple[str, int, str, int | str | None]:
	"""Normalize and validate user-controlled components of output names."""
	pid = str(pdb_id).strip().upper()
	if not _PDB_ID_RE.fullmatch(pid):
		raise ValueError("Invalid PDB id for output layout; expected four alphanumeric characters.")
	try:
		model = int(model_number)
	except (TypeError, ValueError) as exc:
		raise ValueError("model_number must be a positive integer for output layout.") from exc
	if model < 1:
		raise ValueError("model_number must be a positive integer for output layout.")
	chain = str(chain_id).strip()
	if not chain or chain in {".", ".."} or any(c in chain for c in "/\\") or any(ord(c) < 32 for c in chain):
		raise ValueError("chain_id contains an invalid path component.")
	if order_id is None:
		order = None
	elif isinstance(order_id, str):
		order = order_id.strip()
		if order.isdigit() and int(order) > 0:
			order = int(order)
		else:
			raise ValueError("order_id contains an invalid output-name component.")
	else:
		if isinstance(order_id, bool):
			raise ValueError("order_id must be a positive integer.")
		order = int(order_id)
		if order < 1:
			raise ValueError("order_id must be positive.")
	return pid, model, chain, order


def pdb_output_dir(out_dir: str | Path, pdb_id: str) -> Path:
	"""Return the per-PDB output directory."""
	pid, _, _, _ = normalize_layout_components(pdb_id, 1, "A")
	return Path(out_dir) / pid


def parser_file_suffix(pdb_id: str, model_number: int, chain_id: str) -> str:
	"""Return the shared parser-stage filename suffix."""
	pid, model, chain, _ = normalize_layout_components(pdb_id, model_number, chain_id)
	return f"{pid}_model{model}_chain{chain}.dat"


def parser_structure_file(out_dir: str | Path, pdb_id: str, model_number: int, chain_id: str) -> Path:
	"""Return the parser-stage `X_*` file path."""
	return pdb_output_dir(out_dir, pdb_id) / f"X_{parser_file_suffix(pdb_id, model_number, chain_id)}"


def filtered_atoms_file(out_dir: str | Path, pdb_id: str, model_number: int, chain_id: str) -> Path:
	"""Return the persistent filtered-coordinate TSV path."""
	return pdb_output_dir(out_dir, pdb_id) / f"filtered_atoms_{parser_file_suffix(pdb_id, model_number, chain_id)}"


def parser_metadata_file(out_dir: str | Path, pdb_id: str, model_number: int, chain_id: str) -> Path:
	"""Return the parser-stage metadata sidecar path."""
	pid, model, chain, _ = normalize_layout_components(pdb_id, model_number, chain_id)
	return pdb_output_dir(out_dir, pid) / f"metadata_{pid}_model{model}_chain{chain}.json"


def parser_angular_file(out_dir: str | Path, pdb_id: str, model_number: int, chain_id: str) -> Path:
	"""Return the parser-stage `A_*` file path."""
	return pdb_output_dir(out_dir, pdb_id) / f"A_{parser_file_suffix(pdb_id, model_number, chain_id)}"


def parser_distance_file(out_dir: str | Path, pdb_id: str, model_number: int, chain_id: str) -> Path:
	"""Return the parser-stage `I_*` file path."""
	return pdb_output_dir(out_dir, pdb_id) / f"I_{parser_file_suffix(pdb_id, model_number, chain_id)}"


def remapped_output_dir(out_dir: str | Path, pdb_id: str) -> Path:
	"""Return the remapped output directory for one PDB identifier."""
	return pdb_output_dir(out_dir, pdb_id) / "remapped"


def remapped_file_suffix(
	pdb_id: str,
	model_number: int,
	chain_id: str,
	order_id: int | str,
	atom_layout: str = "backbone_plus_hydrogens",
	order_type: str = DDGP_ORDER_TYPE,
) -> str:
	"""Return the shared remapped-stage filename suffix."""
	pdb_id, model_number, chain_id, order_id = normalize_layout_components(
		pdb_id, model_number, chain_id, order_id
	)
	if atom_layout not in {"backbone", "backbone_plus_hydrogens"}:
		raise ValueError("Unsupported atom layout for remapped output.")
	if atom_layout == "backbone":
		# Backbone outputs always use the fixed DMDGP order; order_id and order_type are ignored.
		return f"{pdb_id}_model{model_number}_chain{chain_id}_dmdgp_order.dat"
	if order_type == DMDGP_ORDER_TYPE:
		raise ValueError("order_type: dmdgp applies only to backbone outputs.")
	if order_type == REDMDGP_ORDER_TYPE:
		return f"{pdb_id}_model{model_number}_chain{chain_id}_redmdgp_hc_order_{order_id}.dat"
	return f"{pdb_id}_model{model_number}_chain{chain_id}_ddgp_order_{order_id}.dat"


def remapped_structure_file(
	out_dir: str | Path,
	pdb_id: str,
	model_number: int,
	chain_id: str,
	order_id: int | str,
	atom_layout: str = "backbone_plus_hydrogens",
	order_type: str = DDGP_ORDER_TYPE,
) -> Path:
	"""Return the remapped-stage `X_*` file path."""
	return remapped_output_dir(out_dir, pdb_id) / f"X_{remapped_file_suffix(pdb_id, model_number, chain_id, order_id, atom_layout, order_type)}"


def remapped_distance_file(
	out_dir: str | Path,
	pdb_id: str,
	model_number: int,
	chain_id: str,
	order_id: int | str,
	atom_layout: str = "backbone_plus_hydrogens",
	order_type: str = DDGP_ORDER_TYPE,
) -> Path:
	"""Return the remapped-stage `I_*` file path."""
	return remapped_output_dir(out_dir, pdb_id) / f"I_{remapped_file_suffix(pdb_id, model_number, chain_id, order_id, atom_layout, order_type)}"


def remapped_clique_file(
	out_dir: str | Path,
	pdb_id: str,
	model_number: int,
	chain_id: str,
	order_id: int | str,
	atom_layout: str = "backbone_plus_hydrogens",
	order_type: str = DDGP_ORDER_TYPE,
) -> Path:
	"""Return the remapped-stage `T_*` file path."""
	return remapped_output_dir(out_dir, pdb_id) / f"T_{remapped_file_suffix(pdb_id, model_number, chain_id, order_id, atom_layout, order_type)}"

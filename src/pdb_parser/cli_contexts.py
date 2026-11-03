"""Shared context-loading helpers for the public CLI entry points."""

from __future__ import annotations

from pathlib import Path

from pdb_parser.pipeline import validate_public_parser_params
from pdb_parser.public_config import (
	parse_public_order_choice,
	parse_public_order_types,
	parse_public_remap_order_choices,
	validate_public_remap_params,
)
from pdb_parser.utils import ensure_dir, ensure_existing_dir, read_params, read_pdb_ids


def load_public_pdb_ids(path: Path) -> list[str]:
	"""Load and validate the public PDB-id list input."""
	return read_pdb_ids(path)


def load_public_parser_params(path: Path) -> dict[str, str]:
	"""Load and validate the public parser configuration file."""
	params = read_params(path)
	validate_public_parser_params(params)
	return params


def load_public_remap_params(path: Path) -> dict[str, str]:
	"""Load and validate the public remap configuration file."""
	params = read_params(path)
	validate_public_remap_params(params, path)
	return params


def ensure_public_parser_dirs(
	pdb_data_dir: Path,
	seed_dir: Path,
	out_dir: Path,
) -> tuple[Path, Path, Path]:
	"""Create and return the directories required by the public parser CLI."""
	return (
		ensure_dir(pdb_data_dir),
		ensure_dir(seed_dir),
		ensure_dir(out_dir),
	)


def load_public_remap_order_choice(params: dict[str, str]) -> int | str:
	"""Return the validated public remap order choice."""
	return parse_public_order_choice(params["order_id"])


def load_public_remap_order_choices(params: dict[str, str]) -> list[int]:
	"""Return the validated public remap order ids (one value, a list or ranges)."""
	return parse_public_remap_order_choices(params)


def load_public_remap_order_types(params: dict[str, str]) -> list[str]:
	"""Return the validated public remap order types."""
	return parse_public_order_types(params)


def ensure_public_remap_out_dir(out_dir: Path) -> Path:
	"""Require and return the public remap output directory."""
	return ensure_existing_dir(out_dir, "Output directory")

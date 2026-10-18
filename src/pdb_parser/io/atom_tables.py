"""Normalized tabular inputs shared by distance-constraint writers."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from pdb_parser.io.sequence import residue_sequence_from_ids

REQUIRED_FILTERED_ATOM_COLUMNS = ("atom_id", "atom_name", "resid", "resname", "x", "y", "z")


def load_filtered_atoms_table(tsv_structure_file: str | Path) -> pd.DataFrame:
	"""Load a filtered-atom TSV and normalize its column names."""
	path = Path(tsv_structure_file)
	if not path.exists():
		raise ValueError(f"Filtered input file not found: {path}")
	dataframe = pd.read_csv(path, sep="\t", dtype={"atom_name": str, "resname": str}, float_precision="round_trip")
	dataframe.columns = dataframe.columns.str.strip()
	_validate_filtered_atoms_table(dataframe, path)
	return dataframe


def _validate_filtered_atoms_table(dataframe: pd.DataFrame, path: Path) -> None:
	"""Validate the semantic contract shared by all filtered-table consumers."""
	missing = [column for column in REQUIRED_FILTERED_ATOM_COLUMNS if column not in dataframe.columns]
	if missing:
		raise ValueError(
			f"Invalid filtered atom table {path}: missing required columns "
			f"{', '.join(missing)}."
		)
	if dataframe.empty:
		return
	for column in ("atom_name", "resname"):
		if dataframe[column].isna().any() or dataframe[column].astype(str).str.strip().eq("").any():
			raise ValueError(f"Invalid filtered atom table {path}: '{column}' cannot be empty.")
	for column in ("atom_id", "resid", "x", "y", "z"):
		values = pd.to_numeric(dataframe[column], errors="coerce")
		if values.isna().any() or not values.map(float).map(math.isfinite).all():
			raise ValueError(f"Invalid filtered atom table {path}: '{column}' must be finite numeric data.")
	if any((pd.to_numeric(dataframe[column]) <= 0).any() for column in ("atom_id", "resid")):
		raise ValueError(f"Invalid filtered atom table {path}: atom IDs and residue IDs must be positive.")
	for column in ("atom_id", "resid"):
		values = pd.to_numeric(dataframe[column])
		if not values.map(lambda value: float(value).is_integer()).all():
			raise ValueError(f"Invalid filtered atom table {path}: '{column}' must contain integers.")
	if dataframe["atom_id"].duplicated().any():
		raise ValueError(f"Invalid filtered atom table {path}: atom IDs must be unique.")
	identity = dataframe["resid"].astype(int).astype(str) + "\x00" + dataframe["atom_name"].astype(str).str.strip().str.upper()
	if identity.duplicated().any():
		raise ValueError(f"Invalid filtered atom table {path}: duplicate atom for a residue identity.")
	residue_sequence_from_ids(dataframe["resid"], context=f"filtered atom table {path}")


def validate_filtered_atoms_table(
    dataframe: pd.DataFrame,
    source: str | Path = "<in-memory>",
    *,
    require_nonempty: bool = False,
) -> pd.DataFrame:
	"""Validate an already loaded filtered atom table."""
	_validate_filtered_atoms_table(dataframe, Path(source))
	if require_nonempty and dataframe.empty:
		raise ValueError(
			f"Invalid filtered atom table {source}: no residues were provided; "
			"a valid residue sequence is required."
		)
	return dataframe

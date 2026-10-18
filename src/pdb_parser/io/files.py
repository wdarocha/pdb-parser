from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pdb_parser.messages import print_ok

DistanceRowFrame = pd.DataFrame
StructureRowFrame = pd.DataFrame

# Distances (I) and angles (T) use 12 decimals: values below 1000 keep at most 15 significant
# digits, all of which a double carries faithfully (DBL_DIG); more decimals would only add
# rounding noise. Readers split on whitespace, so files written with other widths or
# precisions (e.g. %20.16f before 0.3.0) remain readable.
DISTANCE_ROW_FORMAT = "%5d %5d %6d %6d %16.12f %16.12f %4.4s %4.4s %s %s\n"
# Minimum widths align the columns; the explicit spaces keep them separable at any magnitude.
# Coordinates use the parser-stage X width (%9.3f covers the PDB range -999.999..9999.999).
COORDINATE_ROW_FORMAT = "%9.3f %9.3f %9.3f\n"
# T: four vertex ids, sign of tau (-1, 0, 1), |tau| and its half-width in degrees (<= 180),
# with the same precision as the distances in I.
CLIQUE_MATRIX_FORMAT = "%5d %5d %5d %5d %2d %16.12f %16.12f"
POSITIONAL_SCHEMAS = {
	"X": ("atom_id", "atom_name", "resid", "resname", "x", "y", "z"),
	"I": ("atom_id_i", "atom_id_j", "resid_i", "resid_j", "d_l", "d_u", "atom_name_i", "atom_name_j", "resname_i", "resname_j"),
	"A": ("resid", "resname", "omega_center", "omega_radius", "phi_center", "phi_radius", "psi_center", "psi_radius"),
}
POSITIONAL_NUMERIC_COLUMNS = {
	"X": (0, 2, 4, 5, 6),
	"I": (0, 1, 2, 3, 4, 5),
	"A": (0, 2, 3, 4, 5, 6, 7),
}
POSITIONAL_TEXT_COLUMNS = {
	"X": (1, 3),
	"I": (6, 7, 8, 9),
	"A": (1,),
}
DISTANCE_TABLE_COLUMNS = [
	"atom_id_i", "atom_id_j", "resid_i", "resid_j", "d_l", "d_u",
	"atom_name_i", "atom_name_j", "resname_i", "resname_j",
]


def distance_table(records=None, **columns) -> pd.DataFrame:
	"""Return an intermediate distance table with the columns of ``distance_constraints_*.dat``.

	Pass either ``records`` (one tuple per row, in column order) or one sequence per column.
	"""
	if columns:
		return pd.DataFrame({name: columns[name] for name in DISTANCE_TABLE_COLUMNS})
	return pd.DataFrame(list(records or []), columns=DISTANCE_TABLE_COLUMNS)


def _normalized_text(values: pd.Series) -> pd.Series:
	"""Return ``str(value).strip().upper()`` for each element, computed once per distinct value."""
	return values.map({value: str(value).strip().upper() for value in pd.unique(values)})


def _resolve_output_path(filepath: str | Path) -> Path:
	"""Return the normalized output path and ensure its parent directory exists."""
	output_path = Path(filepath)
	output_path.parent.mkdir(parents=True, exist_ok=True)
	return output_path


def read_space_separated_file(input_file: str | Path) -> pd.DataFrame:
	"""Read a space-separated file without header into a DataFrame."""
	input_path = Path(input_file)

	if not input_path.exists():
		raise ValueError(f"Input file not found: {input_path}")

	# round_trip parses each value to the nearest double; the default parser can be off by one
	# ulp, so even a value written with every digit (e.g. %20.16f) would not come back exactly.
	return pd.read_csv(
		input_path,
		sep=r"\s+",
		header=None,
		engine="c",
		float_precision="round_trip",
	)


def read_validated_positional_file(
	input_file: str | Path,
	schema: str,
) -> pd.DataFrame:
	"""Read one X/I/A file and validate it before positional consumers run."""
	input_path = Path(input_file)
	try:
		df = read_space_separated_file(input_path)
	except pd.errors.EmptyDataError as exc:
		expected = POSITIONAL_SCHEMAS.get(schema)
		format_hint = ", ".join(expected) if expected else "a known positional schema"
		raise ValueError(
			f"Invalid {schema} schema in {input_path}: file is empty; expected format is "
			f"{format_hint}."
		) from exc
	return validate_positional_schema(df, schema, input_path)


def validate_positional_schema(
	df: pd.DataFrame, schema: str, input_file: str | Path
) -> pd.DataFrame:
	"""Validate a fixed-column X/I/A table before positional consumers use it."""
	try:
		expected = POSITIONAL_SCHEMAS[schema]
	except KeyError as exc:
		raise ValueError(f"Unknown positional schema {schema!r}.") from exc
	actual_count = len(df.columns)
	if actual_count != len(expected):
		raise ValueError(
			f"Invalid {schema} schema in {input_file}: expected {len(expected)} columns "
			f"({', '.join(expected)}), found {actual_count}."
		)
	for column_index in POSITIONAL_NUMERIC_COLUMNS[schema]:
		column = expected[column_index]
		try:
			pd.to_numeric(df.iloc[:, column_index], errors="raise")
		except (TypeError, ValueError) as exc:
				raise ValueError(
				f"Invalid {schema} schema in {input_file}: column {column_index + 1} "
				f"'{column}' must be numeric; expected format is "
				f"{', '.join(expected)}."
				) from exc
	for column_index in POSITIONAL_TEXT_COLUMNS[schema]:
		if pd.api.types.infer_dtype(df.iloc[:, column_index], skipna=True) not in {"string", "empty"}:
			column = expected[column_index]
			raise ValueError(
				f"Invalid {schema} schema in {input_file}: column {column_index + 1} "
				f"'{column}' must contain text; expected format is {', '.join(expected)}."
			)
	_validate_semantics(df, schema, input_file)
	return df


def validate_positional_cross_references(
	df_X: pd.DataFrame,
	df_A: pd.DataFrame,
	df_I: pd.DataFrame,
	input_file: str | Path = "<positional tables>",
) -> None:
	"""Validate the relationships between one parser X/A/I table set.

	The positional schemas validate each file in isolation.  This second layer
	ensures that metadata cannot drift between files and that every distance
	constraint endpoint has a corresponding atom in ``X``.
	"""
	validate_positional_schema(df_X, "X", input_file)
	validate_positional_schema(df_A, "A", input_file)
	validate_positional_schema(df_I, "I", input_file)

	x_by_id: dict[int, tuple[int, str, str]] = {}
	for row in df_X.itertuples(index=False, name=None):
		atom_id = int(row[0])
		identity = (int(row[2]), str(row[1]).strip().upper(), str(row[3]).strip().upper())
		if atom_id in x_by_id:
			raise ValueError(f"Invalid positional tables in {input_file}: X atom_id values must be unique.")
		x_by_id[atom_id] = identity

	x_residue_order = list(dict.fromkeys(int(value) for value in df_X.iloc[:, 2]))
	a_residues = [int(value) for value in df_A.iloc[:, 0]]
	if a_residues != x_residue_order:
		raise ValueError(
			f"Invalid positional tables in {input_file}: A must contain exactly one row "
			"for each X residue, in order."
		)
	for row in df_A.itertuples(index=False, name=None):
		resid = int(row[0])
		resname = str(row[1]).strip().upper()
		if any(identity[0] == resid and identity[2] != resname for identity in x_by_id.values()):
			raise ValueError(
				f"Invalid positional tables in {input_file}: A residue metadata is inconsistent with X for resid {resid}."
			)

	# Each I endpoint (atom_id, resid, atom_name, resname) must match its X atom.
	x_identity = pd.DataFrame.from_dict(
		x_by_id, orient="index", columns=["resid", "atom_name", "resname"]
	)
	endpoint_ids, endpoint_bad = [], []
	for id_column, resid_column, name_column, resname_column in ((0, 2, 6, 8), (1, 3, 7, 9)):
		atom_ids = df_I.iloc[:, id_column].astype(int).to_numpy()
		expected = x_identity.reindex(atom_ids)
		endpoint_ids.append(atom_ids)
		endpoint_bad.append(
			expected["resid"].isna().to_numpy()
			| (expected["resid"].to_numpy() != df_I.iloc[:, resid_column].astype(int).to_numpy())
			| (expected["atom_name"].to_numpy() != _normalized_text(df_I.iloc[:, name_column]).to_numpy())
			| (expected["resname"].to_numpy() != _normalized_text(df_I.iloc[:, resname_column]).to_numpy())
		)
	bad_rows = np.flatnonzero(endpoint_bad[0] | endpoint_bad[1])
	if bad_rows.size:
		row = bad_rows[0]
		atom_id = endpoint_ids[0][row] if endpoint_bad[0][row] else endpoint_ids[1][row]
		raise ValueError(
			f"Invalid positional tables in {input_file}: I references atom_id {atom_id} "
			"absent from X or with inconsistent residue/atom metadata."
		)


def _validate_semantics(df: pd.DataFrame, schema: str, input_file: str | Path) -> None:
	"""Validate physical and relational invariants of a positional table."""
	def numeric(index: int) -> pd.Series:
		return pd.to_numeric(df.iloc[:, index], errors="raise")

	def require_finite(indexes: tuple[int, ...]) -> None:
		if not all(np.isfinite(numeric(index)).all() for index in indexes):
			raise ValueError(f"Invalid {schema} data in {input_file}: numeric values must be finite.")

	def require_positive(indexes: tuple[int, ...]) -> None:
		if any((numeric(index) <= 0).any() for index in indexes):
			raise ValueError(f"Invalid {schema} data in {input_file}: identifiers and residue IDs must be positive.")

	def require_positive_integers(indexes: tuple[int, ...]) -> None:
		for index in indexes:
			values = numeric(index)
			if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
				column = POSITIONAL_SCHEMAS[schema][index]
				raise ValueError(
					f"Invalid {schema} data in {input_file}: column '{column}' must contain positive integers."
				)
		require_positive(indexes)

	def require_consistent_mapping(
		key_indexes: tuple[int, ...], value_indexes: tuple[int, ...]
	) -> None:
		mapping: dict[object, tuple[str, ...]] = {}
		for row in df.itertuples(index=False, name=None):
			key = tuple(row[index] for index in key_indexes)
			value = tuple(str(row[index]).strip().upper() for index in value_indexes)
			previous = mapping.setdefault(key, value)
			if previous != value:
				raise ValueError(
					f"Invalid {schema} data in {input_file}: repeated identifiers must keep consistent residue/atom metadata."
				)

	for index in range(len(df.columns)):
		if index not in POSITIONAL_NUMERIC_COLUMNS[schema]:
			column = df.iloc[:, index]
			if column.isna().any() or any(str(value).strip() == "" for value in pd.unique(column)):
				raise ValueError(f"Invalid {schema} data in {input_file}: text columns cannot be empty.")

	if schema == "X":
		require_finite((4, 5, 6))
		require_positive_integers((0, 2))
		if numeric(0).duplicated().any():
			raise ValueError(
				f"Invalid X data in {input_file}: atom_id values must be unique."
			)
		atom_identity = [
			(int(row[2]), str(row[1]).strip().upper())
			for row in df.itertuples(index=False, name=None)
		]
		if len(set(atom_identity)) != len(atom_identity):
			raise ValueError(
				f"Invalid X data in {input_file}: duplicate atom for a residue identity."
			)
		require_consistent_mapping((0,), (1, 3))
		require_consistent_mapping((2,), (3,))
	elif schema == "I":
		require_finite((4, 5))
		require_positive_integers((0, 1, 2, 3))
		if (numeric(4) < 0).any() or (numeric(5) < 0).any():
			raise ValueError(
				f"Invalid I data in {input_file}: distance bounds must be non-negative."
			)
		if (numeric(0) == numeric(1)).any():
			raise ValueError(
				f"Invalid I data in {input_file}: constraint endpoints must use distinct atoms."
			)
		if (numeric(4) > numeric(5)).any():
			raise ValueError(f"Invalid I data in {input_file}: every lower distance bound must be <= its upper bound.")
		same_residue = numeric(2) == numeric(3)
		if (_normalized_text(df.loc[same_residue, 8]) != _normalized_text(df.loc[same_residue, 9])).any():
			raise ValueError(
				f"Invalid I data in {input_file}: endpoints in one residue must use the same residue name."
			)
	elif schema == "A":
		require_finite((2, 3, 4, 5, 6, 7))
		require_positive_integers((0,))
		if any((numeric(index) < 0).any() for index in (3, 5, 7)):
			raise ValueError(f"Invalid A data in {input_file}: angular radii must be non-negative.")
		require_consistent_mapping((0,), (1,))


def iter_distance_rows(df_I: DistanceRowFrame):
	"""Yield the rows of a positional I table formatted with ``DISTANCE_ROW_FORMAT``.

	Each column is converted once (ids to ``int``, bounds to ``float``, names to ``str``)
	instead of row by row.
	"""
	columns = [df_I.iloc[:, index] for index in range(10)]
	values = zip(
		*(column.astype(int).tolist() for column in columns[:4]),
		*(column.astype(float).tolist() for column in columns[4:6]),
		*(column.astype(str).tolist() for column in columns[6:]),
	)
	for row in values:
		yield DISTANCE_ROW_FORMAT % row


def save_distances_from_df_structure(
	df_I: DistanceRowFrame,
	filepath: str | Path,
) -> None:
	"""Save one distance-constraint DataFrame using the current fixed-width layout."""
	output_path = _resolve_output_path(filepath)
	with output_path.open("w", encoding="utf-8") as output_stream:
		output_stream.writelines(iter_distance_rows(df_I))

	print_ok(f"Distance constraint file saved to: {output_path}")


def save_coordinates_from_df_structure(
	df_X: StructureRowFrame,
	filepath: str | Path,
) -> None:
	"""Save the trailing `x y z` columns of one structure DataFrame."""
	output_path = _resolve_output_path(filepath)
	with output_path.open("w", encoding="utf-8") as output_stream:
		for row in df_X.itertuples(index=False, name=None):
			x, y, z = row[-3:]
			output_stream.write(COORDINATE_ROW_FORMAT % (float(x), float(y), float(z)))

	print_ok(f"Coordinate file saved to: {output_path}")


def save_cliques_from_matrix_T(
	T: np.ndarray,
	filepath: str | Path,
) -> None:
	"""Save one clique matrix `T` using the current seven-column numeric format."""
	output_path = _resolve_output_path(filepath)
	try:
		matrix = np.asarray(T, dtype=float)
	except (TypeError, ValueError) as exc:
		raise ValueError("Invalid T matrix: all values must be numeric.") from exc
	if matrix.ndim != 2 or matrix.shape[1] != 7:
		raise ValueError(f"Invalid T matrix: expected a two-dimensional matrix with 7 columns, got shape {matrix.shape}.")
	if not np.isfinite(matrix).all():
		raise ValueError("Invalid T matrix: all values must be finite.")
	if not np.equal(matrix[:, :5], np.floor(matrix[:, :5])).all():
		raise ValueError("Invalid T matrix: the first five columns must contain integers.")
	np.savetxt(output_path, matrix, fmt=CLIQUE_MATRIX_FORMAT)
	print_ok(f"Cliques file saved to: {output_path}")

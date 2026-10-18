"""Merge and tighten distance-constraint files."""

from __future__ import annotations

from pathlib import Path
from typing import TypedDict

import numpy as np
import pandas as pd

from pdb_parser.io.files import (
	DISTANCE_TABLE_COLUMNS,
	iter_distance_rows,
	read_validated_positional_file,
	save_distances_from_df_structure,
)
from pdb_parser.messages import print_ok


class MergedDistanceRecord(TypedDict):
	atom_id_i: int
	atom_id_j: int
	resid_i: int
	resid_j: int
	atom_name_i: str
	atom_name_j: str
	resname_i: str
	resname_j: str
	d_l: float
	d_u: float


def merge_distance_constraint_files(
	files_dir: str | Path,
	instance_file: str | Path,
	*,
	tighten_upper_bounds: bool = False,
) -> None:
	"""Merge the ``distance_constraints_*.dat`` files of ``files_dir`` into ``instance_file``.

	See :func:`merge_distance_constraint_tables`; the parser passes its tables in memory.
	"""
	files_dir = Path(files_dir)
	input_files = sorted(files_dir.glob("distance_constraints_*.dat"))
	if not input_files:
		raise ValueError(f"No distance_constraints_*.dat files found in: {files_dir}")
	all_entries: list[pd.DataFrame] = []
	for input_file in input_files:
		df = pd.read_csv(input_file, sep="\t", float_precision="round_trip")
		df.columns = df.columns.str.strip()
		missing_columns = [column for column in DISTANCE_TABLE_COLUMNS if column not in df.columns]
		if missing_columns:
			raise ValueError(f"Missing required columns in {input_file}: {', '.join(missing_columns)}")
		all_entries.append(df)
	merge_distance_constraint_tables(all_entries, instance_file, tighten_upper_bounds=tighten_upper_bounds)


def merge_distance_constraint_tables(
	tables: list[pd.DataFrame],
	instance_file: str | Path,
	*,
	tighten_upper_bounds: bool = False,
) -> None:
	"""Merge intermediate constraint tables, preserving precise values when available.

	The tables hold the full-precision values (the ``tmp/`` files are only a rounded
	record of them). With ``tighten_upper_bounds`` the shortest-path tightening runs on
	the merged values before they are written, so each bound is rounded once.
	"""
	if not tables:
		raise ValueError("No distance constraint tables to merge.")
	all_entries = []
	for table in tables:
		table = table[DISTANCE_TABLE_COLUMNS].copy()
		table["d_l"] = table["d_l"].astype(float)
		table["d_u"] = table["d_u"].astype(float)
		all_entries.append(table)

	group_keys = [
		"atom_id_i", "atom_id_j", "resid_i", "resid_j",
		"atom_name_i", "atom_name_j", "resname_i", "resname_j",
	]
	# Per pair: the smallest precise value wins; otherwise intersect the intervals.
	entries = pd.concat(all_entries, ignore_index=True)
	entries["fixed_value"] = entries["d_l"].where(np.isclose(entries["d_l"], entries["d_u"]))
	merged = (
		entries.groupby(group_keys, dropna=False)
		.agg(fixed_value=("fixed_value", "min"), d_l=("d_l", "max"), d_u=("d_u", "min"))
		.reset_index()
	)
	has_fixed = merged["fixed_value"].notna()
	merged.loc[has_fixed, "d_l"] = merged.loc[has_fixed, "fixed_value"]
	merged.loc[has_fixed, "d_u"] = merged.loc[has_fixed, "fixed_value"]
	inconsistent = merged[merged["d_l"] > merged["d_u"]]
	if not inconsistent.empty:
		pair_data = inconsistent.iloc[0]
		raise ValueError(
			"Inconsistent merged interval for pair "
			f"(atom_id_i={pair_data['atom_id_i']}, atom_id_j={pair_data['atom_id_j']}, "
			f"resid_i={pair_data['resid_i']}, resid_j={pair_data['resid_j']}, "
			f"atom_name_i={pair_data['atom_name_i']}, atom_name_j={pair_data['atom_name_j']}, "
			f"resname_i={pair_data['resname_i']}, resname_j={pair_data['resname_j']}): "
			f"d_l={pair_data['d_l']:.16f} > d_u={pair_data['d_u']:.16f}."
		)

	final_df = (
		merged[[*group_keys, "d_l", "d_u"]]
		.sort_values(by=["atom_id_i", "atom_id_j", "resid_i", "resid_j"])
		.reset_index(drop=True)
	)
	if tighten_upper_bounds:
		final_df["d_u"] = _shortest_path_upper_bounds(
			final_df["atom_id_i"].astype(int).tolist(),
			final_df["atom_id_j"].astype(int).tolist(),
			final_df["d_u"].astype(float).tolist(),
			instance_file,
		)
	positional_columns = [
		"atom_id_i", "atom_id_j", "resid_i", "resid_j", "d_l", "d_u",
		"atom_name_i", "atom_name_j", "resname_i", "resname_j",
	]
	with Path(instance_file).open("w", encoding="utf-8") as output_file:
		output_file.writelines(iter_distance_rows(final_df[positional_columns]))
	print_ok(f"Merged distance constraints saved to: {instance_file}")


def tightening_upper_bounds(instance_file: str | Path) -> None:
	"""Replace every upper bound of an I file with its shortest path over the finite bounds."""
	instance_file = Path(instance_file)
	df_instance = read_validated_positional_file(instance_file, "I")
	if df_instance.empty:
		save_distances_from_df_structure(df_instance, instance_file)
		return
	if df_instance.shape[1] < 6:
		raise ValueError(f"Expected at least 6 columns in distance file {instance_file}, found {df_instance.shape[1]}.")

	df_instance[5] = _shortest_path_upper_bounds(
		df_instance[0].astype(int).tolist(),
		df_instance[1].astype(int).tolist(),
		df_instance[5].astype(float).tolist(),
		instance_file,
	)
	save_distances_from_df_structure(df_instance, instance_file)


def _shortest_path_upper_bounds(
	atom_ids_i: list[int],
	atom_ids_j: list[int],
	upper_bounds: list[float],
	source: str | Path,
) -> list[float]:
	"""Return, for every pair, its shortest path in the graph of the finite (<= 998) upper bounds.

	Pairs whose own bound is not finite (the 999 van der Waals rows) get that value too; a pair
	without a finite path raises ``ValueError``.
	"""
	# Imported here so that importing the package (e.g. by the remap stage) does not load SciPy.
	from scipy.sparse import csr_matrix
	from scipy.sparse.csgraph import dijkstra as csgraph_dijkstra

	ids_i = np.asarray(atom_ids_i, dtype=np.int64)
	ids_j = np.asarray(atom_ids_j, dtype=np.int64)
	bounds = np.asarray(upper_bounds, dtype=float)
	if not (len(ids_i) == len(ids_j) == len(bounds)):
		raise ValueError("atom_ids_i, atom_ids_j and upper_bounds must have the same length.")
	vertex_order = np.unique(np.concatenate([ids_i, ids_j]))
	index_i = np.searchsorted(vertex_order, ids_i)
	index_j = np.searchsorted(vertex_order, ids_j)

	# One edge per unordered pair, weighted by its smallest finite upper bound.
	finite = bounds <= 998.0
	low = np.minimum(index_i, index_j)[finite]
	high = np.maximum(index_i, index_j)[finite]
	weights = bounds[finite]
	edge_key = low * len(vertex_order) + high
	by_key_then_weight = np.lexsort((weights, edge_key))
	sorted_keys = edge_key[by_key_then_weight]
	smallest = by_key_then_weight[np.r_[True, sorted_keys[1:] != sorted_keys[:-1]]] if len(sorted_keys) else by_key_then_weight
	low, high, weights = low[smallest], high[smallest], weights[smallest]

	graph = csr_matrix(
		(np.concatenate([weights, weights]), (np.concatenate([low, high]), np.concatenate([high, low]))),
		shape=(len(vertex_order), len(vertex_order)),
		dtype=float,
	)
	all_pairs_shortest_paths = csgraph_dijkstra(graph, directed=False, return_predecessors=False)
	tightened = all_pairs_shortest_paths[index_i, index_j]
	unreachable = np.flatnonzero(~np.isfinite(tightened))
	if unreachable.size:
		first = unreachable[0]
		raise ValueError(
			f"No finite shortest-path upper bound found for pair ({int(ids_i[first])}, {int(ids_j[first])}) in {source}."
		)
	return tightened.tolist()

"""Long-range distance constraints derived from NMR atom selections."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from pdb_parser.io.atom_tables import load_filtered_atoms_table, validate_filtered_atoms_table
from pdb_parser.io.distance_intervals import sample_centered_intervals
from pdb_parser.io.files import distance_table
from pdb_parser.io.sequence import residue_sequence_from_ids
from pdb_parser.io.vdw import lower_bounds_below
from pdb_parser.messages import print_ok
from pdb_parser.remapping.clique_patterns import ATOM_ALIAS_MAP


def _get_nmr_candidate_indices(atom_name: np.ndarray, atom_selection: str) -> np.ndarray:
	"""Return backbone proxies or explicit hydrogen atoms eligible for NMR restraints."""
	selection = (atom_selection or "").strip().lower()
	if selection == "backbone":
		return np.flatnonzero(np.isin(atom_name, ["N", "C"]))
	if selection in {"backbone_plus_hydrogens", "full_chain", "backbone_plus_neighbors"}:
		return np.flatnonzero(np.char.startswith(atom_name.astype(str), "H"))
	raise ValueError(f"Unsupported atom_selection: {selection!r}.")


def _discretization_hh_pairs(
	atom_name: np.ndarray, resid: np.ndarray, residue_ids: tuple[int, ...],
) -> set[frozenset[int]]:
	"""Return the H--H pairs used as interval discretization edges by the remap stage.

	The DDGP cliques take the interval torsion from ``HN_i--HA_i`` (pattern 6 and the
	first residue) and from ``HA_{i-1}--HN_i`` (pattern 10). ``HN`` and ``HA`` are
	resolved exactly as in the remap stage: the first available name of each alias list.
	"""
	names_by_residue: dict[int, dict[str, int]] = {}
	for index, (name, residue) in enumerate(zip(atom_name, resid)):
		names_by_residue.setdefault(int(residue), {})[str(name)] = index

	def resolve(residue: int, alias: str) -> int | None:
		names = names_by_residue.get(residue, {})
		return next((names[name] for name in ATOM_ALIAS_MAP[alias] if name in names), None)

	pairs: set[frozenset[int]] = set()
	for position, residue in enumerate(residue_ids):
		hn, ha = resolve(residue, "HN"), resolve(residue, "HA")
		if hn is not None and ha is not None:
			pairs.add(frozenset((hn, ha)))
		if position > 0 and hn is not None:
			previous_ha = resolve(residue_ids[position - 1], "HA")
			if previous_ha is not None:
				pairs.add(frozenset((previous_ha, hn)))
	return pairs


def nmr_distance_constraints(
	tsv_structure_file: str | Path, output_distance_file: str | Path,
	distance_constraints: str, atom_selection: str,
	epsilon_short: float, epsilon_long: float, max_distance: float,
	noe_strong: float, noe_medium: float, noe_weak: float,
	vdw_threshold: float, rng: np.random.Generator | None = None,
	*, hh_distance_percentage: float = 100.0,
	rng_hh_distance_subset: np.random.Generator | None = None,
):
	"""Write long-range NMR restraints using backbone proxies or explicit H atoms; return them.

	The H--H discretization edges used by the remap stage (``HN_i--HA_i`` and
	``HA_{i-1}--HN_i``) are always kept. Of the other restraints, only
	``hh_distance_percentage`` percent are kept, drawn uniformly without replacement by
	``rng_hh_distance_subset``. The draw runs after all intervals are built, so the kept
	restraints are identical to those of a 100% run.
	"""
	if not 0.0 <= hh_distance_percentage <= 100.0:
		raise ValueError("hh_distance_percentage must be between 0 and 100.")
	mode = (distance_constraints or "").strip().lower()
	if mode not in {"precise", "interval_centered", "interval_experimental"}:
		raise ValueError(f"Unsupported distance_constraints mode: {mode!r}.")
	df = load_filtered_atoms_table(tsv_structure_file)
	validate_filtered_atoms_table(df, tsv_structure_file, require_nonempty=True)
	header = "atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n"
	if df.empty:
		Path(output_distance_file).write_text(header, encoding="utf-8")
		print_ok(f"NMR distance table (mode={mode}) saved to (empty): {output_distance_file}")
		return distance_table()
	coords = df[["x", "y", "z"]].to_numpy(dtype=float)
	atom_id, resid = df["atom_id"].to_numpy(), df["resid"].to_numpy()
	residue_ids = tuple(residue_sequence_from_ids(resid, context="filtered atom table"))
	resname = df["resname"].str.strip().str.upper().to_numpy()
	atom_name = df["atom_name"].str.strip().str.upper().to_numpy()
	candidates = _get_nmr_candidate_indices(atom_name, atom_selection)
	rng = rng if rng is not None else np.random.default_rng()
	discretization_pairs = _discretization_hh_pairs(atom_name, resid, residue_ids)
	# Every candidate pair (i after j in candidate order), as in a double loop over candidates.
	pos_i, pos_j = np.tril_indices(len(candidates), -1)
	index_i, index_j = candidates[pos_i], candidates[pos_j]
	delta = coords[index_i] - coords[index_j]
	dist = np.sqrt(delta[:, 0] * delta[:, 0] + delta[:, 1] * delta[:, 1] + delta[:, 2] * delta[:, 2])
	d_l, d_u, kept = _nmr_intervals(
		mode, dist, resid[index_i], resid[index_j], epsilon_short, epsilon_long, max_distance,
		noe_strong, noe_medium, noe_weak, vdw_threshold, rng,
	)
	index_i, index_j, d_l, d_u = index_i[kept], index_j[kept], d_l[kept], d_u[kept]

	n_atoms = len(df)
	edge_keys = np.array(
		sorted(min(pair) * n_atoms + max(pair) for pair in map(tuple, discretization_pairs)), dtype=np.int64
	)
	protected = np.isin(np.minimum(index_i, index_j) * n_atoms + np.maximum(index_i, index_j), edge_keys).tolist()
	ids, resids = atom_id.astype(int), resid.astype(int)
	row_format = "%d\t%d\t%d\t%d\t%.12f\t%.12f\t%s\t%s\t%s\t%s\n"
	rows = [
		row_format % (int(ids[i]), int(ids[j]), int(resids[i]), int(resids[j]), lower, upper,
			atom_name[i], atom_name[j], resname[i], resname[j])
		for i, j, lower, upper in zip(index_i.tolist(), index_j.tolist(), d_l.tolist(), d_u.tolist())
	]
	n_generated = len(rows)
	kept_rows = _select_restraint_indices(len(rows), hh_distance_percentage, rng_hh_distance_subset, protected)
	rows = [rows[index] for index in kept_rows]
	index_i, index_j, d_l, d_u = index_i[kept_rows], index_j[kept_rows], d_l[kept_rows], d_u[kept_rows]
	with Path(output_distance_file).open("w", encoding="utf-8") as stream:
		stream.write(header)
		stream.writelines(rows)
	if len(rows) < n_generated:
		print_ok(
			f"percentage_hh_distance_constraints={hh_distance_percentage:g}: "
			f"kept {len(rows)} of {n_generated} NMR restraints "
			f"({sum(protected)} discretization edges always kept)."
		)
	print_ok(f"NMR-derived distances (mode={mode}) saved to: {output_distance_file}")
	return distance_table(
		atom_id_i=ids[index_i], atom_id_j=ids[index_j], resid_i=resids[index_i], resid_j=resids[index_j],
		d_l=d_l, d_u=d_u, atom_name_i=atom_name[index_i], atom_name_j=atom_name[index_j],
		resname_i=resname[index_i], resname_j=resname[index_j],
	)


def _nmr_intervals(
	mode: str, dist: np.ndarray, resid_i: np.ndarray, resid_j: np.ndarray,
	epsilon_short: float, epsilon_long: float, max_distance: float,
	noe_strong: float, noe_medium: float, noe_weak: float, vdw_threshold: float,
	rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
	"""Vectorized ``get_precise_interval`` / ``get_centered_interval`` / ``get_experimental_interval``.

	Returns ``(d_l, d_u, kept)``; a pair is kept when its scalar builder returns an interval.
	"""
	vdw = np.full(len(dist), float(vdw_threshold))
	if mode == "precise":
		return dist, dist, dist < max_distance
	if mode == "interval_experimental":
		upper = np.select(
			[dist < noe_strong, dist < noe_medium, dist < noe_weak], [noe_strong, noe_medium, noe_weak], np.nan
		)
		kept = ~np.isnan(upper)
		return lower_bounds_below(dist, vdw), upper, kept

	# interval_centered: residues are a validated consecutive sequence, so adjacency is |di| <= 1.
	epsilon = np.where(np.abs(resid_i.astype(np.int64) - resid_j.astype(np.int64)) <= 1, epsilon_short, epsilon_long)
	sampled = np.isfinite(dist) & np.isfinite(epsilon) & (epsilon > 0.0)
	if not (math.isfinite(vdw_threshold) and math.isfinite(max_distance) and vdw_threshold > 0.0 and max_distance > 0.0):
		sampled[:] = False
	sampled &= max_distance > dist
	d_l = np.zeros(len(dist))
	d_u = np.zeros(len(dist))
	kept = np.zeros(len(dist), dtype=bool)
	rows = np.flatnonzero(sampled)
	d_l[rows], d_u[rows], kept[rows] = sample_centered_intervals(
		dist[rows], epsilon[rows], lower_bounds_below(dist[rows], vdw[rows]), max_distance, rng
	)
	return d_l, d_u, kept


def _select_restraint_indices(
	count: int, percentage: float, rng: np.random.Generator | None,
	protected: list[bool] | None = None,
) -> list[int]:
	"""Return the kept row positions: every protected row plus round(percentage% of the others).

	The unprotected rows are chosen uniformly at random without replacement; the positions
	are returned in their original order.
	"""
	protected = protected if protected is not None else [False] * count
	free = [index for index, keep in enumerate(protected) if not keep]
	n_keep = int(round(len(free) * percentage / 100.0))
	if n_keep >= len(free):
		return list(range(count))
	if n_keep == 0:
		chosen: set[int] = set()
	else:
		if rng is None:
			raise ValueError("A random generator is required to select a subset of NMR restraints.")
		chosen = {free[index] for index in rng.choice(len(free), size=n_keep, replace=False)}
	return [index for index in range(count) if protected[index] or index in chosen]

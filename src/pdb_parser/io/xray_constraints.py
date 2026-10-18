"""Long-range CA--CA distance constraints for X-ray structures."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from pdb_parser.io.atom_tables import load_filtered_atoms_table, validate_filtered_atoms_table
from pdb_parser.io.distance_intervals import sample_centered_intervals
from pdb_parser.io.files import distance_table
from pdb_parser.io.sequence import residue_sequence_from_ids
from pdb_parser.io.vdw import VDW_CARBON_RADIUS, lower_bounds_below
from pdb_parser.messages import print_ok


def xray_ca_distance_constraints(
	tsv_structure_file: str | Path,
	output_distance_file: str | Path,
	distance_constraints: str,
	epsilon_long: float | None,
	max_distance: float,
	rng: np.random.Generator | None = None,
):
	"""Write non-adjacent long-range CA--CA constraints for X-ray structures and return them.

	Every pair of CA atoms (the later one first, in table order) is considered, except pairs
	in the same or adjacent residues. The intervals are those of ``get_precise_interval`` and
	``get_xray_ca_centered_interval``; ``interval_centered`` consumes the random stream as
	one ``get_xray_ca_centered_interval`` call per pair (see ``sample_centered_intervals``).
	"""
	mode = distance_constraints.strip().lower()
	if mode not in {"precise", "interval_centered"}:
		raise ValueError("X-ray CA--CA constraints require 'precise' or 'interval_centered'.")

	df = load_filtered_atoms_table(tsv_structure_file)
	validate_filtered_atoms_table(df, tsv_structure_file, require_nonempty=True)
	df = df[df["atom_name"].astype(str).str.strip().str.upper() == "CA"]
	# A validated consecutive sequence: adjacent residues differ by exactly 1.
	residue_sequence_from_ids(df["resid"], context="filtered atom table")

	atom_id = df["atom_id"].astype(int).to_numpy()
	resid = df["resid"].astype(int).to_numpy()
	resname = np.array([str(name).upper() for name in df["resname"]], dtype=object)
	coords = df[["x", "y", "z"]].to_numpy(dtype=float)

	index_i, index_j = np.tril_indices(len(df), -1)
	apart = np.abs(resid[index_i] - resid[index_j]) > 1
	index_i, index_j = index_i[apart], index_j[apart]
	delta = coords[index_i] - coords[index_j]
	dist = np.sqrt(delta[:, 0] * delta[:, 0] + delta[:, 1] * delta[:, 1] + delta[:, 2] * delta[:, 2])

	if mode == "precise":
		kept = dist < max_distance
		d_l = d_u = dist
	else:
		if epsilon_long is None:
			if len(dist):
				raise ValueError("epsilon_long is required for X-ray interval_centered constraints.")
			epsilon_long = 0.0
		d_l, d_u = np.zeros(len(dist)), np.zeros(len(dist))
		kept = np.zeros(len(dist), dtype=bool)
		if math.isfinite(epsilon_long) and math.isfinite(max_distance) and epsilon_long > 0.0 and max_distance > 0.0:
			rows = np.flatnonzero(np.isfinite(dist) & (max_distance > dist))
			rng = rng if rng is not None else np.random.default_rng()
			d_l[rows], d_u[rows], kept[rows] = sample_centered_intervals(
				dist[rows],
				np.full(len(rows), float(epsilon_long)),
				lower_bounds_below(dist[rows], np.full(len(rows), 2.0 * VDW_CARBON_RADIUS)),
				max_distance,
				rng,
			)

	index_i, index_j, d_l, d_u = index_i[kept], index_j[kept], d_l[kept], d_u[kept]
	row_format = "%d\t%d\t%d\t%d\t%.12f\t%.12f\tCA\tCA\t%s\t%s\n"
	with Path(output_distance_file).open("w", encoding="utf-8") as stream:
		stream.write("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n")
		stream.writelines(
			row_format % (atom_id[i], atom_id[j], resid[i], resid[j], lower, upper, resname[i], resname[j])
			for i, j, lower, upper in zip(index_i.tolist(), index_j.tolist(), d_l.tolist(), d_u.tolist())
		)
	print_ok(f"X-ray CA--CA distances (mode={mode}) saved to: {output_distance_file}")
	return distance_table(
		atom_id_i=atom_id[index_i], atom_id_j=atom_id[index_j], resid_i=resid[index_i], resid_j=resid[index_j],
		d_l=d_l, d_u=d_u, atom_name_i=np.full(len(index_i), "CA", dtype=object),
		atom_name_j=np.full(len(index_i), "CA", dtype=object),
		resname_i=resname[index_i], resname_j=resname[index_j],
	)

"""Van der Waals distance-constraint helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pdb_parser.io.atom_tables import load_filtered_atoms_table, validate_filtered_atoms_table
from pdb_parser.io.files import distance_table
from pdb_parser.io.sequence import are_adjacent_residues, residue_sequence_from_ids
from pdb_parser.io.topology import AtomKey, AtomPair, TopologyData
from pdb_parser.io.vdw import VDW_RADII, lower_bounds_below
from pdb_parser.messages import print_ok


def _is_omega_path(
	path_keys: tuple[AtomKey, AtomKey, AtomKey, AtomKey],
	residue_ids: tuple[int, ...],
) -> bool:
	"""Return True if the path matches backbone omega."""
	(resid_1, atom_1), (resid_2, atom_2), (resid_3, atom_3), (resid_4, atom_4) = path_keys

	return (
		atom_1 == "C"
		and atom_2 == "N"
		and atom_3 == "CA"
		and atom_4 == "C"
		and resid_2 == resid_3 == resid_4
		and are_adjacent_residues(residue_ids, resid_1, resid_2, context="topology residue sequence")
	)


def _is_phi_path(
	path_keys: tuple[AtomKey, AtomKey, AtomKey, AtomKey],
	residue_ids: tuple[int, ...],
) -> bool:
	"""Return True if the path matches backbone phi."""
	(resid_1, atom_1), (resid_2, atom_2), (resid_3, atom_3), (resid_4, atom_4) = path_keys

	return (
		atom_1 == "N"
		and atom_2 == "CA"
		and atom_3 == "C"
		and atom_4 == "N"
		and resid_1 == resid_2 == resid_3
		and are_adjacent_residues(residue_ids, resid_1, resid_4, context="topology residue sequence")
	)


def _is_psi_path(
	path_keys: tuple[AtomKey, AtomKey, AtomKey, AtomKey],
	residue_ids: tuple[int, ...],
) -> bool:
	"""Return True if the path matches backbone psi."""
	(resid_1, atom_1), (resid_2, atom_2), (resid_3, atom_3), (resid_4, atom_4) = path_keys

	return (
		atom_1 == "CA"
		and atom_2 == "C"
		and atom_3 == "N"
		and atom_4 == "CA"
		and resid_1 == resid_2
		and resid_3 == resid_4
		and are_adjacent_residues(residue_ids, resid_1, resid_3, context="topology residue sequence")
	)


def is_backbone_dihedral_path(
	path_keys: tuple[AtomKey, AtomKey, AtomKey, AtomKey],
	residue_ids: tuple[int, ...] | None = None,
) -> bool:
	"""Return True if the path matches backbone omega, phi, or psi."""
	if residue_ids is None:
		residue_ids = tuple(residue_sequence_from_ids((resid for resid, _ in path_keys), context="topology residue sequence"))
	reversed_path_keys = tuple(reversed(path_keys))
	return (
		_is_omega_path(path_keys, residue_ids)
		or _is_phi_path(path_keys, residue_ids)
		or _is_psi_path(path_keys, residue_ids)
		or _is_omega_path(reversed_path_keys, residue_ids)
		or _is_phi_path(reversed_path_keys, residue_ids)
		or _is_psi_path(reversed_path_keys, residue_ids)
	)


def build_vdw_torsion_pairs(
	topology: TopologyData,
) -> tuple[set[AtomPair], set[AtomPair]]:
	"""Build general relaxed 1-4 pairs and excluded backbone-dihedral pairs."""
	atom_index = topology["atom_index"]
	residue_ids = tuple(topology.get("residue_ids", residue_sequence_from_ids((resid for resid, _ in atom_index), context="topology residue sequence")))
	logical_graph = topology["logical_graph"]
	bonded_pairs = topology["bonded_pairs"]
	angle_pairs = topology["angle_pairs"]

	relaxed_torsion_pairs: set[AtomPair] = set()
	backbone_dihedral_pairs: set[AtomPair] = set()

	for start_key, first_neighbor_keys in logical_graph.items():
		start_index = atom_index.get(start_key)
		if start_index is None:
			continue

		for first_key in first_neighbor_keys:
			for second_key in logical_graph.get(first_key, set()):
				if second_key == start_key:
					continue

				for third_key in logical_graph.get(second_key, set()):
					if third_key in {start_key, first_key, second_key}:
						continue

					end_index = atom_index.get(third_key)
					if end_index is None or end_index == start_index:
						continue

					pair = (start_index, end_index) if start_index > end_index else (end_index, start_index)
					if pair in bonded_pairs or pair in angle_pairs:
						continue

					path_keys = (start_key, first_key, second_key, third_key)
					if is_backbone_dihedral_path(path_keys, residue_ids):
						backbone_dihedral_pairs.add(pair)
					else:
						relaxed_torsion_pairs.add(pair)

	return relaxed_torsion_pairs, backbone_dihedral_pairs


def _pairs_by_row(pairs: set[AtomPair]) -> dict[int, list[int]]:
	"""Group ``(i, j)`` index pairs (``i > j``) by their first index."""
	by_row: dict[int, list[int]] = {}
	for i, j in pairs:
		by_row.setdefault(i, []).append(j)
	return by_row


def vdw_distances(
	tsv_structure_file: str | Path,
	topology: TopologyData,
	output_distance_file: str | Path,
):
	"""Write van der Waals lower-bound constraints and return them as a table."""
	output_path = Path(output_distance_file)
	df = load_filtered_atoms_table(tsv_structure_file)
	validate_filtered_atoms_table(df, tsv_structure_file, require_nonempty=True)

	coords = df[["x", "y", "z"]].to_numpy(dtype=float)
	atom_id = df["atom_id"].to_numpy()
	resid = df["resid"].to_numpy()
	resname = df["resname"].astype(str).str.strip().str.upper().to_numpy()
	atom_name = df["atom_name"].astype(str).str.strip().str.upper().to_numpy()

	relaxed_torsion_pairs, backbone_dihedral_pairs = build_vdw_torsion_pairs(topology)
	excluded_pairs = topology["bonded_pairs"].union(topology["angle_pairs"]).union(backbone_dihedral_pairs)

	default_radius = 1.50
	hydrogen_soft_radius = 0.90
	torsion_scale = 0.8

	# Per-atom element radius; an H--H pair uses the soft radius for both atoms.
	elements = [name[0].upper() if name else "X" for name in atom_name]
	radius = np.array([VDW_RADII.get(element, default_radius) for element in elements], dtype=float)
	is_hydrogen = np.array([element == "H" for element in elements], dtype=bool)
	excluded_by_row = _pairs_by_row(excluded_pairs)
	relaxed_by_row = _pairs_by_row(relaxed_torsion_pairs)

	atom_ids = [int(value) for value in atom_id]
	resids = [int(value) for value in resid]
	names = [str(value) for value in atom_name]
	resnames = [str(value) for value in resname]
	row_format = "%d\t%d\t%d\t%d\t%.12f\t999.000000000000\t%s\t%s\t%s\t%s\n"
	table_i: list[np.ndarray] = []
	table_j: list[np.ndarray] = []
	table_lower: list[np.ndarray] = []

	with output_path.open("w", encoding="utf-8") as stream:
		stream.write(
			"atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\t"
			"atom_name_i\tatom_name_j\tresname_i\tresname_j\n"
		)

		for i in range(len(df)):
			# All partners j < i of atom i at once, in increasing j as in the pairwise loop.
			keep = np.ones(i, dtype=bool)
			keep[excluded_by_row.get(i, [])] = False
			partners = np.flatnonzero(keep)
			if partners.size == 0:
				continue

			delta = coords[i] - coords[partners]
			distance = np.sqrt(delta[:, 0] * delta[:, 0] + delta[:, 1] * delta[:, 1] + delta[:, 2] * delta[:, 2])
			both_hydrogen = is_hydrogen[i] & is_hydrogen[partners]
			radius_sum = np.where(
				both_hydrogen,
				hydrogen_soft_radius + hydrogen_soft_radius,
				radius[i] + radius[partners],
			)
			relaxed = np.zeros(i, dtype=bool)
			relaxed[relaxed_by_row.get(i, [])] = True
			base = np.where(relaxed[partners], torsion_scale * radius_sum, radius_sum)
			lower_bounds = lower_bounds_below(distance, base)
			table_i.append(np.full(partners.size, i))
			table_j.append(partners)
			table_lower.append(lower_bounds)

			stream.writelines(
				row_format % (
					atom_ids[i], atom_ids[j], resids[i], resids[j], lower_bound,
					names[i], names[j], resnames[i], resnames[j],
				)
				for j, lower_bound in zip(partners.tolist(), lower_bounds.tolist())
			)

	print_ok(f"Van der Waals distances saved to: {output_path}")
	rows_i = np.concatenate(table_i) if table_i else np.zeros(0, dtype=int)
	rows_j = np.concatenate(table_j) if table_j else np.zeros(0, dtype=int)
	lower = np.concatenate(table_lower) if table_lower else np.zeros(0)
	return distance_table(
		atom_id_i=atom_id[rows_i].astype(int), atom_id_j=atom_id[rows_j].astype(int),
		resid_i=resid[rows_i].astype(int), resid_j=resid[rows_j].astype(int),
		d_l=lower, d_u=np.full(len(lower), 999.0),
		atom_name_i=atom_name[rows_i], atom_name_j=atom_name[rows_j],
		resname_i=resname[rows_i], resname_j=resname[rows_j],
	)

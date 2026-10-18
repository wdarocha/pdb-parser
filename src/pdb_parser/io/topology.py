"""Covalent and planar topology helpers for distance-constraint generation."""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import TypedDict

import numpy as np

from pdb_parser.io.atom_tables import load_filtered_atoms_table, validate_filtered_atoms_table
from pdb_parser.io.files import distance_table
from pdb_parser.io.sequence import residue_sequence_from_ids
from pdb_parser.messages import print_ok

AtomKey = tuple[int, str]
AtomPair = tuple[int, int]
LogicalGraph = dict[AtomKey, set[AtomKey]]


class TopologyData(TypedDict):
	atom_index: dict[AtomKey, int]
	residue_ids: tuple[int, ...]
	logical_graph: LogicalGraph
	bonded_pairs: set[AtomPair]
	angle_pairs: set[AtomPair]
	excluded_pairs: set[AtomPair]


# Covalent bond templates for the 20 standard amino acids.
# Pairs are (atom_name1, atom_name2) within the same residue.
RESIDUE_BONDS = {
	"ALA": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB1"), ("CB", "HB2"), ("CB", "HB3"), ("CA", "C"), ("C", "O")},
	"GLY": {("N", "H"), ("N", "CA"), ("CA", "HA2"), ("CA", "HA3"), ("CA", "C"), ("C", "O")},
	"PRO": {("N", "CD"), ("CD", "HD2"), ("CD", "HD3"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "HG2"), ("CG", "HG3"), ("CG", "CD"), ("CA", "C"), ("C", "O")},
	"SER": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "OG"), ("OG", "HG"), ("CB", "HB2"), ("CB", "HB3"), ("CA", "C"), ("C", "O")},
	"THR": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "OG1"), ("OG1", "HG1"), ("CB", "CG2"), ("CB", "HB"), ("CG2", "HG21"), ("CG2", "HG22"), ("CG2", "HG23"), ("CA", "C"), ("C", "O")},
	"VAL": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CB", "HB"), ("CG1", "HG11"), ("CG1", "HG12"), ("CG1", "HG13"), ("CG2", "HG21"), ("CG2", "HG22"), ("CG2", "HG23"), ("CA", "C"), ("C", "O")},
	"LEU": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD11"), ("CD1", "HD12"), ("CD1", "HD13"), ("CD2", "HD21"), ("CD2", "HD22"), ("CD2", "HD23"), ("CG", "HG"), ("CA", "C"), ("C", "O")},
	"ILE": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CB", "HB"), ("CG1", "CD1"), ("CG1", "HG12"), ("CG1", "HG13"), ("CD1", "HD11"), ("CD1", "HD12"), ("CD1", "HD13"), ("CG2", "HG21"), ("CG2", "HG22"), ("CG2", "HG23"), ("CA", "C"), ("C", "O")},
	"MET": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "HG2"), ("CG", "HG3"), ("CG", "SD"), ("SD", "CE"), ("CE", "HE1"), ("CE", "HE2"), ("CE", "HE3"), ("CA", "C"), ("C", "O")},
	"CYS": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "SG"), ("SG", "HG"), ("CA", "C"), ("C", "O")},
	"ASN": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2"), ("ND2", "HD21"), ("ND2", "HD22"), ("CA", "C"), ("C", "O")},
	"GLN": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "HG2"), ("CG", "HG3"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2"), ("NE2", "HE21"), ("NE2", "HE22"), ("CA", "C"), ("C", "O")},
	"ASP": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2"), ("CA", "C"), ("C", "O")},
	"GLU": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "HG2"), ("CG", "HG3"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2"), ("CA", "C"), ("C", "O")},
	"HIS": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "HD1"), ("CD2", "HD2"), ("CE1", "HE1"), ("NE2", "HE2"), ("ND1", "CE1"), ("CE1", "NE2"), ("CD2", "NE2"), ("CA", "C"), ("C", "O")},
	"PHE": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD2", "HD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CE1", "HE1"), ("CE2", "HE2"), ("CZ", "HZ"), ("CA", "C"), ("C", "O")},
	"TYR": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD2", "HD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CE1", "HE1"), ("CE2", "HE2"), ("CZ", "OH"), ("OH", "HH"), ("CA", "C"), ("C", "O")},
	"TRP": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD1", "NE1"), ("NE1", "HE1"), ("NE1", "CE2"), ("CD2", "CE2"), ("CD2", "CE3"), ("CE2", "CZ2"), ("CE3", "HE3"), ("CE3", "CZ3"), ("CZ2", "HZ2"), ("CZ2", "CH2"), ("CZ3", "HZ3"), ("CZ3", "CH2"), ("CH2", "HH2"), ("CA", "C"), ("C", "O")},
	"LYS": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "HG2"), ("CG", "HG3"), ("CG", "CD"), ("CD", "HD2"), ("CD", "HD3"), ("CD", "CE"), ("CE", "HE2"), ("CE", "HE3"), ("CE", "NZ"), ("NZ", "HZ1"), ("NZ", "HZ2"), ("NZ", "HZ3"), ("CA", "C"), ("C", "O")},
	"ARG": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "HG2"), ("CG", "HG3"), ("CG", "CD"), ("CD", "HD2"), ("CD", "HD3"), ("CD", "NE"), ("NE", "HE"), ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2"), ("NH1", "HH11"), ("NH1", "HH12"), ("NH2", "HH21"), ("NH2", "HH22"), ("CA", "C"), ("C", "O")},
}


def build_covalent_and_planar_topology(
	tsv_structure_file: str | Path,
) -> TopologyData:
	"""Build reusable covalent/planar topology structures from a filtered atom table."""
	df = load_filtered_atoms_table(tsv_structure_file)
	validate_filtered_atoms_table(df, tsv_structure_file, require_nonempty=True)
	unknown_residues = sorted(
		set(df["resname"].astype(str).str.strip().str.upper()) - set(RESIDUE_BONDS)
	)
	if unknown_residues:
		raise ValueError(
			"Filtered atom table contains protein residues without covalent topology "
			f"templates: {', '.join(unknown_residues)}."
		)

	atom_index: dict[AtomKey, int] = {
		(int(residue_id), str(atom_name).strip().upper()): int(row_index)
		for row_index, residue_id, atom_name in zip(df.index, df["resid"], df["atom_name"], strict=False)
	}

	logical_graph: LogicalGraph = defaultdict(set)
	resid_array = df["resid"].to_numpy()

	for residue_id, group in df.groupby("resid"):
		residue_id = int(residue_id)
		residue_name = str(group.iloc[0]["resname"]).strip().upper()
		names_present = set(group["atom_name"].astype(str).str.strip().str.upper())

		pairs = set(RESIDUE_BONDS.get(residue_name, set()))

		for hydrogen_name in ("H", "H1", "H2", "H3", "HD2", "HD3"):
			if hydrogen_name in names_present:
				pairs.add(("N", hydrogen_name))

		if "OXT" in names_present:
			pairs.add(("C", "OXT"))

		for atom_name_1, atom_name_2 in pairs:
			key_1 = (residue_id, atom_name_1)
			key_2 = (residue_id, atom_name_2)

			if key_1 in atom_index and key_2 in atom_index:
				logical_graph[key_1].add(key_2)
				logical_graph[key_2].add(key_1)

	all_residues = residue_sequence_from_ids(resid_array, context="filtered atom table")

	for residue_1, residue_2 in zip(all_residues, all_residues[1:]):
		# The filtered-table loader has already checked continuity.  Keep this
		# explicit guard close to the peptide-bond construction as a defense for
		# callers that provide a topology object assembled elsewhere.
		if residue_2 != residue_1 + 1:
			raise ValueError(
				f"Cannot infer covalent peptide bond across residue gap {residue_1}->{residue_2}."
			)
		key_1 = (residue_1, "C")
		key_2 = (residue_2, "N")

		if key_1 in atom_index and key_2 in atom_index:
			logical_graph[key_1].add(key_2)
			logical_graph[key_2].add(key_1)

	bonded_pairs: set[AtomPair] = set()

	for atom_key, neighbor_keys in logical_graph.items():
		i = atom_index.get(atom_key)
		if i is None:
			continue

		for neighbor_key in neighbor_keys:
			j = atom_index.get(neighbor_key)
			if j is None:
				continue

			pair = (i, j) if i > j else (j, i)
			bonded_pairs.add(pair)

	angle_pairs: set[AtomPair] = set()

	for neighbor_keys in logical_graph.values():
		neighbor_indices = [atom_index[key] for key in neighbor_keys if key in atom_index]

		for i, j in combinations(sorted(neighbor_indices), 2):
			pair = (i, j) if i > j else (j, i)

			if pair not in bonded_pairs:
				angle_pairs.add(pair)

	excluded_pairs = bonded_pairs.union(angle_pairs)

	return {
		"atom_index": atom_index,
		"residue_ids": tuple(all_residues),
		"logical_graph": logical_graph,
		"bonded_pairs": bonded_pairs,
		"angle_pairs": angle_pairs,
		"excluded_pairs": excluded_pairs,
	}


def covalent_and_planar_distances(
	tsv_structure_file: str | Path,
	topology: TopologyData,
	output_distance_file: str | Path,
):
	"""Write exact distances for covalent and planar pairs and return them as a table."""
	output_path = Path(output_distance_file)
	df = load_filtered_atoms_table(tsv_structure_file)
	validate_filtered_atoms_table(df, tsv_structure_file, require_nonempty=True)

	coords = df[["x", "y", "z"]].to_numpy(dtype=float)
	atom_id = df["atom_id"].to_numpy()
	resid = df["resid"].to_numpy()
	resname = df["resname"].astype(str).str.strip().str.upper().to_numpy()
	atom_name = df["atom_name"].astype(str).str.strip().str.upper().to_numpy()

	all_pairs = sorted(topology["bonded_pairs"].union(topology["angle_pairs"]))

	with output_path.open("w", encoding="utf-8") as file:
		file.write(
			"atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\t"
			"atom_name_i\tatom_name_j\tresname_i\tresname_j\n"
		)

		records = []
		for i, j in all_pairs:
			distance = float(np.linalg.norm(coords[i] - coords[j]))
			records.append((int(atom_id[i]), int(atom_id[j]), int(resid[i]), int(resid[j]), distance, distance,
				atom_name[i], atom_name[j], resname[i], resname[j]))
			file.write(
				f"{int(atom_id[i])}\t{int(atom_id[j])}\t{int(resid[i])}\t{int(resid[j])}\t"
				f"{distance:.12f}\t{distance:.12f}\t{atom_name[i]}\t{atom_name[j]}\t"
				f"{resname[i]}\t{resname[j]}\n"
			)

	print_ok(f"Covalent/planar pair distances saved to: {output_path}")
	return distance_table(records)

"""Distance constraints induced by peptide-plane geometry."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np

from pdb_parser.io.atom_tables import load_filtered_atoms_table, validate_filtered_atoms_table
from pdb_parser.io.files import distance_table
from pdb_parser.io.sequence import residue_sequence_from_ids
from pdb_parser.io.topology import AtomPair
from pdb_parser.messages import print_ok


def planar_peptide_distances(
	tsv_structure_file: str | Path,
	output_distance_file: str | Path,
	atom_selection: str,
):
	"""Generate exact peptide-plane distance pairs from a filtered TSV."""
	output_path = Path(output_distance_file)
	df = load_filtered_atoms_table(tsv_structure_file)
	validate_filtered_atoms_table(df, tsv_structure_file, require_nonempty=True)

	coords = df[["x", "y", "z"]].to_numpy(dtype=float)
	atom_id = df["atom_id"].to_numpy()
	resid = df["resid"].to_numpy()
	resname = df["resname"].str.strip().str.upper().to_numpy()
	atom_name = df["atom_name"].str.strip().str.upper().to_numpy()
	atom_index = {
		(int(residue_id), str(name).strip().upper()): int(row_index)
		for row_index, residue_id, name in zip(
			df.index, df["resid"], df["atom_name"], strict=False
		)
	}

	def find_any(residue_id: int, names: Iterable[str]) -> set[int]:
		return {
			index
			for name in names
			if (index := atom_index.get((int(residue_id), name))) is not None
		}

	pairs: set[AtomPair] = set()

	def add_pair(
		name1: Iterable[str] | str,
		res1: int,
		name2: Iterable[str] | str,
		res2: int,
	) -> None:
		names1 = (name1,) if isinstance(name1, str) else tuple(name1)
		names2 = (name2,) if isinstance(name2, str) else tuple(name2)
		for index1 in find_any(res1, names1):
			for index2 in find_any(res2, names2):
				first_atom = int(atom_id[index1])
				second_atom = int(atom_id[index2])
				if first_atom != second_atom:
					pairs.add((first_atom, second_atom) if first_atom > second_atom else (second_atom, first_atom))

	all_residues = residue_sequence_from_ids(resid, context="filtered atom table")
	for residue1, residue2 in zip(all_residues, all_residues[1:], strict=False):
		residue2_name = df.loc[df["resid"] == residue2, "resname"].iloc[0].strip().upper()
		add_pair("O", residue1, "CA", residue2)
		if residue2_name != "PRO":
			add_pair("CA", residue1, "H", residue2)
			add_pair("O", residue1, "H", residue2)
		elif atom_selection == "backbone_plus_hydrogens":
			for hydrogen_name in ("HD2", "HD3"):
				add_pair("CA", residue1, hydrogen_name, residue2)
				add_pair("O", residue1, hydrogen_name, residue2)
		else:
			add_pair("CA", residue1, "CD", residue2)
			add_pair("O", residue1, "CD", residue2)

	with output_path.open("w", encoding="utf-8") as output_file:
		output_file.write("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n")
		records = []
		for first_atom, second_atom in sorted(pairs):
			index1 = int(np.flatnonzero(atom_id == first_atom)[0])
			index2 = int(np.flatnonzero(atom_id == second_atom)[0])

			distance = float(np.linalg.norm(coords[index1] - coords[index2]))
			records.append((int(first_atom), int(second_atom), int(resid[index1]), int(resid[index2]), distance,
				distance, atom_name[index1], atom_name[index2], resname[index1], resname[index2]))
			output_file.write(
				f"{first_atom}\t{second_atom}\t{int(resid[index1])}\t{int(resid[index2])}\t"
				f"{distance:.12f}\t{distance:.12f}\t{atom_name[index1]}\t{atom_name[index2]}\t"
				f"{resname[index1]}\t{resname[index2]}\n"
			)

	print_ok(f"Peptide-plane distances saved to: {output_distance_file}")
	return distance_table(records)

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Iterable, Tuple, Set, Dict
from pdb_parser.utils import ensure_dir
from pdb_parser.geometry.distance_geometry import torsion_angle_2_endpoint_distance

import numpy as np
import pandas as pd
import MDAnalysis as mda

from MDAnalysis.lib.distances import calc_dihedrals

# -----------------------------------------------------------------------------------------------------
# Covalent bond templates for the 20 standard amino acids
# Pairs are (atom_name1, atom_name2) within the *same* residue.
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
	"HIS": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "HD1"), ("CD2", "HD2"), ("CE1", "HE1"), ("NE2", "HE2"), ("CA", "C"), ("C", "O")},
	"PHE": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD2", "HD2"), ("CE1", "HE1"), ("CE2", "HE2"), ("CZ", "HZ"), ("CA", "C"), ("C", "O")},
	"TYR": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD2", "HD2"), ("CE1", "HE1"), ("CE2", "HE2"), ("CZ", "OH"), ("OH", "HH"), ("CA", "C"), ("C", "O")},
	"TRP": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "CD1"), ("CD1", "HD1"), ("CD2", "CE2"), ("CE2", "HE2"), ("CE3", "HE3"), ("CZ2", "HZ2"), ("CZ3", "HZ3"), ("CH2", "HH2"), ("CA", "C"), ("C", "O")},
	"LYS": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "HG2"), ("CG", "HG3"), ("CG", "CD"), ("CD", "HD2"), ("CD", "HD3"), ("CD", "CE"), ("CE", "HE2"), ("CE", "HE3"), ("CE", "NZ"), ("NZ", "HZ1"), ("NZ", "HZ2"), ("NZ", "HZ3"), ("CA", "C"), ("C", "O")},
	"ARG": {("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "HG2"), ("CG", "HG3"), ("CG", "CD"), ("CD", "HD2"), ("CD", "HD3"), ("CD", "NE"), ("NE", "HE"), ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2"), ("NH1", "HH11"), ("NH1", "HH12"), ("NH2", "HH21"), ("NH2", "HH22"), ("CA", "C"), ("C", "O")}
}
# -----------------------------------------------------------------------------------------------------
def load_filtered_atoms_table(
	tsv_structure_file: str | Path,
) -> pd.DataFrame:
	"""
	Load and normalize a filtered atom TSV file.
	"""
	tsv_structure_file = Path(tsv_structure_file)

	if not tsv_structure_file.exists():
		raise ValueError(f"Filtered input file not found: {tsv_structure_file}")

	dataframe = pd.read_csv(tsv_structure_file, sep="\t", dtype={"atom_name": str, "resname": str})
	dataframe.columns = dataframe.columns.str.strip()

	return dataframe
# -----------------------------------------------------------------------------------------------------
def build_covalent_and_planar_topology(
	tsv_structure_file: str | Path,
) -> dict[str, object]:
	"""
	Build reusable covalent/planar topology structures from a filtered atom table.

	Returned dictionary keys
	------------------------
	atom_index
		Mapping (resid, atom_name) -> df row index.

	logical_graph
		Covalent adjacency graph using keys of the form (resid, atom_name).

	bonded_pairs
		Set of atom pairs separated by exactly 1 covalent bond, stored as
		(max_row_index, min_row_index).

	angle_pairs
		Set of atom pairs separated by exactly 2 covalent bonds, stored as
		(max_row_index, min_row_index).

	excluded_pairs
		Union of bonded_pairs and angle_pairs.
	"""
	df = load_filtered_atoms_table(tsv_structure_file)
	
	atom_index: dict[tuple[int, str], int] = {
		(int(row.resid), str(row.atom_name).strip().upper()): i
		for i, row in df.iterrows()
	}

	logical_graph: dict[tuple[int, str], set[tuple[int, str]]] = defaultdict(set)

	resid_array = df["resid"].to_numpy()

	for residue_id, group in df.groupby("resid"):
		residue_id = int(residue_id)
		residue_name = str(group.iloc[0]["resname"]).strip().upper()
		names_present = set(group["atom_name"].astype(str).str.strip().str.upper())

		pairs = set(RESIDUE_BONDS.get(residue_name, set()))

		# N-terminus H variants if present
		for hydrogen_name in ("H", "H1", "H2", "H3", "HD2"): # HD2 for PRO
			if hydrogen_name in names_present:
				pairs.add(("N", hydrogen_name))

		# C-terminus OXT if present
		if "OXT" in names_present:
			pairs.add(("C", "OXT"))

		for atom_name_1, atom_name_2 in pairs:
			key_1 = (residue_id, atom_name_1)
			key_2 = (residue_id, atom_name_2)

			if key_1 in atom_index and key_2 in atom_index:
				logical_graph[key_1].add(key_2)
				logical_graph[key_2].add(key_1)

	all_residues = sorted(set(int(value) for value in resid_array))

	for residue_1, residue_2 in zip(all_residues, all_residues[1:]):
		key_1 = (residue_1, "C")
		key_2 = (residue_2, "N")

		if key_1 in atom_index and key_2 in atom_index:
			logical_graph[key_1].add(key_2)
			logical_graph[key_2].add(key_1)

	bonded_pairs: set[tuple[int, int]] = set()

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

	angle_pairs: set[tuple[int, int]] = set()

	for middle_key, neighbor_keys in logical_graph.items():
		neighbor_indices = [atom_index[key] for key in neighbor_keys if key in atom_index]

		for i, j in combinations(sorted(neighbor_indices), 2):
			pair = (i, j) if i > j else (j, i)

			if pair not in bonded_pairs:
				angle_pairs.add(pair)

	excluded_pairs = bonded_pairs.union(angle_pairs)

	return {
		"atom_index": atom_index,
		"logical_graph": logical_graph,
		"bonded_pairs": bonded_pairs,
		"angle_pairs": angle_pairs,
		"excluded_pairs": excluded_pairs,
	}
# -----------------------------------------------------------------------------------------------------
def covalent_and_planar_distances(
	tsv_structure_file: str | Path,
	topology: dict[str, object],
	output_distance_file: str | Path,
) -> None:
	"""
	Write the exact distances for covalent (1-bond) and planar (2-bond) pairs.

	Input df columns (required)
	----------------------------------
	atom_id, atom_name, resid, resname, x, y, z
	"""
	df = load_filtered_atoms_table(tsv_structure_file)
	
	if df.empty:
		print(f"[OK] Covalent/planar pair distances saved to (empty): {output_distance_file}")
		output_distance_file.write_text("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n", encoding="utf-8")
		return

	coords = df[["x", "y", "z"]].to_numpy(dtype=float)
	atom_id = df["atom_id"].to_numpy()
	resid = df["resid"].to_numpy()
	resname = df["resname"].astype(str).str.strip().str.upper().to_numpy()
	atom_name = df["atom_name"].astype(str).str.strip().str.upper().to_numpy()

	bonded_pairs = topology["bonded_pairs"]
	angle_pairs = topology["angle_pairs"]
	all_pairs = sorted(bonded_pairs.union(angle_pairs))

	with output_distance_file.open("w", encoding="utf-8") as file:
		file.write("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n")

		for i, j in all_pairs:
			distance = float(np.linalg.norm(coords[i] - coords[j]))

			file.write(f"{int(atom_id[i])}\t{int(atom_id[j])}\t{int(resid[i])}\t{int(resid[j])}\t{distance:.16f}\t{distance:.16f}\t{atom_name[i]}\t{atom_name[j]}\t{resname[i]}\t{resname[j]}\n")

	print(f"[OK] Covalent/planar pair distances saved to: {output_distance_file}")
# -----------------------------------------------------------------------------------------------------
def vdw_distances(
	tsv_structure_file: str | Path,
	topology: dict[str, object],
	output_distance_file: str | Path,
) -> None:
	"""
	Write van der Waals lower-bound constraints for pairs not separated by
	1 or 2 covalent bonds.

	Input df columns (required)
	----------------------------------
	atom_id, atom_name, resid, resname, x, y, z
	"""
	df = load_filtered_atoms_table(tsv_structure_file)
	
	if df.empty:
		print(f"[OK] Van der Waals distance table saved to (empty): {output_distance_file}")
		output_distance_file.write_text("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n", encoding="utf-8")
		return

	coords = df[["x", "y", "z"]].to_numpy(dtype=float)
	atom_id = df["atom_id"].to_numpy()
	resid = df["resid"].to_numpy()
	resname = df["resname"].astype(str).str.strip().str.upper().to_numpy()
	atom_name = df["atom_name"].astype(str).str.strip().str.upper().to_numpy()

	excluded_pairs = topology["excluded_pairs"]

	vdw_radii = {
		"H": 1.20,
		"C": 1.70,
		"N": 1.55,
		"O": 1.52,
	}
	default_radius = 1.50
	hydrogen_soft_radius = 0.90

	with output_distance_file.open("w", encoding="utf-8") as f:
		f.write("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n")

		for i in range(len(df)):
			for j in range(i):
				if (i, j) in excluded_pairs:
					continue

				distance = float(np.linalg.norm(coords[i] - coords[j]))

				element_i = atom_name[i][0] if atom_name[i] else "X"
				element_j = atom_name[j][0] if atom_name[j] else "X"
				element_i = element_i.upper()
				element_j = element_j.upper()

				if element_i == "H" and element_j == "H":
					radius_i = hydrogen_soft_radius
					radius_j = hydrogen_soft_radius
				else:
					radius_i = vdw_radii.get(element_i, default_radius)
					radius_j = vdw_radii.get(element_j, default_radius)

				radius_sum = radius_i + radius_j
				lower_bound = radius_sum

				if distance < lower_bound:
					for percentage in range(90, 0, -10):
						candidate = (percentage / 100.0) * radius_sum
						if candidate < distance:
							lower_bound = candidate
							break

				upper_bound = 999.0

				f.write(f"{int(atom_id[i])}\t{int(atom_id[j])}\t{int(resid[i])}\t{int(resid[j])}\t{lower_bound:.16f}\t{upper_bound:.16f}\t{atom_name[i]}\t{atom_name[j]}\t{resname[i]}\t{resname[j]}\n")

	print(f"[OK] Van der Waals distances saved to: {output_distance_file}")
# -----------------------------------------------------------------------------------------------------
def planar_peptide_distances(
	tsv_structure_file: str | Path,
	output_distance_file: str | Path
) -> None:
	"""
	Generate planar peptide-group distance pairs from a filtered TSV.

	Input TSV columns (required):
		atom_id, atom_name, resid, resname, x, y, z

	Rules
	-----
	For consecutive residues r1 -> r2, include (if present):
	  - For all resname include:
	  	(O_r1, CA_r2)
	  - If resname(r2) != PRO include: 
	  	(CA_r1, H_r2), (O_r1, H_r2)
	  - If resname(r2) == PRO include:
	  	(CA_r1, CD_r2), (O_r1, CD_r2)
	
	Distances are written as exact: d_l == d_u == ||ri - rj|| (Å).
	"""
	df = load_filtered_atoms_table(tsv_structure_file)

	if df.empty:
		# Nothing to do
		print(f"[OK] Peptide-plane distances saved to (empty): {output_distance_file}")
		output_distance_file.write_text("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n", encoding="utf-8")
		return

	# Extract arrays
	coords = df[["x", "y", "z"]].to_numpy(dtype=float)
	atom_id = df["atom_id"].to_numpy()
	resid = df["resid"].to_numpy()
	resname = df["resname"].str.strip().str.upper().to_numpy()
	atom_name = df["atom_name"].str.strip().str.upper().to_numpy()

	# Index helpers
	# (resid, atom_name) -> row index
	atom_index: dict[tuple[int, str], int] = {
		(int(row.resid), row.atom_name.strip().upper()): i
		for i, row in df.iterrows()
	}

	def _find_any(residue_id: int, names: Iterable[str]) -> set[int]:
		"""
		Return the set of row indices whose atom names match any candidate name.
		"""
		found: set[int] = set()

		for name in names:
			idx = atom_index.get((int(residue_id), name))
			if idx is not None:
				found.add(idx)

		return found

	pairs: set[tuple[int, int]] = set()

	def _add_pair(
		name1: Iterable[str] | str,
		res1: int,
		name2: Iterable[str] | str,
		res2: int,
	) -> None:
		"""
		Add all existing atom pairs between two residue/name specifications.
		"""
		names1 = (name1,) if isinstance(name1, str) else tuple(name1)
		names2 = (name2,) if isinstance(name2, str) else tuple(name2)

		idx1 = _find_any(res1, names1)
		idx2 = _find_any(res2, names2)

		if not idx1 or not idx2:
			return

		for i1 in idx1:
			for i2 in idx2:
				a1 = int(atom_id[i1])
				a2 = int(atom_id[i2])

				if a1 == a2:
					continue

				pair = (a1, a2) if a1 > a2 else (a2, a1)
				pairs.add(pair)

	# Iterate consecutive residues
	all_residues = sorted(set(int(value) for value in resid))

	for r1, r2 in zip(all_residues, all_residues[1:]):
		r2_name = df.loc[df["resid"] == r2, "resname"].iloc[0].strip().upper()
		
		# Base peptide-plane related pairs
		_add_pair("O", r1, "CA", r2)
		# nonPRO-specific pairs involving residue r2
		if r2_name != "PRO":
			_add_pair("CA", r1, "H", r2)
			_add_pair("O", r1, "H", r2)
		# PRO-specific pairs involving residue r2
		if r2_name == "PRO":
			_add_pair("CA", r1, "CD", r2)
			_add_pair("O", r1, "CD", r2)
		
	# Write output file: distances as tight bounds (d_l = d_u = distance)
	with output_distance_file.open("w", encoding="utf-8") as f:
		f.write("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n")

		for a1, a2 in sorted(pairs):
			i1 = int(a1) - 1
			i2 = int(a2) - 1

			if not (0 <= i1 < len(atom_id)) or int(atom_id[i1]) != a1:
				i1 = int(np.where(atom_id == a1)[0][0])

			if not (0 <= i2 < len(atom_id)) or int(atom_id[i2]) != a2:
				i2 = int(np.where(atom_id == a2)[0][0])

			distance = float(np.linalg.norm(coords[i1] - coords[i2]))

			f.write(f"{a1}\t{a2}\t{int(resid[i1])}\t{int(resid[i2])}\t{distance:.16f}\t{distance:.16f}\t{atom_name[i1]}\t{atom_name[i2]}\t{resname[i1]}\t{resname[i2]}\n")

	print(f"[OK] Peptide-plane distances saved to: {output_distance_file}")
# -----------------------------------------------------------------------------------------------------
def _get_nmr_candidate_indices(atom_name: np.ndarray, atom_selection: str) -> np.ndarray:
	"""Return indices of atoms eligible for NMR-derived distance constraints."""
	atom_selection = (atom_selection or "").strip().lower()

	if atom_selection == "backbone":
		# In backbone-only mode, hydrogens are not present.
		# We map HN -> N and HA -> C.
		mask = np.isin(atom_name, ["N", "C"])
	elif atom_selection in {"backbone_plus_hydrogens", "full_chain", "backbone_plus_neighbors"}:
		# In these modes, NMR restraints are applied directly to hydrogens.
		mask = np.char.startswith(atom_name.astype(str), "H")
	else:
		raise ValueError(f"Unsupported atom_selection: {atom_selection!r}. Expected 'full_chain', 'backbone', 'backbone_plus_hydrogens', or 'backbone_plus_neighbors'.")

	return np.flatnonzero(mask)

def _get_centered_interval(
	dist: float,
	resid_i: int,
	resid_j: int,
	epsilon_short: float,
	epsilon_long: float,
	vdw_threshold: float,
	max_distance: float,
) -> tuple[float, float] | None:
	"""Return a centered synthetic interval around the reference distance."""
	residue_gap = abs(int(resid_i) - int(resid_j))
	epsilon = epsilon_short if residue_gap <= 1 else epsilon_long

	sampled_distance = float(np.random.normal(loc=dist, scale=epsilon / 8.0))

	d_l = max(sampled_distance - epsilon / 2.0, vdw_threshold)
	d_u = min(sampled_distance + epsilon / 2.0, max_distance)

	if d_l > d_u:
		return None

	return d_l, d_u

def _get_experimental_interval(
	dist: float,
	noe_strong: float,
	noe_medium: float,
	noe_weak: float,
	vdw_threshold: float,
) -> tuple[float, float] | None:
	"""Return an experimental-style interval based on NOE intensity classes."""
	if dist < noe_strong:
		return vdw_threshold, noe_strong
	if dist < noe_medium:
		return vdw_threshold, noe_medium
	if dist < noe_weak:
		return vdw_threshold, noe_weak

	return None

def _get_precise_interval(
	dist: float,
	max_distance: float,
) -> tuple[float, float] | None:
	"""Return a degenerate interval [d, d] if the pair is within the cutoff."""
	if dist >= max_distance:
		return None

	return dist, dist

def nmr_distance_constraints(
	tsv_structure_file: str | Path,
	output_distance_file: str | Path,
	distance_constraints: str,
	atom_selection: str,
	epsilon_short: float,
	epsilon_long: float,
	max_distance: float,
	noe_strong: float,
	noe_medium: float,
	noe_weak: float,
	vdw_threshold: float,
) -> None:
	"""Generate NMR-derived distance constraints from a TSV structure file.

	Input TSV must have columns:
		atom_id, atom_name, resid, resname, x, y, z

	Distance modes
	--------------
	- precise
		Use exact geometric distances, i.e., [d, d].

	- interval_centered
		Use synthetic intervals around the reference distance d_ij:
			D_ij = [max(d*_ij - epsilon_ij / 2, vdw_threshold),
			        min(d*_ij + epsilon_ij / 2, max_distance)]
		where
			d*_ij ~ N(d_ij, (epsilon_ij / 8)^2)

	- interval_experimental
		Use NOESY-style intervals:
			strong: [vdw_threshold, noe_strong]
			medium: [vdw_threshold, noe_medium]
			weak: [vdw_threshold, noe_weak]

	Atom-selection behavior
	-----------------------
	- If atom_selection == "backbone":
		NMR restraints are generated between mapped heavy atoms:
			N plays the role of HN
			C plays the role of HA
		thus producing N-N, N-C, and C-C pairs.

	- Otherwise:
		NMR restraints are generated directly between hydrogen atoms.
	"""
	distance_constraints = (distance_constraints or "").strip().lower()
	valid_modes = {"precise", "interval_centered", "interval_experimental"}

	if distance_constraints not in valid_modes:
		raise ValueError(f"Unsupported distance_constraints mode: {distance_constraints!r}. Expected one of {sorted(valid_modes)}.")

	df = load_filtered_atoms_table(tsv_structure_file)

	if df.empty:
		print(f"[OK] NMR distance table (mode={distance_constraints}) saved to (empty): {output_distance_file}")
		output_distance_file.write_text("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n", encoding="utf-8")
		return

	coords = df[["x", "y", "z"]].to_numpy(dtype=float)
	atom_id = df["atom_id"].to_numpy()
	resid = df["resid"].to_numpy()
	resname = df["resname"].str.strip().str.upper().to_numpy()
	atom_name = df["atom_name"].str.strip().str.upper().to_numpy()

	candidate_indices = _get_nmr_candidate_indices(atom_name, atom_selection)
	
	with output_distance_file.open("w", encoding="utf-8") as f:
		f.write("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n")

		if len(candidate_indices) == 0:
			print(f"[OK] NMR distance table (mode={distance_constraints}) saved to (empty): {out_path}")
			return

		for pos_i in range(len(candidate_indices)):
			i = int(candidate_indices[pos_i])

			for pos_j in range(pos_i):
				j = int(candidate_indices[pos_j])

				dist = float(np.linalg.norm(coords[i] - coords[j]))

				if distance_constraints == "precise":
					interval = _get_precise_interval(dist, max_distance)
					
				elif distance_constraints == "interval_centered":
					if dist >= max_distance:
						continue

					interval = _get_centered_interval(dist, int(resid[i]), int(resid[j]), epsilon_short, epsilon_long, vdw_threshold, max_distance)
				
				else:
					interval = _get_experimental_interval(dist, noe_strong, noe_medium, noe_weak, vdw_threshold)

				if interval is None:
					continue

				d_l, d_u = interval

				f.write(f"{int(atom_id[i])}\t{int(atom_id[j])}\t{int(resid[i])}\t{int(resid[j])}\t{d_l:.16f}\t{d_u:.16f}\t{atom_name[i]}\t{atom_name[j]}\t{resname[i]}\t{resname[j]}\n")

	print(f"[OK] NMR-derived distances (mode={distance_constraints}) saved to: {output_distance_file}")
# -----------------------------------------------------------------------------------------------------
def wrap_angle_deg(angle: float) -> float:
	"""
	Wrap an angle (in degrees) to the range [-180, 180).
	"""
	return (angle + 180.0) % 360.0 - 180.0


def full_range_angle() -> tuple[float, float]:
	"""
	Return a full-range torsion interval representing [-180, 180].
	"""
	return 0.0, 180.0

def talos_n_like(
	pdb_file: str | Path,
	chosen_model: int,
	chosen_chain: str,
	output_angular_file: str | Path,
	distance_constraints: str,
	omega_angle_width: float,
	phi_angle_width: float,
	psi_angle_width: float,
	percentage_backbone_torsion_angles: float,
) -> None:
	"""
	Generate TALOS-N-like omega/phi/psi intervals from a PDB structure using
	MDAnalysis, computing backbone dihedrals manually from N/CA/C coordinates.

	For each backbone residue:
	- omega is always treated as detected whenever it can be computed
	- phi and psi are randomly split into detected/full-range according to
	  percentage_backbone_torsion_angles
	- If a valid torsion can be computed, the central value is sampled from
	  Normal(mean = structural angle [deg], sigma = angle_width / 8),
	  and the interval radius is angle_width / 2
	- If a given torsion cannot be computed, that angle is encoded as a
	  full-range interval

	The distance_constraints mode is interpreted as follows:
	- precise:
		all torsion widths are forced to 0°, producing zero-width intervals
	- interval_centered / interval_experimental:
		the provided torsion widths are used as given

	Output format (TSV)
	-------------------
	resid
	resname
	omega_center
	omega_radius
	phi_center
	phi_radius
	psi_center
	psi_radius
	"""
	def sanitize_angle_width(angle_width: float, angle_name: str) -> float:
		"""
		Convert the input width to a valid non-negative float.
		"""
		try:
			angle_width = float(angle_width)
		except (TypeError, ValueError):
			angle_width = 360.0

		if angle_width < 0.0:
			angle_width = 360.0
			print(f"[WARN] Negative {angle_name} detected; using full range (360°).")

		return angle_width

	def sample_angle_interval(base_angle: float | None, angle_width: float) -> tuple[float, float]:
		"""
		Build a torsion interval from a structural torsion angle.
		"""
		if base_angle is None:
			return full_range_angle()

		if angle_width <= 0.0:
			return base_angle, 0.0

		angle_radius = angle_width / 2.0
		angle_sigma = angle_width / 8.0
		angle_center = wrap_angle_deg(float(np.random.normal(loc=base_angle, scale=angle_sigma)))

		return angle_center, angle_radius

	# ---------------------------
	# Parameter checks
	# ---------------------------
	distance_constraints = (distance_constraints or "").strip().lower()
	valid_modes = {"precise", "interval_centered", "interval_experimental"}

	if distance_constraints not in valid_modes:
		raise ValueError(f"Unsupported distance_constraints value: {distance_constraints!r}. Expected one of {sorted(valid_modes)}.")

	omega_angle_width = sanitize_angle_width(omega_angle_width, "omega_angle_width")
	phi_angle_width = sanitize_angle_width(phi_angle_width, "phi_angle_width")
	psi_angle_width = sanitize_angle_width(psi_angle_width, "psi_angle_width")

	if distance_constraints == "precise":
		omega_angle_width = 0.0
		phi_angle_width = 0.0
		psi_angle_width = 0.0

	try:
		percentage_backbone_torsion_angles = float(percentage_backbone_torsion_angles)
	except (TypeError, ValueError):
		percentage_backbone_torsion_angles = 100.0

	if percentage_backbone_torsion_angles < 0.0:
		print("[WARN] percentage_backbone_torsion_angles < 0; clamped to 0%.")
		percentage_backbone_torsion_angles = 0.0
	elif percentage_backbone_torsion_angles > 100.0:
		print("[WARN] percentage_backbone_torsion_angles > 100; clamped to 100%.")
		percentage_backbone_torsion_angles = 100.0

	# ---------------------------
	# Load PDB and select chain
	# ---------------------------
	universe = mda.Universe(str(pdb_file))

	if chosen_model is not None:
		frame_index = int(chosen_model) - 1
		if frame_index < 0 or frame_index >= len(universe.trajectory):
			raise ValueError(f"chosen_model={chosen_model} is out of range for trajectory with {len(universe.trajectory)} frame(s).")
		universe.trajectory[frame_index]

	chain_residues = None
	selection_candidates: list[str] = []

	if chosen_chain:
		selection_candidates.append(f"protein and segid {chosen_chain}")
		selection_candidates.append(f"protein and chainid {chosen_chain}")

	selection_candidates.append("protein")

	for selection in selection_candidates:
		atom_group = universe.select_atoms(selection)
		if len(atom_group.residues) > 0:
			chain_residues = atom_group.residues
			break

	if chain_residues is None or len(chain_residues) == 0:
		raise ValueError(f"No residues found for chain '{chosen_chain}' in structure '{pdb_file}'.")

	# ---------------------------
	# Build backbone coordinate arrays
	# ---------------------------
	residue_ids: list[int] = []
	residue_names: list[str] = []
	n_positions: list[np.ndarray | None] = []
	ca_positions: list[np.ndarray | None] = []
	c_positions: list[np.ndarray | None] = []

	for residue in chain_residues:
		residue_ids.append(int(residue.resid))
		residue_names.append(str(residue.resname).upper())

		n_atoms = residue.atoms.select_atoms("name N")
		ca_atoms = residue.atoms.select_atoms("name CA")
		c_atoms = residue.atoms.select_atoms("name C")

		n_positions.append(n_atoms[0].position.copy() if len(n_atoms) == 1 else None)
		ca_positions.append(ca_atoms[0].position.copy() if len(ca_atoms) == 1 else None)
		c_positions.append(c_atoms[0].position.copy() if len(c_atoms) == 1 else None)

	n_residues = len(residue_ids)

	# ---------------------------
	# Compute backbone torsions
	# ---------------------------
	omega_angles: list[float | None] = [None] * n_residues
	phi_angles: list[float | None] = [None] * n_residues
	psi_angles: list[float | None] = [None] * n_residues

	for i in range(n_residues):
		if i > 0:
			ca_prev = ca_positions[i - 1]
			c_prev = c_positions[i - 1]
			n_i = n_positions[i]
			ca_i = ca_positions[i]

			if ca_prev is not None and c_prev is not None and n_i is not None and ca_i is not None:
				omega_rad = calc_dihedrals(
					ca_prev.reshape(1, 3),
					c_prev.reshape(1, 3),
					n_i.reshape(1, 3),
					ca_i.reshape(1, 3),
				)[0]
				omega_angles[i] = wrap_angle_deg(float(np.degrees(omega_rad)))

		if i > 0:
			c_prev = c_positions[i - 1]
			n_i = n_positions[i]
			ca_i = ca_positions[i]
			c_i = c_positions[i]

			if c_prev is not None and n_i is not None and ca_i is not None and c_i is not None:
				phi_rad = calc_dihedrals(
					c_prev.reshape(1, 3),
					n_i.reshape(1, 3),
					ca_i.reshape(1, 3),
					c_i.reshape(1, 3),
				)[0]
				phi_angles[i] = wrap_angle_deg(float(np.degrees(phi_rad)))

		if i < n_residues - 1:
			n_i = n_positions[i]
			ca_i = ca_positions[i]
			c_i = c_positions[i]
			n_next = n_positions[i + 1]

			if n_i is not None and ca_i is not None and c_i is not None and n_next is not None:
				psi_rad = calc_dihedrals(
					n_i.reshape(1, 3),
					ca_i.reshape(1, 3),
					c_i.reshape(1, 3),
					n_next.reshape(1, 3),
				)[0]
				psi_angles[i] = wrap_angle_deg(float(np.degrees(psi_rad)))

	# ---------------------------
	# Random split for phi/psi only
	# ---------------------------
	n_detected = int(round(n_residues * (percentage_backbone_torsion_angles / 100.0)))
	n_full_range = n_residues - n_detected

	all_indices = np.arange(n_residues, dtype=int)

	if n_full_range > 0:
		full_range_indices = set(np.random.choice(all_indices, size=n_full_range, replace=False))
	else:
		full_range_indices = set()

	# ---------------------------
	# Write output
	# ---------------------------
	with output_angular_file.open("w", encoding="utf-8") as f:
		f.write("resid\tresname\tomega_center\tomega_radius\tphi_center\tphi_radius\tpsi_center\tpsi_radius\n")

		for i in range(n_residues):
			residue_id = residue_ids[i]
			residue_name = residue_names[i]

			omega_center, omega_radius = sample_angle_interval(omega_angles[i], omega_angle_width)

			if i in full_range_indices:
				phi_center, phi_radius = full_range_angle()
				psi_center, psi_radius = full_range_angle()
			else:
				phi_center, phi_radius = sample_angle_interval(phi_angles[i], phi_angle_width)
				psi_center, psi_radius = sample_angle_interval(psi_angles[i], psi_angle_width)

			f.write(f"{residue_id}\t{residue_name}\t{omega_center:.6f}\t{omega_radius:.6f}\t{phi_center:.6f}\t{phi_radius:.6f}\t{psi_center:.6f}\t{psi_radius:.6f}\n")

	print(f"[OK] TALOS-N-like phi/psi intervals saved to: {output_angular_file}")
# -----------------------------------------------------------------------------------------------------
def wrap_angle_deg(angle: float) -> float:
	"""Wrap an angle in degrees to the interval [-180, 180)."""
	return (angle + 180.0) % 360.0 - 180.0

def circular_distance_deg(angle_a: float, angle_b: float) -> float:
	"""Return the shortest circular distance between two angles in degrees."""
	return abs(wrap_angle_deg(angle_a - angle_b))

def angle_belongs_to_interval(target: float, center: float, radius: float, tol: float = 1e-12) -> bool:
	"""Check whether a target angle belongs to the circular interval centered at center."""
	return circular_distance_deg(target, center) <= radius + tol

def torsion_interval_distance_bounds(
	x1: np.ndarray,
	x2: np.ndarray,
	x3: np.ndarray,
	x4: np.ndarray,
	center_deg: float,
	radius_deg: float,
) -> tuple[float, float]:
	"""
	Return the distance interval induced by the circular torsion interval
	[center - radius, center + radius].

	The distance extrema are evaluated at the interval endpoints and at the
	critical angles 0 and ±180 whenever they belong to the interval.
	"""
	center_deg = wrap_angle_deg(float(center_deg))
	radius_deg = float(radius_deg)

	if radius_deg < 0.0:
		raise ValueError(f"radius_deg must be nonnegative, got {radius_deg}.")

	if radius_deg >= 180.0:
		candidate_angles_deg = [0.0, 180.0]
	else:
		left_deg = wrap_angle_deg(center_deg - radius_deg)
		right_deg = wrap_angle_deg(center_deg + radius_deg)

		candidate_angles_deg = [left_deg, right_deg]

		if angle_belongs_to_interval(0.0, center_deg, radius_deg):
			candidate_angles_deg.append(0.0)

		if angle_belongs_to_interval(180.0, center_deg, radius_deg):
			candidate_angles_deg.append(180.0)

		if angle_belongs_to_interval(-180.0, center_deg, radius_deg):
			candidate_angles_deg.append(-180.0)

	distances = [
		torsion_angle_2_endpoint_distance(x1, x2, x3, x4, np.deg2rad(angle_deg))
		for angle_deg in candidate_angles_deg
	]

	return min(distances), max(distances)

def build_backbone_residue_map(df_atoms: pd.df) -> dict[int, dict[str, object]]:
	"""
	Build a residue-indexed dictionary with N, CA, and C atom ids and coordinates.
	"""
	residue_map: dict[int, dict[str, object]] = {}

	for _, row in df_atoms.iterrows():
		resid = int(row["resid"])
		atom_id = int(row["atom_id"])
		atom_name = str(row["atom_name"]).strip().upper()
		resname = str(row["resname"]).strip().upper()
		coord = np.array([float(row["x"]), float(row["y"]), float(row["z"])], dtype=float)

		if atom_name not in {"N", "CA", "C"}:
			continue

		if resid not in residue_map:
			residue_map[resid] = {
				"resname": resname,
				"atoms": {},
			}

		residue_map[resid]["atoms"][atom_name] = {
			"atom_id": atom_id,
			"coord": coord,
		}

	return residue_map

def get_backbone_atom(
	residue_map: dict[int, dict[str, object]],
	resid: int,
	atom_name: str,
) -> dict[str, object] | None:
	"""Return the stored atom record for a given residue and backbone atom name."""
	residue_data = residue_map.get(int(resid))
	if residue_data is None:
		return None

	atoms = residue_data.get("atoms", {})
	return atoms.get(atom_name)

def omega_interval_to_distance_record(
	residue_map: dict[int, dict[str, object]],
	resid: int,
	resname: str,
	omega_center: float,
	omega_radius: float,
) -> tuple[int, int, int, int, float, float, str, str, str, str] | None:
	"""
	Convert the omega interval of residue i into a CA_{i-1} -- CA_i distance interval.
	"""
	CA_prev = get_backbone_atom(residue_map, resid - 1, "CA")
	C_prev = get_backbone_atom(residue_map, resid - 1, "C")
	N_i = get_backbone_atom(residue_map, resid, "N")
	CA_i = get_backbone_atom(residue_map, resid, "CA")

	if CA_prev is None or C_prev is None or N_i is None or CA_i is None:
		return None

	d_l, d_u = torsion_interval_distance_bounds(
		CA_prev["coord"],
		C_prev["coord"],
		N_i["coord"],
		CA_i["coord"],
		omega_center,
		omega_radius,
	)

	resname_prev = str(residue_map[resid - 1]["resname"]).upper()

	return (
		int(CA_i["atom_id"]),
		int(CA_prev["atom_id"]),
		int(resid),
		int(resid - 1),
		d_l,
		d_u,
		"CA",
		"CA",
		str(resname).upper(),
		resname_prev,
	)

def phi_interval_to_distance_record(
	residue_map: dict[int, dict[str, object]],
	resid: int,
	resname: str,
	phi_center: float,
	phi_radius: float,
) -> tuple[int, int, int, int, float, float, str, str, str, str] | None:
	"""
	Convert the phi interval of residue i into a C_{i-1} -- C_i distance interval.
	"""
	C_prev = get_backbone_atom(residue_map, resid - 1, "C")
	N_i = get_backbone_atom(residue_map, resid, "N")
	CA_i = get_backbone_atom(residue_map, resid, "CA")
	C_i = get_backbone_atom(residue_map, resid, "C")

	if C_prev is None or N_i is None or CA_i is None or C_i is None:
		return None

	d_l, d_u = torsion_interval_distance_bounds(
		C_prev["coord"],
		N_i["coord"],
		CA_i["coord"],
		C_i["coord"],
		phi_center,
		phi_radius,
	)

	resname_prev = str(residue_map[resid - 1]["resname"]).upper()

	return (
		int(C_i["atom_id"]),
		int(C_prev["atom_id"]),
		int(resid),
		int(resid - 1),
		d_l,
		d_u,
		"C",
		"C",
		str(resname).upper(),
		resname_prev,
	)

def psi_interval_to_distance_record(
	residue_map: dict[int, dict[str, object]],
	resid: int,
	resname: str,
	psi_center: float,
	psi_radius: float,
) -> tuple[int, int, int, int, float, float, str, str, str, str] | None:
	"""
	Convert the psi interval of residue i into an N_i -- N_{i+1} distance interval.
	"""
	N_i = get_backbone_atom(residue_map, resid, "N")
	CA_i = get_backbone_atom(residue_map, resid, "CA")
	C_i = get_backbone_atom(residue_map, resid, "C")
	N_next = get_backbone_atom(residue_map, resid + 1, "N")

	if N_i is None or CA_i is None or C_i is None or N_next is None:
		return None

	d_l, d_u = torsion_interval_distance_bounds(
		N_i["coord"],
		CA_i["coord"],
		C_i["coord"],
		N_next["coord"],
		psi_center,
		psi_radius,
	)

	resname_next = str(residue_map[resid + 1]["resname"]).upper()

	return (
		int(N_next["atom_id"]),
		int(N_i["atom_id"]),
		int(resid + 1),
		int(resid),
		d_l,
		d_u,
		"N",
		"N",
		resname_next,
		str(resname).upper(),
	)

def backbone_angular_interval_to_distance_interval(
	tsv_structure_file: str | Path,
	tsv_angular_file: str | Path,
	output_distance_file: str | Path
) -> None:
	"""
	Convert backbone torsion-angle intervals into backbone distance intervals.

	Inputs
	------
	tsv_structure_file
		TSV file containing filtered backbone atoms. Expected columns:
		atom_id, atom_name, resid, resname, x, y, z

	tsv_angular_file
		TSV file containing torsion intervals. Expected columns:
		resid, resname, phi_center, phi_radius, psi_center, psi_radius,
		omega_center, omega_radius

	output_distance_file
		Output TSV file containing the induced distance intervals.

	Generated distance constraints
	------------------------------
	- omega_i -> distance between CA_{i-1} and CA_i
	- phi_i   -> distance between C_{i-1} and C_i
	- psi_i   -> distance between N_i and N_{i+1}

	For each angular interval [center - radius, center + radius], the code
	evaluates the induced distance extrema at the interval endpoints and at
	the critical angles 0 and ±180 whenever they belong to the interval.
	"""
	
	df_atoms = load_filtered_atoms_table(tsv_structure_file)
	
	if not tsv_angular_file.exists():
		raise ValueError(f"Angular-interval file not found: {tsv_angular_file}")
	df_angles = pd.read_csv(tsv_angular_file, sep="\t", dtype={"resname": str})
	df_angles.columns = df_angles.columns.str.strip()

	residue_map = build_backbone_residue_map(df_atoms)

	records: list[tuple[int, int, int, int, float, float, str, str, str, str]] = []

	for _, row in df_angles.iterrows():
		resid = int(row["resid"])
		resname = str(row["resname"]).strip().upper()

		omega_center = float(row["omega_center"])
		omega_radius = float(row["omega_radius"])
		phi_center = float(row["phi_center"])
		phi_radius = float(row["phi_radius"])
		psi_center = float(row["psi_center"])
		psi_radius = float(row["psi_radius"])

		omega_record = omega_interval_to_distance_record(residue_map, resid, resname, omega_center, omega_radius)
		if omega_record is not None:
			records.append(omega_record)

		phi_record = phi_interval_to_distance_record(residue_map, resid, resname, phi_center, phi_radius)
		if phi_record is not None:
			records.append(phi_record)

		psi_record = psi_interval_to_distance_record(residue_map, resid, resname, psi_center, psi_radius)
		if psi_record is not None:
			records.append(psi_record)

	with output_distance_file.open("w", encoding="utf-8") as f:
		f.write("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n")

		for record in records:
			f.write(f"{record[0]}\t{record[1]}\t{record[2]}\t{record[3]}\t{record[4]:.16f}\t{record[5]:.16f}\t{record[6]}\t{record[7]}\t{record[8]}\t{record[9]}\n")

	print(f"[OK] Backbone angular intervals converted to distance intervals saved to: {output_distance_file}")
# -----------------------------------------------------------------------------------------------------
def merge_distance_constraint_files(files_dir: str | Path, instance_file: str | Path) -> Path:
	"""
	Merge all files named 'distance_constraints_*.dat' in output_dir into a single
	consistent interval table.

	Parameters
	----------
	files_dir
		Directory containing intermediate constraint files such as
		distance_constraints_1.dat, distance_constraints_2.dat, ...

	output_instance_file
		Name of the final merged file, saved inside output_dir.

	Behavior
	--------
	- All matching files are concatenated.
	- Intervals referring to the same atom pair are merged.
	- If any constraint is precise (d_l == d_u), the smallest precise value dominates.
	- Otherwise the intersection rule is used:

		d_l = max(d_l)
		d_u = min(d_u)

	A consistency check is performed after merging each group. If the merged
	interval is empty (d_l > d_u), a ValueError is raised.
	"""
	
	input_files = sorted(files_dir.glob("distance_constraints_*.dat"))
	if not input_files:
		raise ValueError(f"No distance_constraints_*.dat files found in: {files_dir}")

	all_entries: list[pd.df] = []

	for input_file in input_files:
		df = pd.read_csv(input_file, sep="\t")
		df.columns = df.columns.str.strip()

		required_columns = [
			"atom_id_i",
			"atom_id_j",
			"resid_i",
			"resid_j",
			"d_l",
			"d_u",
			"atom_name_i",
			"atom_name_j",
			"resname_i",
			"resname_j",
		]
		missing_columns = [column for column in required_columns if column not in df.columns]
		if missing_columns:
			raise ValueError(f"Missing required columns in {input_file}: {', '.join(missing_columns)}")

		df["d_l"] = df["d_l"].astype(float)
		df["d_u"] = df["d_u"].astype(float)

		all_entries.append(df[required_columns])

	merged_df = pd.concat(all_entries, ignore_index=True)

	group_keys = [
		"atom_id_i",
		"atom_id_j",
		"resid_i",
		"resid_j",
		"atom_name_i",
		"atom_name_j",
		"resname_i",
		"resname_j",
	]

	results: list[dict[str, object]] = []

	for key, group in merged_df.groupby(group_keys, dropna=False):
		fixed_intervals = group[np.isclose(group["d_l"], group["d_u"])]

		if not fixed_intervals.empty:
			d_fixed = float(fixed_intervals["d_l"].min())
			d_l = d_fixed
			d_u = d_fixed
		else:
			d_l = float(group["d_l"].max())
			d_u = float(group["d_u"].min())

		if d_l > d_u:
			pair_data = dict(zip(group_keys, key))
			raise ValueError(
				"Inconsistent merged interval for pair "
				f"(atom_id_i={pair_data['atom_id_i']}, atom_id_j={pair_data['atom_id_j']}, "
				f"resid_i={pair_data['resid_i']}, resid_j={pair_data['resid_j']}, "
				f"atom_name_i={pair_data['atom_name_i']}, atom_name_j={pair_data['atom_name_j']}, "
				f"resname_i={pair_data['resname_i']}, resname_j={pair_data['resname_j']}): "
				f"d_l={d_l:.16f} > d_u={d_u:.16f}."
			)

		record = dict(zip(group_keys, key))
		record["d_l"] = d_l
		record["d_u"] = d_u
		results.append(record)

	final_df = pd.DataFrame(results)
	final_df = final_df.sort_values(by=["atom_id_i", "atom_id_j", "resid_i", "resid_j"]).reset_index(drop=True)

	with instance_file.open("w", encoding="utf-8") as f:
		#f.write("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n")
		for _, row in final_df.iterrows():
			f.write(
				f"{int(row['atom_id_i']):5d} {int(row['atom_id_j']):5d} "
				f"{int(row['resid_i']):6d} {int(row['resid_j']):6d} "
				f"{row['d_l']:20.16f} {row['d_u']:20.16f} "
				f"{row['atom_name_i']:>4s} {row['atom_name_j']:>4s} "
				f"{row['resname_i']} {row['resname_j']}\n"
			)
	
	print(f"[OK] Merged distance constraints saved to: {instance_file}")
# -----------------------------------------------------------------------------------------------------


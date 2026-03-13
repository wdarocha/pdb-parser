"""
vertex_ordering.py

Utilities for computing and applying vertex orderings in the distance
geometry pipeline.

This module contains functions used to:

- build vertex orderings from atom/residue information
- sort atom lists
- sort distance constraint files
- support DDGP-compatible orderings

The functions defined here operate only on index lists and ordering
vectors; they do not perform geometry computations.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional
import numpy as np
import pandas as pd

from pdb_parser.utils import *
from pdb_parser.io import *
from pdb_parser.geometry import *

# -----------------------------------------------------------------------------------------------------
def get_missing_backbone_plus_hydrogen_atoms(
	atom_names: set[str],
	resname: str,
	resnum: int,
) -> list[str]:
	"""
	Return the list of missing required atoms for one residue.
	"""
	resname = resname.strip().upper()
	missing_atoms: list[str] = []

	# --- backbone atoms always required
	for atom_name in ("N", "CA", "C"):
		if atom_name not in atom_names:
			missing_atoms.append(atom_name)

	# --- backbone nitrogen hydrogen requirement
	if resname != "PRO":
		if resnum == 1:
			if not any(h in atom_names for h in ("H", "H1", "H2", "H3")):
				missing_atoms.append("H/H1/H2/H3")
		else:
			if "H" not in atom_names:
				missing_atoms.append("H")

	# --- residue-specific side/backbone hydrogen requirements
	if resname == "GLY":
		if ("HA2" not in atom_names) and ("HA3" not in atom_names):
			missing_atoms.append("HA2/HA3")

	elif resname == "PRO":
		if "HA" not in atom_names:
			missing_atoms.append("HA")

		if ("HD2" not in atom_names) and ("HD3" not in atom_names):
			missing_atoms.append("HD2/HD3")

	else:
		if "HA" not in atom_names:
			missing_atoms.append("HA")

	return missing_atoms
# -----------------------------------------------------------------------------------------------------
def validate_backbone_plus_hydrogens_residue(
	df_X_i,
	resnum_i: int,
	pdb_id: str,
	chosen_model: int,
	chosen_chain: str,
) -> int:
	"""
	Validate one residue for the 'backbone_plus_hydrogens' atom selection.

	Returns
	-------
	int
		1 if the structure must be skipped, 0 otherwise.
	"""
	if df_X_i.empty:
		print(
			f"[skipping] {pdb_id}, model {chosen_model}, chain {chosen_chain} "
			f"does not contain residue {resnum_i}"
		)
		return 1

	atom_names_i = set(
		df_X_i[1].astype(str).str.strip().str.upper().tolist()
	)

	resname_i = str(df_X_i.iat[0, 3]).strip().upper()

	missing_atoms = get_missing_backbone_plus_hydrogen_atoms(
		atom_names_i,
		resname_i,
		resnum_i,
	)

	if missing_atoms:
		print(
			f"[skipping] {pdb_id} (model {chosen_model}, chain {chosen_chain}) "
			f"- option 'backbone_plus_hydrogens': missing atom(s) {missing_atoms} "
			f"in residue {resnum_i} ({resname_i})."
		)
		return 1

	return 0
# -----------------------------------------------------------------------------------------------------
def get_nonterminal_h_group(resname: str) -> tuple[str, ...]:
	"""
	Return the allowed atom names for the backbone-N hydrogen position
	in non-terminal residues.
	"""
	resname = resname.strip().upper()

	if resname == "PRO":
		return ("HD2", "HD3")

	return ("H",)

def get_ha_group(resname: str) -> tuple[str, ...]:
	"""
	Return the allowed atom names for the alpha-hydrogen position.
	"""
	resname = resname.strip().upper()

	if resname == "GLY":
		return ("HA2", "HA3")

	return ("HA",)

def get_internal_residue_ordering_template(
	ordering_id: int,
	resname: str,
) -> list[str]:
	"""
	Return the concrete atom ordering template for a non-terminal residue.

	The function resolves hydrogen alternatives (H/HD2/HD3 and HA/HA2/HA3)
	based on the residue name.
	"""

	h_group = get_nonterminal_h_group(resname)
	ha_group = get_ha_group(resname)

	# choose the first valid atom name from each group
	h_atom = h_group[0]
	ha_atom = ha_group[0]

	ordering_templates: dict[int, list[str]] = {
		1:  ["N", "CA", "C", h_atom, ha_atom],
		2:  ["N", "CA", "C", ha_atom, h_atom],
		3:  ["N", "CA", h_atom, ha_atom, "C"],
		4:  ["N", "CA", h_atom, "C", ha_atom],
		5:  ["N", h_atom, "CA", ha_atom, "C"],
		6:  ["N", h_atom, "CA", "C", ha_atom],
		7:  [h_atom, "N", "CA", ha_atom, "C"],
		8:  [h_atom, "N", "CA", "C", ha_atom],
		9:  [h_atom, "CA", "N", ha_atom, "C"],
		10: [h_atom, "CA", "N", "C", ha_atom],
	}

	if ordering_id not in ordering_templates:
		raise ValueError(f"Invalid ordering_id: {ordering_id}. Expected an integer from 1 to 10.")

	return ordering_templates[ordering_id]

def get_numeric_atom_order_from_named_order(
	df_X_i,
	ordered_atom_names: list[str],
) -> list[int]:
	"""
	Return the ordered list of atom ids for one residue from an ordered
	list of atom names.

	The input DataFrame must contain exactly one residue and follow the
	column convention:
	- column 0: atom_id
	- column 1: atom_name
	"""
	if df_X_i.empty:
		raise ValueError("The residue DataFrame is empty.")

	atom_name_to_atom_id = {
		str(row[1]).strip().upper(): int(row[0])
		for _, row in df_X_i.iterrows()
	}

	try:
		return [atom_name_to_atom_id[atom_name] for atom_name in ordered_atom_names]
	except KeyError as exc:
		raise ValueError(
			f"Atom {exc} is missing in residue {int(df_X_i.iat[0, 2])} "
			f"({str(df_X_i.iat[0, 3]).strip().upper()})."
		) from exc

def get_internal_residue_numeric_order(
	df_X_i,
	ordering_id: int,
) -> list[tuple[int, str]]:
	"""
	Return the ordered list of (atom_id, atom_name) pairs for one
	non-terminal residue.
	"""
	if df_X_i.empty:
		raise ValueError("The residue DataFrame is empty.")

	resname = str(df_X_i.iat[0, 3]).strip().upper()

	ordered_atom_names = get_internal_residue_ordering_template(
		ordering_id=ordering_id,
		resname=resname,
	)

	ordered_atom_ids = get_numeric_atom_order_from_named_order(
		df_X_i=df_X_i,
		ordered_atom_names=ordered_atom_names,
	)

	return [
		(atom_id, atom_name)
		for atom_id, atom_name in zip(ordered_atom_ids, ordered_atom_names)
	]
# -----------------------------------------------------------------------------------------------------
def choose_first_existing_atom(
	atom_names_present: set[str],
	candidates: tuple[str, ...],
	excluded_atoms: set[str] | None = None,
) -> str | None:
	"""
	Return the first candidate atom name that exists in the residue and is
	not excluded. Return None if no candidate is available.
	"""
	if excluded_atoms is None:
		excluded_atoms = set()

	for atom_name in candidates:
		if atom_name in atom_names_present and atom_name not in excluded_atoms:
			return atom_name

	return None

def get_first_residue_ha_atom(
	atom_names_present: set[str],
) -> str:
	"""
	Return the alpha-hydrogen atom to be used in the first residue.
	"""
	ha_atom = choose_first_existing_atom(
		atom_names_present,
		("HA", "HA2", "HA3"),
	)

	if ha_atom is None:
		raise ValueError("The first residue does not contain any HA/HA2/HA3 atom.")

	return ha_atom


def first_residue_order(
	df_X_i,
) -> list[tuple[int, str]]:
	"""
	Return the ordered list of (atom_id, atom_name) pairs for the first residue.

	Supported patterns
	------------------
	1. H3, H2, h, N, CA, ha, C
	2. H3, h, H1, N, CA, ha, C
	3. h, H2, H1, N, CA, ha, C
	4. H2, h, N, CA, ha, C
	5. H3, h, N, CA, ha, C
	6. N, CA, C, ha, h

	where h can be H, H1, HD2, or HD3, and ha can be HA, HA2, or HA3.
	"""
	if df_X_i.empty:
		raise ValueError("The first-residue DataFrame is empty.")

	atom_name_to_atom_id = {
		str(row[1]).strip().upper(): int(row[0])
		for _, row in df_X_i.iterrows()
	}

	atom_names_present = set(atom_name_to_atom_id.keys())

	required_backbone = {"N", "CA", "C"}
	missing_backbone = sorted(required_backbone - atom_names_present)

	if missing_backbone:
		raise ValueError(f"The first residue is missing required backbone atom(s): {missing_backbone}.")

	ha_atom = get_first_residue_ha_atom(atom_names_present)

	def atom_pair(atom_name: str) -> tuple[int, str]:
		return (atom_name_to_atom_id[atom_name], atom_name)

	# Candidate atoms that may play the role of h
	h_candidates = ("H", "H1", "HD2", "HD3")

	# --------------------------------------------------------------
	# Pattern 1: H3, H2, h, N, CA, ha, C
	# --------------------------------------------------------------
	if "H3" in atom_names_present and "H2" in atom_names_present:
		h_atom = choose_first_existing_atom(
			atom_names_present,
			h_candidates,
			excluded_atoms={"H3", "H2"},
		)
		if h_atom is not None:
			return [
				atom_pair("H3"),
				atom_pair("H2"),
				atom_pair(h_atom),
				atom_pair("N"),
				atom_pair("CA"),
				atom_pair(ha_atom),
				atom_pair("C"),
			]

	# --------------------------------------------------------------
	# Pattern 2: H3, h, H1, N, CA, ha, C
	# --------------------------------------------------------------
	if "H3" in atom_names_present and "H1" in atom_names_present:
		h_atom = choose_first_existing_atom(
			atom_names_present,
			h_candidates,
			excluded_atoms={"H3", "H1"},
		)
		if h_atom is not None:
			return [
				atom_pair("H3"),
				atom_pair(h_atom),
				atom_pair("H1"),
				atom_pair("N"),
				atom_pair("CA"),
				atom_pair(ha_atom),
				atom_pair("C"),
			]

	# --------------------------------------------------------------
	# Pattern 3: h, H2, H1, N, CA, ha, C
	# --------------------------------------------------------------
	if "H2" in atom_names_present and "H1" in atom_names_present:
		h_atom = choose_first_existing_atom(
			atom_names_present,
			h_candidates,
			excluded_atoms={"H2", "H1"},
		)
		if h_atom is not None:
			return [
				atom_pair(h_atom),
				atom_pair("H2"),
				atom_pair("H1"),
				atom_pair("N"),
				atom_pair("CA"),
				atom_pair(ha_atom),
				atom_pair("C"),
			]

	# --------------------------------------------------------------
	# Pattern 4: H2, h, N, CA, ha, C
	# --------------------------------------------------------------
	if "H2" in atom_names_present:
		h_atom = choose_first_existing_atom(
			atom_names_present,
			h_candidates,
			excluded_atoms={"H2"},
		)
		if h_atom is not None:
			return [
				atom_pair("H2"),
				atom_pair(h_atom),
				atom_pair("N"),
				atom_pair("CA"),
				atom_pair(ha_atom),
				atom_pair("C"),
			]

	# --------------------------------------------------------------
	# Pattern 5: H3, h, N, CA, ha, C
	# --------------------------------------------------------------
	if "H3" in atom_names_present:
		h_atom = choose_first_existing_atom(
			atom_names_present,
			h_candidates,
			excluded_atoms={"H3"},
		)
		if h_atom is not None:
			return [
				atom_pair("H3"),
				atom_pair(h_atom),
				atom_pair("N"),
				atom_pair("CA"),
				atom_pair(ha_atom),
				atom_pair("C"),
			]

	# --------------------------------------------------------------
	# Pattern 6: N, CA, C, ha, h
	# Covers the remaining single-hydrogen cases.
	# --------------------------------------------------------------
	h_atom = choose_first_existing_atom(
		atom_names_present,
		h_candidates,
	)
	if h_atom is not None:
		return [
			atom_pair("N"),
			atom_pair("CA"),
			atom_pair("C"),
			atom_pair(ha_atom),
			atom_pair(h_atom),
		]

	raise ValueError("Could not match any valid ordering pattern for the first residue.")
# -----------------------------------------------------------------------------------------------------
def sort_structure_dataframe(
	df,
	new_order: list[int],
):
	df_sorted = df.set_index(0).loc[new_order].reset_index()

	df_sorted[0] = range(1, len(df_sorted) + 1)

	return df_sorted
# -----------------------------------------------------------------------------------------------------
def sort_distance_dataframe(
	df_D,
	atom_ids_new_order: list[int],
):
	"""
	Sort the atom indices in a distance-constraint DataFrame.

	Parameters
	----------
	df_D
	    Distance constraint DataFrame.

	atom_ids_new_order
	    Ordered list of old atom ids.
	    The position in this list defines the new atom id (1-based).

	Column convention
	-----------------
	0 : atom_id_i
	1 : atom_id_j
	2 : resid_i
	3 : resid_j
	4 : d_l
	5 : d_u
	6 : atom_name_i
	7 : atom_name_j
	8 : resname_i
	9 : resname_j
	"""

	df_D_new = df_D.copy()

	# --------------------------------------------------------------
	# Build old_atom_id -> new_atom_id lookup
	# --------------------------------------------------------------
	old_to_new = {
		old_atom_id: new_atom_id
		for new_atom_id, old_atom_id in enumerate(atom_ids_new_order, start=1)
	}

	# --------------------------------------------------------------
	# Apply remapping to the first two columns
	# --------------------------------------------------------------
	df_D_new[0] = df_D_new[0].map(old_to_new)
	df_D_new[1] = df_D_new[1].map(old_to_new)

	# Ensure integers
	df_D_new[0] = df_D_new[0].astype(int)
	df_D_new[1] = df_D_new[1].astype(int)

	# --------------------------------------------------------------
	# Enforce atom_id_i >= atom_id_j
	# --------------------------------------------------------------
	swap_mask = df_D_new[0] < df_D_new[1]

	if swap_mask.any():

		# swap atom indices
		df_D_new.loc[swap_mask, [0, 1]] = df_D_new.loc[swap_mask, [1, 0]].to_numpy()

		# swap residue indices
		df_D_new.loc[swap_mask, [2, 3]] = df_D_new.loc[swap_mask, [3, 2]].to_numpy()

		# swap atom names
		df_D_new.loc[swap_mask, [6, 7]] = df_D_new.loc[swap_mask, [7, 6]].to_numpy()

		# swap residue names
		df_D_new.loc[swap_mask, [8, 9]] = df_D_new.loc[swap_mask, [9, 8]].to_numpy()

	# --------------------------------------------------------------
	# Final ordering
	# --------------------------------------------------------------
	df_D_new = df_D_new.sort_values(by=[0, 1]).reset_index(drop=True)

	return df_D_new
# -----------------------------------------------------------------------------------------------------
CLIQUE_PATTERN_MAP: dict[int, tuple[tuple[int, str], ...]] = {
	1: ((-1, "N"), (-1, "CA"), (-1, "C"), (0, "N")),
	2: ((-1, "CA"), (-1, "C"), (0, "N"), (0, "CA")),
	3: ((-1, "C"), (0, "N"), (0, "CA"), (0, "C")),
	4: ((-1, "C"), (0, "N"), (0, "CA"), (0, "HN")),
	5: ((0, "N"), (0, "CA"), (0, "C"), (0, "HA")),
	6: ((0, "HN"), (0, "N"), (0, "CA"), (0, "HA")),
	7: ((0, "N"), (0, "CA"), (0, "HA"), (0, "C")),
	8: ((-1, "CA"), (-1, "C"), (0, "N"), (0, "HN")),
	9: ((-1, "C"), (0, "N"), (0, "HN"), (0, "CA")),
	10: ((-1, "HA"), (-1, "CA"), (-1, "C"), (0, "HN")),
	11: ((-1, "CA"), (-1, "C"), (0, "HN"), (0, "N")),
	12: ((-1, "CA"), (-1, "C"), (0, "HN"), (0, "CA")),
	13: ((-1, "C"), (0, "HN"), (0, "CA"), (0, "N")),
}

DDGP_ORDER_PATTERN_MAP: dict[int, tuple[int, int, int, int, int]] = {
	1: (1, 2, 3, 4, 5),
	2: (1, 2, 3, 5, 4),
	3: (1, 2, 4, 6, 7),
	4: (1, 2, 4, 3, 5),
	5: (1, 8, 9, 6, 7),
	6: (1, 8, 9, 3, 5),
	7: (10, 11, 9, 6, 7),
	8: (10, 11, 9, 3, 5),
	9: (10, 12, 13, 6, 7),
	10: (10, 12, 13, 3, 5),
}

ATOM_ALIAS_MAP: dict[str, tuple[str, ...]] = {
	"HN": ("H1", "H", "HD2", "HD3"),
	"HA": ("HA", "HA2", "HA3"),
}


def build_available_atoms(df_Xreord) -> list[tuple[int, str]]:
	"""
	Build the available atom list as (residue_id, atom_name).

	Parameters
	----------
	df_Xreord : pd.DataFrame
		Reordered atom dataframe.

	Returns
	-------
	list[tuple[int, str]]
		List of (residue_id, atom_name).
	"""
	available_atoms: list[tuple[int, str]] = []

	for row in df_Xreord.itertuples(index=False):
		atom_name = str(row[1])
		residue_id = int(row[2])
		available_atoms.append((residue_id, atom_name))

	return available_atoms


def get_atoms_of_residue(
	residue_id: int,
	available_atoms: list[tuple[int, str]],
) -> set[str]:
	"""
	Return the atom names that belong to one residue.

	Parameters
	----------
	residue_id : int
		Target residue ID.
	available_atoms : list[tuple[int, str]]
		List of (residue_id, atom_name).

	Returns
	-------
	set[str]
		Set of atom names present in the residue.
	"""
	return {
		atom_name
		for current_residue_id, atom_name in available_atoms
		if current_residue_id == residue_id
	}


def resolve_atom_name(
	atom_name: str,
	residue_id: int,
	available_atoms: list[tuple[int, str]] | None = None,
) -> str:
	"""
	Resolve atom aliases according to the atoms available in the target residue.

	HN -> H1, H, HD2, HD3
	HA -> HA, HA2, HA3

	Parameters
	----------
	atom_name : str
		Base atom label from the clique pattern.
	residue_id : int
		Residue ID where the atom must be resolved.
	available_atoms : list[tuple[int, str]] | None
		List of (residue_id, atom_name).

	Returns
	-------
	str
		Resolved atom name.
	"""
	if available_atoms is None:
		return atom_name

	candidate_atoms = ATOM_ALIAS_MAP.get(atom_name)
	if candidate_atoms is None:
		return atom_name

	residue_atoms = get_atoms_of_residue(residue_id, available_atoms)

	for candidate_atom in candidate_atoms:
		if candidate_atom in residue_atoms:
			return candidate_atom

	return atom_name


def get_clique_pattern_ids(ddgp_pattern_id: int) -> tuple[int, int, int, int, int]:
	"""
	Return the 5 clique pattern IDs associated with a DDGP order pattern.
	"""
	if ddgp_pattern_id not in DDGP_ORDER_PATTERN_MAP:
		raise ValueError(
			f"DDGP order pattern {ddgp_pattern_id} is not implemented. "
			f"Valid values are 1 to 10."
		)

	return DDGP_ORDER_PATTERN_MAP[ddgp_pattern_id]


def map_clique_type(pattern_id: int) -> int:
	"""
	Map a clique pattern ID to its corresponding class.

	Parameters
	----------
	pattern_id : int
		Clique pattern ID.

	Returns
	-------
	int
		Mapped clique type.
	"""
	if pattern_id in {2, 4, 5, 7, 8, 9, 11, 12, 13}:
		return 1

	if pattern_id in {6, 10}:
		return 2

	if pattern_id in {1, 3}:
		return 3

	raise ValueError(f"Invalid pattern id: {pattern_id}")


def build_clique_pairs(
	residue_id: int,
	clique_pattern_id: int,
	available_atoms: list[tuple[int, str]] | None = None,
) -> list[tuple[int, str]]:
	"""
	Build the 4 ordered pairs (residue_id, atom_name) for one clique pattern.

	Parameters
	----------
	residue_id : int
		Reference residue index i.
	clique_pattern_id : int
		Clique pattern selector.
	available_atoms : list[tuple[int, str]] | None
		List of (residue_id, atom_name).

	Returns
	-------
	list[tuple[int, str]]
		List with 4 ordered pairs.
	"""
	if clique_pattern_id not in CLIQUE_PATTERN_MAP:
		raise ValueError(
			f"Clique pattern {clique_pattern_id} is not implemented. "
			f"Valid values are 1 to 13."
		)

	pairs: list[tuple[int, str]] = []

	for residue_offset, atom_name in CLIQUE_PATTERN_MAP[clique_pattern_id]:
		current_residue_id = residue_id + residue_offset
		resolved_atom_name = resolve_atom_name(
			atom_name=atom_name,
			residue_id=current_residue_id,
			available_atoms=available_atoms,
		)
		pairs.append((current_residue_id, resolved_atom_name))

	return pairs


def build_ddgp_pattern_entries(
	ddgp_pattern_id: int,
	residue_id: int,
	available_atoms: list[tuple[int, str]] | None = None,
) -> list[tuple[list[tuple[int, str]], int]]:
	"""
	Build the 5 entries associated with one DDGP pattern.

	Parameters
	----------
	ddgp_pattern_id : int
		DDGP order pattern selector.
	residue_id : int
		Reference residue index i.
	available_atoms : list[tuple[int, str]] | None
		List of (residue_id, atom_name).

	Returns
	-------
	list[tuple[list[tuple[int, str]], int]]
		List with 5 entries in the form (pairs, clique_type).
	"""
	clique_pattern_ids = get_clique_pattern_ids(ddgp_pattern_id)

	entries: list[tuple[list[tuple[int, str]], int]] = []

	for clique_pattern_id in clique_pattern_ids:
		pairs = build_clique_pairs(
			residue_id=residue_id,
			clique_pattern_id=clique_pattern_id,
			available_atoms=available_atoms,
		)
		clique_type = map_clique_type(clique_pattern_id)
		entries.append((pairs, clique_type))

	return entries
# -----------------------------------------------------------------------------------------------------
ResidueAtomPair = tuple[int, str]
Residue1Row = tuple[list[Optional[ResidueAtomPair]], int]

HN_CANDIDATES: tuple[str, ...] = ("H1", "H", "HD2", "HD3")
HA_CANDIDATES: tuple[str, ...] = ("HA", "HA2", "HA3")
HM_CANDIDATES: tuple[str, ...] = ("H3", "H2")

def _get_first_residue_atom_names(
	available_atoms: list[tuple[int, str]],
) -> set[str]:
	"""
	Extract the atom names that belong to residue 1.

	Parameters
	----------
	available_atoms : list[tuple[int, str]]
		List of (residue_id, atom_name).

	Returns
	-------
	set[str]
		Set of atom names found in residue 1.
	"""
	first_residue_atom_names: set[str] = set()

	for residue_id, atom_name in available_atoms:
		if residue_id == 1:
			first_residue_atom_names.add(atom_name)

	return first_residue_atom_names

def _pick_first_available_atom(
	available_atom_names: set[str],
	candidates: tuple[str, ...],
	label: str,
) -> str:
	"""
	Return the first atom name from candidates that exists in available_atom_names.

	Parameters
	----------
	available_atom_names : set[str]
		Set of atom names available in the residue.
	candidates : tuple[str, ...]
		Priority-ordered candidate atom names.
	label : str
		Logical atom label used only for error messages.

	Returns
	-------
	str
		Resolved atom name.

	Raises
	------
	ValueError
		If none of the candidate atom names exists.
	"""
	for atom_name in candidates:
		if atom_name in available_atom_names:
			return atom_name

	raise ValueError(
		f"Could not resolve {label}. Expected one of {candidates}, "
		f"but none was found in {sorted(available_atom_names)}."
	)

def _build_pair(atom_name: str) -> ResidueAtomPair:
	"""
	Build the pair for residue 1.

	Parameters
	----------
	atom_name : str
		Atom name.

	Returns
	-------
	tuple[int, str]
		Pair in the form (1, atom_name).
	"""
	return (1, atom_name)

def _build_row(
	atom_names: list[Optional[str]],
	row_id: int,
) -> Residue1Row:
	"""
	Build one row with 4 fixed positions plus the row ID.

	Parameters
	----------
	atom_names : list[Optional[str]]
		List with exactly 4 atom names or None values.
	row_id : int
		Integer associated with the row.

	Returns
	-------
	Residue1Row
		Row in the form ([pair_or_none, ..., pair_or_none], row_id).
	"""
	if len(atom_names) != 4:
		raise ValueError(f"atom_names must have length 4, got {len(atom_names)}.")

	pairs: list[Optional[ResidueAtomPair]] = []

	for atom_name in atom_names:
		if atom_name is None:
			pairs.append(None)
		else:
			pairs.append(_build_pair(atom_name))

	return pairs, row_id

def build_first_residue_pattern(
	available_atoms: list[tuple[int, str]],
) -> list[Residue1Row]:
	"""
	Build the initialization pattern for residue 1 according to the number
	of atoms present in residue 1.

	Supported cases
	---------------
	- 7 atoms
	- 6 atoms
	- 5 atoms

	Resolution rules
	----------------
	- HN is resolved from: H1, H, HD2, HD3
	- HA is resolved from: HA, HA2, HA3
	- HM is resolved from: H3, H2

	Parameters
	----------
	available_atoms : list[tuple[int, str]]
		List of (residue_id, atom_name).

	Returns
	-------
	list[Residue1Row]
		List of rows. Each row is:
		([pair_or_none, pair_or_none, pair_or_none, pair_or_none], int)

	Notes
	-----
	- Missing positions are represented with None.
	- All returned pairs use residue index 1.
	"""
	available_atom_names = _get_first_residue_atom_names(available_atoms)
	atom_count = len(available_atom_names)

	if atom_count == 7:
		hn_atom = _pick_first_available_atom(available_atom_names, HN_CANDIDATES, "HN")
		ha_atom = _pick_first_available_atom(available_atom_names, HA_CANDIDATES, "HA")

		return [
			_build_row([None, None, None, "H3"], 1),
			_build_row([None, None, "H3", "H2"], 1),
			_build_row([None, "H3", "H2", hn_atom], 1),
			_build_row(["H3", "H2", hn_atom, "N"], 1),
			_build_row(["H2", hn_atom, "N", "CA"], 1),
			_build_row([hn_atom, "N", "CA", ha_atom], 2),
			_build_row(["N", "CA", ha_atom, "C"], 1),
		]

	if atom_count == 6:
		hm_atom = _pick_first_available_atom(available_atom_names, HM_CANDIDATES, "HM")
		hn_atom = _pick_first_available_atom(available_atom_names, HN_CANDIDATES, "HN")
		ha_atom = _pick_first_available_atom(available_atom_names, HA_CANDIDATES, "HA")

		return [
			_build_row([None, None, None, hm_atom], 1),
			_build_row([None, None, hm_atom, hn_atom], 1),
			_build_row([None, hm_atom, hn_atom, "N"], 1),
			_build_row([hm_atom, hn_atom, "N", "CA"], 1),
			_build_row([hn_atom, "N", "CA", ha_atom], 2),
			_build_row(["N", "CA", ha_atom, "C"], 1),
		]

	if atom_count == 5:
		ha_atom = _pick_first_available_atom(available_atom_names, HA_CANDIDATES, "HA")
		hn_atom = _pick_first_available_atom(available_atom_names, HN_CANDIDATES, "HN")

		return [
			_build_row([None, None, None, "N"], 1),
			_build_row([None, None,"N", "CA",], 1),
			_build_row([None, "N", "CA", "C"], 1),
			_build_row(["N", "CA", "C", ha_atom], 1),
			_build_row([ha_atom, "CA", "N", hn_atom], 2),
		]

	raise ValueError(
		f"Unsupported number of atoms for residue 1: {atom_count}. "
		f"Expected 5, 6, or 7. Found atoms: {sorted(available_atom_names)}."
	)
# -----------------------------------------------------------------------------------------------------
def build_interval_lookup(
	df_I: pd.DataFrame,
) -> dict[tuple[int, int], tuple[float, float]]:
	"""
	Build a lookup table for interval distance constraints.

	The lookup maps an atom pair `(atom_1_idx, atom_4_idx)` to the
	corresponding lower and upper bounds `(dl, du)`.

	Expected columns in df_I
	------------------------
	column 1 : first atom index
	column 2 : fourth atom index
	column 5 : lower bound dl
	column 6 : upper bound du

	Parameters
	----------
	df_I : pd.DataFrame
		Dataframe containing interval distance constraints.

	Returns
	-------
	dict[tuple[int, int], tuple[float, float]]
		Dictionary mapping `(atom_1_idx, atom_4_idx)` to `(dl, du)`.

	Notes
	-----
	If the same `(atom_1_idx, atom_4_idx)` pair appears multiple times,
	the last occurrence in `df_I` overwrites previous ones.
	"""
	lookup: dict[tuple[int, int], tuple[float, float]] = {}

	for row in df_I.itertuples(index=False):
		atom_1_idx = int(row[0])
		atom_4_idx = int(row[1])
		dl = float(row[4])
		du = float(row[5])

		lookup[(atom_1_idx, atom_4_idx)] = (dl, du)

	return lookup

def get_interval_dl_du(
	atom_1_idx: int,
	atom_4_idx: int,
	interval_lookup: dict[tuple[int, int], tuple[float, float]],
) -> tuple[float, float]:
	"""
	Return the interval bounds `(dl, du)` for a given atom pair.

	Parameters
	----------
	atom_1_idx : int
		Index of the first atom.
	atom_4_idx : int
		Index of the fourth atom.
	interval_lookup : dict[tuple[int, int], tuple[float, float]]
		Lookup table created by `build_interval_lookup`.

	Returns
	-------
	tuple[float, float]
		Lower and upper interval bounds `(dl, du)`.

	Raises
	------
	KeyError
		If the atom pair `(atom_1_idx, atom_4_idx)` is not present in the lookup.
	"""
	key = (atom_1_idx, atom_4_idx)

	if key not in interval_lookup:
		raise KeyError(
			f"Interval not found for atom pair ({atom_1_idx}, {atom_4_idx})."
		)

	return interval_lookup[key]

def update_torsion_matrix(
	T: np.ndarray,
	df_X: pd.DataFrame,
	df_A,
	df_I: pd.DataFrame,
) -> np.ndarray:
	"""
	Update the torsion-related columns of matrix `T`.

	The matrix `T` is assumed to store, in columns 0 to 3, four 1-based
	atom indices defining a local geometric pattern. Column 4 initially
	stores the clique type, and columns 5 and 6 store torsion-related
	information.

	Processing rules by clique type
	-------------------------------
	- clique_type == 1:
		Compute the torsion angle directly from atomic coordinates.
		Then store:
			- T[k, 4] = sign of tau
			- T[k, 5] = absolute value of tau in degrees

	- clique_type == 2:
		Use an interval constraint on distance d14 from `df_I`, together with
		the other exact distances derived from coordinates, to estimate a
		torsion interval. Then store:
			- T[k, 4] = 0.0
			- T[k, 5] = midpoint of the absolute torsion interval in degrees
			- T[k, 6] = half-width of the absolute torsion interval in degrees

	- clique_type == 3:
		Read tabulated phi/psi values from `df_A`.
		If atom 1 is a carbon atom "C", use phi from the current residue row.
		Otherwise, use psi from the previous row.
		Then store:
			- T[k, 4] = sign of phi/psi
			- T[k, 5] = absolute value of phi/psi
			- T[k, 6] = corresponding uncertainty

	Assumptions
	-----------
	- Atom indices in T[:, 0:4] are 1-based.
	- df_X columns:
			0 -> atom index
			1 -> atom name
			2 -> residue id
			4:7 -> Cartesian coordinates x, y, z
	- df_A column 0 stores residue id.
	- df_A columns 4 and 5 store phi and delta_phi.
	- df_A columns 6 and 7 store psi and delta_psi.
	- For psi access in clique type 3, the relevant row is assumed to be
		the previous row in df_A.

	Parameters
	----------
	T : np.ndarray
		Input/output matrix containing atom quadruplets and torsion metadata.
	df_X : pd.DataFrame
		Dataframe containing atom identifiers, names, residue indices, and coordinates.
	df_A : np.ndarray | pd.DataFrame
		Table containing tabulated phi/psi values and uncertainties.
	df_I : pd.DataFrame
		Dataframe containing interval distance constraints.

	Returns
	-------
	np.ndarray
		The updated matrix `T`.

	Raises
	------
	ValueError
		If an atom index or residue id cannot be found in the corresponding lookup,
		or if an unsupported clique type is encountered.
	KeyError
		If an interval constraint is missing for a clique of type 2.
	"""
	
	X = df_X.iloc[:, 4:7].to_numpy()

	if hasattr(df_A, "to_numpy"):
		A = df_A.to_numpy()
	else:
		A = df_A

	# Build lookups once
	atom_to_name: dict[int, str] = {}
	atom_to_residue: dict[int, int] = {}
		
	for row in df_X.itertuples(index=False):
		atom_index = int(row[0])
		atom_name = str(row[1])
		residue_id = int(row[2])

		atom_to_name[atom_index] = atom_name
		atom_to_residue[atom_index] = residue_id
	
	residue_to_Arow = {int(row[0]): i for i, row in enumerate(A)}
	
	interval_lookup = build_interval_lookup(df_I)

	n_rows = T.shape[0]
	
	for k in range(3, n_rows):
		clique_type = int(T[k, 4])

		if clique_type == 1:
			atom_1_idx = int(T[k, 0])
			atom_2_idx = int(T[k, 1])
			atom_3_idx = int(T[k, 2])
			atom_4_idx = int(T[k, 3])

			tau = torsion_angle_with_points(
				X[atom_4_idx - 1],
				X[atom_3_idx - 1],
				X[atom_2_idx - 1],
				X[atom_1_idx - 1],
			)

			T[k, 4] = np.sign(tau)
			T[k, 5] = np.abs(tau) * 180.0 / np.pi

		elif clique_type == 2:
			T[k, 4] = 0.0

			atom_1_idx = int(T[k, 0])
			atom_2_idx = int(T[k, 1])
			atom_3_idx = int(T[k, 2])
			atom_4_idx = int(T[k, 3])

			d14l, d14u = get_interval_dl_du(atom_1_idx, atom_4_idx, interval_lookup)

			x4, x3, x2, x1 = X[
				[atom_1_idx - 1, atom_2_idx - 1, atom_3_idx - 1, atom_4_idx - 1]
			]

			d12 = np.linalg.norm(x1 - x2)
			d13 = np.linalg.norm(x1 - x3)
			d23 = np.linalg.norm(x2 - x3)
			d24 = np.linalg.norm(x2 - x4)
			d34 = np.linalg.norm(x3 - x4)

			abs_tau_l = distances_2_abs_torsion_angle(d12, d13, d14l, d23, d24, d34) * 180.0 / np.pi
			abs_tau_u = distances_2_abs_torsion_angle(d12, d13, d14u, d23, d24, d34) * 180.0 / np.pi
			
			abs_tau_mid = 0.5 * (abs_tau_u + abs_tau_l)
			delta_tau   = 0.5 * (abs_tau_u - abs_tau_l)

			T[k, 5] = abs_tau_mid 
			T[k, 6] = delta_tau

		elif clique_type == 3:
			atom_index = int(T[k, 0])

			if atom_index not in atom_to_name:
				raise ValueError(f"Atom index {atom_index} not found in df_X.")

			if atom_index not in atom_to_residue:
				raise ValueError(f"Residue for atom index {atom_index} not found in df_X.")

			atom_name = atom_to_name[atom_index]
			residue_id = atom_to_residue[atom_index]

			if residue_id not in residue_to_Arow:
				raise ValueError(f"Residue {residue_id} not found in df_A.")

			row_idx = residue_to_Arow[residue_id]
			is_phi = (atom_name == "C")

			if is_phi:
				phi = A[row_idx, 4]
				delta_phi = A[row_idx, 5]

				T[k, 4] = np.sign(phi)
				T[k, 5] = np.abs(phi)
				T[k, 6] = delta_phi

			else:
				psi = A[row_idx - 1, 6]
				delta_psi = A[row_idx - 1, 7]

				T[k, 4] = np.sign(psi)
				T[k, 5] = np.abs(psi)
				T[k, 6] = delta_psi

		else:
			raise ValueError(
				f"Unsupported clique type {clique_type} at row {k}. "
				f"Expected 1, 2 or 3."
			)

	return T
# -----------------------------------------------------------------------------------------------------
def build_atom_lookup(df_X: pd.DataFrame) -> dict[tuple[int, str], int]:
	"""
	Build a lookup table from `(residue_id, atom_name)` to atom index.

	The dataframe `df_X` is assumed to contain, at minimum, the atom index,
	atom name, and residue identifier. This function extracts those three
	fields and creates a dictionary that allows constant-time lookup of the
	atom index associated with a given `(residue_id, atom_name)` pair.

	Expected columns in df_X
	------------------------
	column 0 : atom index
	column 1 : atom name
	column 2 : residue id

	Parameters
	----------
	df_X : pd.DataFrame
		Dataframe containing atomic information.

	Returns
	-------
	dict[tuple[int, str], int]
		Dictionary mapping `(residue_id, atom_name)` to the corresponding
		1-based atom index stored in `df_X`.

	Notes
	-----
	If the same `(residue_id, atom_name)` pair appears multiple times,
	the last occurrence overwrites previous ones.
	"""
	atom_indices = df_X.iloc[:, 0].to_numpy()
	atom_names = df_X.iloc[:, 1].to_numpy()
	residue_ids = df_X.iloc[:, 2].to_numpy()

	atom_lookup: dict[tuple[int, str], int] = {}

	for atom_index, atom_name, residue_id in zip(atom_indices, atom_names, residue_ids):
		atom_lookup[(int(residue_id), str(atom_name))] = int(atom_index)

	return atom_lookup

def build_atom_clique_index_matrix(
	atom_cliques: list[tuple[list[tuple[int, str] | None], int]],
	df_X: pd.DataFrame,
	df_A: pd.DataFrame,
	df_I: pd.DataFrame,
) -> np.ndarray:
	"""
	Convert atom cliques into the numerical clique matrix `T`.

	Each element of `atom_cliques` has the form:

		([pair_1, pair_2, pair_3, pair_4], clique_type)

	where each `pair_j` is either:
	- `(residue_id, atom_name)`, or
	- `None` for missing positions.

	This function performs two steps:
	1. Convert each `(residue_id, atom_name)` pair into the corresponding
	   atom index using `df_X`.
	2. Fill the geometric/torsion-related columns by calling
	   `update_torsion_matrix`.

	Final output column order
	-------------------------
	column 0 : atom index of pair 4
	column 1 : atom index of pair 3
	column 2 : atom index of pair 2
	column 3 : atom index of pair 1
	column 4 : sign or clique-type-derived torsion metadata
	column 5 : central torsion value
	column 6 : torsion deviation / uncertainty

	Parameters
	----------
	atom_cliques : list[tuple[list[tuple[int, str] | None], int]]
		List of clique definitions. Each row contains exactly 4 atom pairs
		(or `None`) and one clique type.
	df_X : pd.DataFrame
		Atom dataframe.
	df_A : pd.DataFrame
		Torsion-angle constraint dataframe.
	df_I : pd.DataFrame
		Distance-interval constraint dataframe.

	Returns
	-------
	np.ndarray
		Matrix `T` with shape `(len(atom_cliques), 7)`.

	Notes
	-----
	- Columns 0 to 3 store atom indices, but the matrix dtype is floating-point
	  because columns 4 to 6 store real-valued angular information.
	- Missing pairs are encoded as zero in columns 0 to 3.
	- Atom indices remain 1-based, consistent with `df_X`.

	Raises
	------
	ValueError
		If a clique row does not contain exactly 4 positions.
	KeyError
		If a `(residue_id, atom_name)` pair is not found in `df_X`.
	"""
	atom_lookup = build_atom_lookup(df_X)

	n_rows = len(atom_cliques)
	T = np.zeros((n_rows, 7), dtype=np.float64)

	lookup_get = atom_lookup.get

	for row_id, (pairs, clique_type) in enumerate(atom_cliques):
		if len(pairs) != 4:
			raise ValueError(
				f"Each clique must contain exactly 4 positions, got {len(pairs)} at row {row_id}."
			)

		pair_1 = pairs[0]
		pair_2 = pairs[1]
		pair_3 = pairs[2]
		pair_4 = pairs[3]

		if pair_4 is not None:
			key_4 = (int(pair_4[0]), str(pair_4[1]))
			value_4 = lookup_get(key_4)
			if value_4 is None:
				raise KeyError(f"Atom pair {key_4} was not found in df_X.")
			T[row_id, 0] = value_4

		if pair_3 is not None:
			key_3 = (int(pair_3[0]), str(pair_3[1]))
			value_3 = lookup_get(key_3)
			if value_3 is None:
				raise KeyError(f"Atom pair {key_3} was not found in df_X.")
			T[row_id, 1] = value_3

		if pair_2 is not None:
			key_2 = (int(pair_2[0]), str(pair_2[1]))
			value_2 = lookup_get(key_2)
			if value_2 is None:
				raise KeyError(f"Atom pair {key_2} was not found in df_X.")
			T[row_id, 2] = value_2

		if pair_1 is not None:
			key_1 = (int(pair_1[0]), str(pair_1[1]))
			value_1 = lookup_get(key_1)
			if value_1 is None:
				raise KeyError(f"Atom pair {key_1} was not found in df_X.")
			T[row_id, 3] = value_1

		T[row_id, 4] = float(clique_type)

	T = update_torsion_matrix(T, df_X, df_A, df_I)

	return T
# -----------------------------------------------------------------------------------------------------


"""
vertex_ordering.py

Utilities for computing and applying vertex orderings in the distance
geometry pipeline.

This module contains functions used to:

- build vertex orderings from atom/residue information
- reorder atom lists
- reorder distance constraint files
- support DDGP-compatible orderings

The functions defined here operate only on index lists and ordering
vectors; they do not perform geometry computations.
"""

from pathlib import Path

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
) -> list[int]:
	"""
	Return the ordered list of atom ids for one non-terminal residue.
	"""
	
	resname = str(df_X_i.iat[0, 3]).strip().upper()

	ordered_atom_names = get_internal_residue_ordering_template(
		ordering_id=ordering_id,
		resname=resname,
	)

	return get_numeric_atom_order_from_named_order(
		df_X_i=df_X_i,
		ordered_atom_names=ordered_atom_names,
	)
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
) -> list[int]:
	"""
	Return the ordered list of atom ids for the first residue.

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
		raise ValueError(
			f"The first residue is missing required backbone atom(s): {missing_backbone}."
		)

	ha_atom = get_first_residue_ha_atom(atom_names_present)

	def atom_id(atom_name: str) -> int:
		return atom_name_to_atom_id[atom_name]

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
				atom_id("H3"),
				atom_id("H2"),
				atom_id(h_atom),
				atom_id("N"),
				atom_id("CA"),
				atom_id(ha_atom),
				atom_id("C"),
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
				atom_id("H3"),
				atom_id(h_atom),
				atom_id("H1"),
				atom_id("N"),
				atom_id("CA"),
				atom_id(ha_atom),
				atom_id("C"),
			]

	# --------------------------------------------------------------
	# Pattern 3: h, H2, H1, N, CA, ha, C
	# This now correctly covers cases like H + H2 + H1.
	# --------------------------------------------------------------
	if "H2" in atom_names_present and "H1" in atom_names_present:
		h_atom = choose_first_existing_atom(
			atom_names_present,
			h_candidates,
			excluded_atoms={"H2", "H1"},
		)
		if h_atom is not None:
			return [
				atom_id(h_atom),
				atom_id("H2"),
				atom_id("H1"),
				atom_id("N"),
				atom_id("CA"),
				atom_id(ha_atom),
				atom_id("C"),
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
				atom_id("H2"),
				atom_id(h_atom),
				atom_id("N"),
				atom_id("CA"),
				atom_id(ha_atom),
				atom_id("C"),
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
				atom_id("H3"),
				atom_id(h_atom),
				atom_id("N"),
				atom_id("CA"),
				atom_id(ha_atom),
				atom_id("C"),
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
			atom_id("N"),
			atom_id("CA"),
			atom_id("C"),
			atom_id(ha_atom),
			atom_id(h_atom),
		]

	raise ValueError("Could not match any valid ordering pattern for the first residue.")
# -----------------------------------------------------------------------------------------------------
def reorder_and_renumber_structure_dataframe(
	df,
	new_order: list[int],
):
	df_reordered = df.set_index(0).loc[new_order].reset_index()

	df_reordered[0] = range(1, len(df_reordered) + 1)

	return df_reordered
# -----------------------------------------------------------------------------------------------------
def build_old_to_new_index_map(
	new_order: list[int],
) -> dict[int, int]:
	"""
	Build the mapping old_atom_id -> new_atom_id from a 1-based ordering list.
	"""
	return {
		old_atom_id: new_atom_id
		for new_atom_id, old_atom_id in enumerate(new_order, start=1)
	}


def reorder_distance_dataframe(
	df_D,
	new_order: list[int],
) :
	"""
	Reorder the atom indices in a distance-constraint DataFrame.

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

	The function:
	1. remaps columns 0 and 1 using old_to_new_index_map
	2. enforces column 0 <= column 1
	3. swaps the associated columns when needed
	4. sorts the final DataFrame by column 0 and then column 1
	"""
	
	old_to_new_index_map = build_old_to_new_index_map(new_order)
	
	df_D_new = df_D.copy()

	# ------------------------------------------------------------------
	# Apply the new atom index map to the first two columns
	# ------------------------------------------------------------------
	df_D_new[0] = df_D_new[0].map(old_to_new_index_map)
	df_D_new[1] = df_D_new[1].map(old_to_new_index_map)

	# ------------------------------------------------------------------
	# Identify rows where the order must be swapped
	# We enforce column 0 >= column 1
	# ------------------------------------------------------------------
	swap_mask = df_D_new[0] < df_D_new[1]

	if swap_mask.any():

		# Swap columns 0 and 1
		df_D_new.loc[swap_mask, [0, 1]] = df_D_new.loc[swap_mask, [1, 0]].to_numpy()

		# Swap columns 2 and 3
		df_D_new.loc[swap_mask, [2, 3]] = df_D_new.loc[swap_mask, [3, 2]].to_numpy()

		# Swap columns 6 and 7
		df_D_new.loc[swap_mask, [6, 7]] = df_D_new.loc[swap_mask, [7, 6]].to_numpy()

		# Swap columns 8 and 9
		df_D_new.loc[swap_mask, [8, 9]] = df_D_new.loc[swap_mask, [9, 8]].to_numpy()

	# ------------------------------------------------------------------
	# Sort by the first column and then by the second column
	# ------------------------------------------------------------------
	df_D_new = df_D_new.sort_values(by=[0, 1]).reset_index(drop=True)

	return df_D_new
# -----------------------------------------------------------------------------------------------------
def reorder_instance(
	params: dict[str, object],
	out_dir: str | Path,
	pdb_id: str,
	ddgp_order_vec: int,
) -> int:
	"""
	Read the structure and distance files, compute the DDGP vertex ordering,
	reorder the data accordingly, and save the reordered files.

	The function returns a skip_flag (0 or 1) indicating whether the
	structure should be skipped.
	"""

	out_dir_i = Path(out_dir) / pdb_id
	chosen_model = int(params["model_number"])
	chosen_chain = str(params["chain_id"])
			
	Xfile = out_dir_i / f"X_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"
	Afile = out_dir_i / f"A_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"
	Ifile = out_dir_i / f"I_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"	
	Tfile = out_dir_i / f"T_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"
				
	# ------------------------------------------------------------------
	# Read input files
	# ------------------------------------------------------------------
	df_X = read_space_separated_file(Xfile)
	df_I = read_space_separated_file(Ifile)

	skip_flag = 0

	# ------------------------------------------------------------------
	# Determine number of residues
	# ------------------------------------------------------------------
	n = df_X.shape[0]
	nres = int(df_X.iat[n - 1, 2])

	new_order: list[int] = []

	# ------------------------------------------------------------------
	# Build the new atom ordering residue by residue
	# ------------------------------------------------------------------
	for k in range(nres):

		resnum_i = k + 1
		df_X_i = df_X[df_X[2] == resnum_i]

		skip_flag = validate_backbone_plus_hydrogens_residue(df_X_i, resnum_i, pdb_id, chosen_model, chosen_chain)

		if skip_flag == 1:
			break

		# First residue uses a special ordering
		if k == 0:
			new_order.extend(first_residue_order(df_X_i))

		# Internal residues use the DDGP ordering
		else:
			new_order.extend(get_internal_residue_numeric_order(df_X_i, ddgp_order_vec[k]))

	# ------------------------------------------------------------------
	# Exit early if the structure is invalid
	# ------------------------------------------------------------------
	if skip_flag == 1:
		return skip_flag

	# ------------------------------------------------------------------
	# Reorder structure and distance data
	# ------------------------------------------------------------------
	df_Xreord = reorder_and_renumber_structure_dataframe(df_X, new_order)

	df_Ireord = reorder_distance_dataframe(df_I, new_order)

	# ------------------------------------------------------------------
	# Save reordered files
	# ------------------------------------------------------------------
	save_distances_from_df_structure(df_Ireord, Ifile)

	save_coordinates_from_df_structure(df_Xreord, Xfile)

	return skip_flag
# -----------------------------------------------------------------------------------------------------

import shutil
from pathlib import Path

import numpy as np

from pdb_parser.utils import *
from pdb_parser.io import *
from pdb_parser.reordering import *

VALID_DISTANCE_MODES = {"precise", "interval_centered", "interval_experimental"}
VALID_ATOM_SELECTIONS = {
	"full_chain",
	"backbone",
	"backbone_plus_hydrogens",
	"backbone_plus_neighbors",
}
VALID_VDW_OPTIONS = {"yes", "no"}


def _resolve_effective_angular_width(
	distance_constraints: str,
	torsion_angle_width: float,
) -> float:
	"""Return the effective torsion width used by TALOS-like sampling."""
	if distance_constraints == "precise":
		return 0.0
	if torsion_angle_width < 0.0:
		return 360.0
	return torsion_angle_width


def _resolve_effective_backbone_torsion_percentage(
	distance_constraints: str,
	percentage_backbone_torsion_angles: float,
) -> float:
	"""Return the effective phi/psi selection percentage after parser policy."""
	if distance_constraints == "precise":
		return 100.0
	return min(max(percentage_backbone_torsion_angles, 0.0), 100.0)


def _required_seed_fields_for_execution(
	distance_constraints: str,
	angular_width: float,
	backbone_torsion_percentage: float,
	n_residues: int,
) -> set[str]:
	"""Return exactly the seed fields that can affect this execution."""
	required_fields: set[str] = set()

	if distance_constraints == "interval_centered":
		required_fields.add("seed_interval_centered_distances")

	if angular_width > 0.0:
		required_fields.add("seed_talos_angle_centers")

	n_detected = int(round(n_residues * (backbone_torsion_percentage / 100.0)))
	n_full_range = n_residues - n_detected
	if 0 < n_full_range < n_residues:
		required_fields.add("seed_phi_psi_mask")

	return required_fields

def _require_param(params: dict[str, object], key: str) -> str:
	"""Return a stripped parameter value or raise a clear configuration error."""
	raw_value = params.get(key)
	if raw_value is None:
		raise ValueError(f"Missing required parser parameter: {key}")

	value = str(raw_value).strip()
	if not value:
		raise ValueError(f"Parser parameter '{key}' cannot be empty.")

	return value


def _parse_int_param(params: dict[str, object], key: str) -> int:
	"""Parse an integer parser parameter with a user-facing error message."""
	value = _require_param(params, key)
	try:
		return int(value)
	except ValueError as exc:
		raise ValueError(f"Parser parameter '{key}' must be an integer, got {value!r}.") from exc


def _parse_float_param(params: dict[str, object], key: str) -> float:
	"""Parse a float parser parameter with a user-facing error message."""
	value = _require_param(params, key)
	try:
		return float(value)
	except ValueError as exc:
		raise ValueError(f"Parser parameter '{key}' must be a number, got {value!r}.") from exc


def _parse_choice_param(
	params: dict[str, object],
	key: str,
	valid_values: set[str],
) -> str:
	"""Parse a normalized choice parameter."""
	value = _require_param(params, key).lower()
	if value not in valid_values:
		raise ValueError(
			f"Unsupported value for parser parameter '{key}': {value!r}. "
			f"Expected one of {sorted(valid_values)}."
		)
	return value


def _preflight_parser_params(params: dict[str, object]) -> dict[str, object]:
	"""Validate parser parameters before creating output artifacts."""
	parsed: dict[str, object] = {
		"model_number": _parse_int_param(params, "model_number"),
		"chain_id": _require_param(params, "chain_id"),
		"atom_selection": _parse_choice_param(params, "atom_selection", VALID_ATOM_SELECTIONS),
		"distance_constraints": _parse_choice_param(params, "distance_constraints", VALID_DISTANCE_MODES),
		"vdw_constraints": _parse_choice_param(params, "vdw_constraints", VALID_VDW_OPTIONS),
		"torsion_angle_width": _parse_float_param(params, "torsion_angle_width"),
		"percentage_backbone_torsion_angles": _parse_float_param(params, "percentage_backbone_torsion_angles"),
	}

	distance_constraints = parsed["distance_constraints"]

	if distance_constraints in {"precise", "interval_centered"}:
		parsed["max_distance"] = _parse_float_param(params, "max_distance")

	if distance_constraints == "interval_centered":
		parsed["epsilon_short"] = _parse_float_param(params, "epsilon_short")
		parsed["epsilon_long"] = _parse_float_param(params, "epsilon_long")
		parsed["vdwr_hh"] = _parse_float_param(params, "vdwr_hh")

	elif distance_constraints == "interval_experimental":
		parsed["noe_strong"] = _parse_float_param(params, "noe_strong")
		parsed["noe_medium"] = _parse_float_param(params, "noe_medium")
		parsed["noe_weak"] = _parse_float_param(params, "noe_weak")
		parsed["vdwr_hh"] = _parse_float_param(params, "vdwr_hh")

	return parsed

def parser(
	params: dict[str, object],
	pdb_data_dir: str | Path,
	seed_dir: str | Path,
	out_dir: str | Path,
	pdb_id: str,
	remove_tmp_dir: bool = False,
) -> None:
	"""
	Run the full parser pipeline for one PDB/model/chain combination.
	"""

	# ------------------------------------------------------------------
	# Validate parameters and prepare protein chain
	# ------------------------------------------------------------------
	pdb_file_name = Path(pdb_data_dir) / f"{pdb_id}.pdb"
	parsed_params = _preflight_parser_params(params)
	chosen_model = int(parsed_params["model_number"])
	chosen_chain = str(parsed_params["chain_id"])
	atom_selection = str(parsed_params["atom_selection"])
	distance_constraints = str(parsed_params["distance_constraints"])
	protein_chain = ensure_nmr_model_chain_ready(pdb_file_name, chosen_model, chosen_chain)
	angular_width = _resolve_effective_angular_width(
		distance_constraints,
		float(parsed_params["torsion_angle_width"]),
	)
	backbone_torsion_percentage = _resolve_effective_backbone_torsion_percentage(
		distance_constraints,
		float(parsed_params["percentage_backbone_torsion_angles"]),
	)
	required_seed_fields = _required_seed_fields_for_execution(
		distance_constraints,
		angular_width,
		backbone_torsion_percentage,
		len(protein_chain.residues),
	)

	out_dir_i = Path(out_dir) / pdb_id
	ensure_dir(out_dir_i)

	tmp_dir = out_dir_i / "tmp"
	ensure_dir(tmp_dir)

	# ------------------------------------------------------------------
	# Save filtered atoms (TSV structure file)
	# ------------------------------------------------------------------
	tsv_structure_file = tmp_dir / "filtered_atoms.dat"

	save_filtered_atoms(tsv_structure_file, protein_chain, atom_selection)

	# ------------------------------------------------------------------
	# Load or create per-PDB random seeds only after preflight validation
	# ------------------------------------------------------------------
	seed_bundle = load_or_create_seed_bundle(
		seed_dir,
		pdb_id,
		chosen_model,
		chosen_chain,
		required_seed_fields,
	)
	rng_interval_centered = (
		np.random.default_rng(seed_bundle["seed_interval_centered_distances"])
		if seed_bundle["seed_interval_centered_distances"] is not None
		else None
	)
	rng_talos_angles = (
		np.random.default_rng(seed_bundle["seed_talos_angle_centers"])
		if seed_bundle["seed_talos_angle_centers"] is not None
		else None
	)
	rng_phi_psi_mask = (
		np.random.default_rng(seed_bundle["seed_phi_psi_mask"])
		if seed_bundle["seed_phi_psi_mask"] is not None
		else None
	)

	# ------------------------------------------------------------------
	# Save structure in fixed-width format (X file)
	# ------------------------------------------------------------------
	xfile = out_dir_i / f"X_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"

	convert_tsv_structure_to_pdb_format(tsv_structure_file, xfile)

	# ------------------------------------------------------------------
	# Build covalent topology
	# ------------------------------------------------------------------
	topology = build_covalent_and_planar_topology(tsv_structure_file)

	# ------------------------------------------------------------------
	# Distance constraints 1: covalent + planar
	# ------------------------------------------------------------------
	dcfile_1 = tmp_dir / "distance_constraints_1.dat"

	covalent_and_planar_distances(tsv_structure_file, topology, dcfile_1)

	# ------------------------------------------------------------------
	# Distance constraints 2: van der Waals
	# ------------------------------------------------------------------
	dcfile_2 = tmp_dir / "distance_constraints_2.dat"

	if parsed_params["vdw_constraints"] == "yes":
		vdw_distances(tsv_structure_file, topology, dcfile_2)

	# ------------------------------------------------------------------
	# Distance constraints 3: planar peptide geometry
	# ------------------------------------------------------------------
	dcfile_3 = tmp_dir / "distance_constraints_3.dat"

	planar_peptide_distances(tsv_structure_file, dcfile_3, atom_selection)

	# ------------------------------------------------------------------
	# Distance constraints 4: NMR (NOE) constraints
	# ------------------------------------------------------------------
	dcfile_4 = tmp_dir / "distance_constraints_4.dat"
	
	if distance_constraints == "interval_centered":
		nmr_distance_constraints(
			tsv_structure_file,
			dcfile_4,
			distance_constraints,
			atom_selection,
			float(parsed_params["epsilon_short"]),
			float(parsed_params["epsilon_long"]),
			float(parsed_params["max_distance"]),
			None,
			None,
			None,
			float(parsed_params["vdwr_hh"]),
			rng_interval_centered,
		)
	elif distance_constraints == "precise":
		nmr_distance_constraints(
			tsv_structure_file,
			dcfile_4,
			distance_constraints,
			atom_selection,
			None,
			None,
			float(parsed_params["max_distance"]),
			None,
			None,
			None,
			0.0,
			None,
		)
	elif distance_constraints == "interval_experimental":
		nmr_distance_constraints(
			tsv_structure_file,
			dcfile_4,
			distance_constraints,
			atom_selection,
			None,
			None,
			None,
			float(parsed_params["noe_strong"]),
			float(parsed_params["noe_medium"]),
			float(parsed_params["noe_weak"]),
			float(parsed_params["vdwr_hh"]),
			None,
		)
	else:
		raise ValueError(
			f"Unsupported distance_constraints mode: {distance_constraints!r}. "
			f"Expected one of {sorted(VALID_DISTANCE_MODES)}."
		)

	# ------------------------------------------------------------------
	# Angular constraints (TALOS-like)
	# ------------------------------------------------------------------
	acfile_1 = tmp_dir / f"angular_constraints_5.dat"

	talos_n_like(
		pdb_file_name,
		chosen_model,
		chosen_chain,
		acfile_1,
		0.0,
		angular_width,
		angular_width,
		backbone_torsion_percentage,
		rng_talos_angles,
		rng_phi_psi_mask,
	)
	
	acfile = out_dir_i / f"A_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"
	convert_tsv_angular_to_pdb_format(acfile_1, acfile)

	# ------------------------------------------------------------------
	# Convert angular intervals to distance intervals
	# ------------------------------------------------------------------
	dcfile_5 = tmp_dir / "distance_constraints_5.dat"

	backbone_angular_interval_to_distance_interval(tsv_structure_file, acfile_1, dcfile_5)

	# ------------------------------------------------------------------
	# Merge all distance constraint files
	# ------------------------------------------------------------------
	dcfile = out_dir_i / f"I_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"

	merge_distance_constraint_files(tmp_dir, dcfile)
	
	# ------------------------------------------------------------------
	# Distance constraints 5: tightening upper bounds
	# ------------------------------------------------------------------
	#if params.get("vdw_constraints", "").strip().lower() == "yes":
	#	tightening_upper_bounds(dcfile)

	# ------------------------------------------------------------------
	# Optional cleanup of temporary directory
	# ------------------------------------------------------------------
	if remove_tmp_dir and tmp_dir.exists():
		shutil.rmtree(tmp_dir)
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

	out_dir_pdb_id = Path(out_dir) / pdb_id
	chosen_model = int(params["model_number"])
	chosen_chain = str(params["chain_id"])
			
	Xfile = out_dir_pdb_id / f"X_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"
	Afile = out_dir_pdb_id / f"A_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"
	Ifile = out_dir_pdb_id / f"I_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"			
	# ------------------------------------------------------------------
	# Read input files
	# ------------------------------------------------------------------
	df_X = read_space_separated_file(Xfile)
	df_I = read_space_separated_file(Ifile)
	df_A = read_space_separated_file(Afile)

	skip_flag = 0

	# ------------------------------------------------------------------
	# Determine number of residues
	# ------------------------------------------------------------------
	n = df_X.shape[0]
	
	res_ids = df_X.iloc[:, 2].unique().tolist()
	# in case the protein chain does not starts from residue 1
	nres = res_ids[len(res_ids) - 1] - res_ids[0] + 1
	new_order: list[tuple[int, str]] = []

	# ------------------------------------------------------------------
	# Build the new atom ordering residue by residue
	# ------------------------------------------------------------------
	for k in range(nres):

		resnum_i = res_ids[k]
		df_X_i = df_X[df_X[2] == resnum_i]

		skip_flag = validate_backbone_plus_hydrogens_residue(df_X_i, resnum_i, pdb_id, chosen_model, chosen_chain, res_ids[0])

		if skip_flag == 1:
			break

		# First residue uses a special ordering
		if k == 0:
			new_order.extend(first_residue_order(df_X_i, res_ids[k]))

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
	atom_ids_new_order = [atom_id for atom_id, _ in new_order]
		
	df_Xreordered = reorder_structure_dataframe(df_X, atom_ids_new_order)

	df_Ireordered = reorder_distance_dataframe(df_I, atom_ids_new_order)
	
	atom_names_new_order = [atom_name for _, atom_name in new_order]

	# ------------------------------------------------------------------
	# Build the atom cliques residue by residue
	# ------------------------------------------------------------------
	atom_cliques: list[tuple[list[tuple[int, str]], int]] = []

	available_atoms = build_available_atoms(df_Xreordered)

	for k in range(nres):
		if k == 0:
			atom_cliques.extend(build_first_residue_pattern(available_atoms, res_ids[k]))
		else:
			atom_cliques.extend(build_ddgp_pattern_entries(ddgp_order_vec[k], res_ids[k], available_atoms))
	
	T = build_atom_clique_index_matrix(atom_cliques, df_Xreordered, df_A, df_Ireordered)

	# ------------------------------------------------------------------
	# Save reordered files
	# ------------------------------------------------------------------
	out_dir_pdb_id_reordered = Path(out_dir_pdb_id) / "reordered"
	ensure_dir(out_dir_pdb_id_reordered)
	chosen_order_id = int(params["order_id"])
	reordered_suffix = f"_ddgpHCorder{chosen_order_id}.dat"

	Xfile = out_dir_pdb_id_reordered / f"X_{pdb_id}_model{chosen_model}_chain{chosen_chain}{reordered_suffix}"
	Ifile = out_dir_pdb_id_reordered / f"I_{pdb_id}_model{chosen_model}_chain{chosen_chain}{reordered_suffix}"
	Tfile = out_dir_pdb_id_reordered / f"T_{pdb_id}_model{chosen_model}_chain{chosen_chain}{reordered_suffix}"

	save_distances_from_df_structure(df_Ireordered, Ifile)
	save_coordinates_from_df_structure(df_Xreordered, Xfile)
	save_cliques_from_matrix_T(T, Tfile)
	
	return skip_flag

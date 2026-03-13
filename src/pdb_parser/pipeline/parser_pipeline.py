import shutil
from pathlib import Path

from pdb_parser.utils import *
from pdb_parser.io import *
from pdb_parser.reordering import *

def parser(
	params: dict[str, object],
	pdb_data_dir: str | Path,
	out_dir: str | Path,
	pdb_id: str,
	remove_tmp_dir: bool = False,
) -> None:
	"""
	Run the full parser pipeline for one PDB/model/chain combination.
	"""
	
	out_dir_i = Path(out_dir) / pdb_id
	ensure_dir(out_dir_i)

	tmp_dir = out_dir_i / "tmp"
	ensure_dir(tmp_dir)

	# ------------------------------------------------------------------
	# Prepare protein chain
	# ------------------------------------------------------------------
	pdb_file_name = Path(pdb_data_dir) / f"{pdb_id}.pdb"
	chosen_model = int(params["model_number"])
	chosen_chain = str(params["chain_id"])

	protein_chain = ensure_nmr_model_chain_ready(pdb_file_name, chosen_model, chosen_chain)

	# ------------------------------------------------------------------
	# Save filtered atoms (TSV structure file)
	# ------------------------------------------------------------------
	tsv_structure_file = tmp_dir / "filtered_atoms.dat"

	save_filtered_atoms(tsv_structure_file, protein_chain, params.get("atom_selection", "").strip().lower())

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

	if params.get("vdw_constraints", "").strip().lower() == "yes":
		vdw_distances(tsv_structure_file, topology, dcfile_2)

	# ------------------------------------------------------------------
	# Distance constraints 3: planar peptide geometry
	# ------------------------------------------------------------------
	dcfile_3 = tmp_dir / "distance_constraints_3.dat"

	planar_peptide_distances(tsv_structure_file, dcfile_3)

	# ------------------------------------------------------------------
	# Distance constraints 4: NMR (NOE) constraints
	# ------------------------------------------------------------------
	dcfile_4 = tmp_dir / "distance_constraints_4.dat"

	nmr_distance_constraints(
		tsv_structure_file,
		dcfile_4,
		params.get("distance_constraints", "").strip().lower(),
		params.get("atom_selection", "").strip().lower(),
		float(params["epsilon_short"]),
		float(params["epsilon_long"]),
		float(params["max_distance"]),
		float(params["noe_strong"]),
		float(params["noe_medium"]),
		float(params["noe_weak"]),
		float(params["vdw_threshold"]),
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
		params.get("distance_constraints", "").strip().lower(),
		0.0,
		float(params["torsion_angle_width"]),
		float(params["torsion_angle_width"]),
		float(params["percentage_backbone_torsion_angles"]),
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
	# Optional cleanup of temporary directory
	# ------------------------------------------------------------------
	if remove_tmp_dir and tmp_dir.exists():
		shutil.rmtree(tmp_dir)
# -----------------------------------------------------------------------------------------------------
def sort_instance(
	params: dict[str, object],
	out_dir: str | Path,
	pdb_id: str,
	ddgp_order_vec: int,
) -> int:
	"""
	Read the structure and distance files, compute the DDGP vertex ordering,
	sort the data accordingly, and save the sorted files.

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
	nres = int(df_X.iat[n - 1, 2])

	new_order: list[tuple[int, str]] = []

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
	# Sort structure and distance data
	# ------------------------------------------------------------------
	atom_ids_new_order = [atom_id for atom_id, _ in new_order]
		
	df_Xreord = sort_structure_dataframe(df_X, atom_ids_new_order)

	df_Ireord = sort_distance_dataframe(df_I, atom_ids_new_order)
	
	atom_names_new_order = [atom_name for _, atom_name in new_order]

	# ------------------------------------------------------------------
	# Build the atom cliques residue by residue
	# ------------------------------------------------------------------
	atom_cliques: list[tuple[list[tuple[int, str]], int]] = []

	available_atoms = build_available_atoms(df_Xreord)

	for k in range(nres):
		if k == 0:
			atom_cliques.extend(build_first_residue_pattern(available_atoms))
		else:
			atom_cliques.extend(build_ddgp_pattern_entries(ddgp_order_vec[k], k + 1, available_atoms))

	T = build_atom_clique_index_matrix(atom_cliques, df_Xreord, df_A, df_Ireord)

	# ------------------------------------------------------------------
	# Save sorted files
	# ------------------------------------------------------------------
	out_dir_pdb_id_sorted = Path(out_dir_pdb_id) / "sorted"
	ensure_dir(out_dir_pdb_id_sorted)
	
	Xfile = out_dir_pdb_id_sorted / f"X_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"
	Ifile = out_dir_pdb_id_sorted / f"I_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"
	Tfile = out_dir_pdb_id_sorted / f"T_{pdb_id}_model{chosen_model}_chain{chosen_chain}.dat"
	
	save_distances_from_df_structure(df_Ireord, Ifile)
	save_coordinates_from_df_structure(df_Xreord, Xfile)
	save_cliques_from_matrix_T(T, Tfile)

	return skip_flag

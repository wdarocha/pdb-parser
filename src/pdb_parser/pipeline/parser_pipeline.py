import shutil
from pathlib import Path

from pdb_parser.utils import *
from pdb_parser.io import *

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

"""Generate TALOS-N-like backbone torsion intervals."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pdb_parser.io.angular_constraints import full_range_angle, wrap_angle_deg
from pdb_parser.io.pdb_ops import extract_model_chain
from pdb_parser.io.sequence import are_adjacent_residues, validate_residue_sequence
from pdb_parser.messages import print_ok, print_warning


def talos_n_like(
	pdb_file: str | Path,
	chosen_model: int,
	chosen_chain: str,
	output_angular_file: str | Path,
	omega_angle_width: float,
	phi_angle_width: float,
	psi_angle_width: float,
	percentage_backbone_torsion_angles: float,
	angle_rng: np.random.Generator | None = None,
	mask_rng: np.random.Generator | None = None,
	drop_nonpositive: bool = True,
) -> None:
	"""Generate TALOS-N-like omega, phi, and psi intervals from a PDB chain."""
	# Imported here so that importing the package (e.g. by the remap stage) does not load MDAnalysis.
	from MDAnalysis.lib.distances import calc_dihedrals

	def sanitize_angle_width(angle_width: float, angle_name: str) -> float:
		try:
			angle_width = float(angle_width)
		except (TypeError, ValueError):
			angle_width = 360.0
		if angle_width < 0.0:
			print_warning(f"Negative {angle_name} detected; using full range (360°).")
			return 360.0
		return angle_width

	def sample_angle_interval(base_angle: float | None, angle_width: float) -> tuple[float, float]:
		if base_angle is None:
			return full_range_angle()
		if angle_width <= 0.0:
			return base_angle, 0.0

		angle_radius = angle_width / 2.0
		base_angle = wrap_angle_deg(base_angle)
		while True:
			angle_center = wrap_angle_deg(
				float(angle_rng.normal(loc=base_angle, scale=angle_width / 8.0))
			)
			if abs(wrap_angle_deg(angle_center - base_angle)) < angle_radius:
				return angle_center, angle_radius

	angle_rng = angle_rng if angle_rng is not None else np.random.default_rng()
	mask_rng = mask_rng if mask_rng is not None else np.random.default_rng()
	omega_angle_width = sanitize_angle_width(omega_angle_width, "omega_angle_width")
	phi_angle_width = sanitize_angle_width(phi_angle_width, "phi_angle_width")
	psi_angle_width = sanitize_angle_width(psi_angle_width, "psi_angle_width")
	try:
		percentage_backbone_torsion_angles = float(percentage_backbone_torsion_angles)
	except (TypeError, ValueError):
		percentage_backbone_torsion_angles = 100.0
	if percentage_backbone_torsion_angles < 0.0:
		print_warning("percentage_backbone_torsion_angles < 0; clamped to 0%.")
		percentage_backbone_torsion_angles = 0.0
	elif percentage_backbone_torsion_angles > 100.0:
		print_warning("percentage_backbone_torsion_angles > 100; clamped to 100%.")
		percentage_backbone_torsion_angles = 100.0

	output_path = Path(output_angular_file)
	chain_residues = extract_model_chain(
		pdb_file,
		int(chosen_model),
		chosen_chain,
		drop_nonpositive=drop_nonpositive,
	).residues
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
	validate_residue_sequence(residue_ids, context=f"chain {chosen_chain}")
	omega_angles: list[float | None] = [None] * n_residues
	phi_angles: list[float | None] = [None] * n_residues
	psi_angles: list[float | None] = [None] * n_residues
	for index in range(n_residues):
		if index > 0:
			ca_prev, c_prev = ca_positions[index - 1], c_positions[index - 1]
			n_i, ca_i = n_positions[index], ca_positions[index]
			adjacent_to_previous = are_adjacent_residues(
				residue_ids, residue_ids[index - 1], residue_ids[index], context=f"chain {chosen_chain}"
			)
			if adjacent_to_previous and ca_prev is not None and c_prev is not None and n_i is not None and ca_i is not None:
				omega_angles[index] = wrap_angle_deg(float(np.degrees(calc_dihedrals(ca_prev.reshape(1, 3), c_prev.reshape(1, 3), n_i.reshape(1, 3), ca_i.reshape(1, 3))[0])))
			c_i = c_positions[index]
			if adjacent_to_previous and c_prev is not None and n_i is not None and ca_i is not None and c_i is not None:
				phi_angles[index] = wrap_angle_deg(float(np.degrees(calc_dihedrals(c_prev.reshape(1, 3), n_i.reshape(1, 3), ca_i.reshape(1, 3), c_i.reshape(1, 3))[0])))
		if index < n_residues - 1:
			n_i, ca_i, c_i, n_next = n_positions[index], ca_positions[index], c_positions[index], n_positions[index + 1]
			adjacent_to_next = are_adjacent_residues(
				residue_ids, residue_ids[index], residue_ids[index + 1], context=f"chain {chosen_chain}"
			)
			if adjacent_to_next and n_i is not None and ca_i is not None and c_i is not None and n_next is not None:
				psi_angles[index] = wrap_angle_deg(float(np.degrees(calc_dihedrals(n_i.reshape(1, 3), ca_i.reshape(1, 3), c_i.reshape(1, 3), n_next.reshape(1, 3))[0])))

	n_detected = int(round(n_residues * (percentage_backbone_torsion_angles / 100.0)))
	n_full_range = n_residues - n_detected
	if n_full_range > 0:
		full_range_indices = set(mask_rng.choice(np.arange(n_residues, dtype=int), size=n_full_range, replace=False))
	else:
		full_range_indices: set[int] = set()

	with output_path.open("w", encoding="utf-8") as output_file:
		output_file.write("resid\tresname\tomega_center\tomega_radius\tphi_center\tphi_radius\tpsi_center\tpsi_radius\n")
		for index, residue_id in enumerate(residue_ids):
			omega_center, omega_radius = sample_angle_interval(omega_angles[index], omega_angle_width)
			if index in full_range_indices:
				phi_center, phi_radius = full_range_angle()
				psi_center, psi_radius = full_range_angle()
			else:
				phi_center, phi_radius = sample_angle_interval(phi_angles[index], phi_angle_width)
				psi_center, psi_radius = sample_angle_interval(psi_angles[index], psi_angle_width)
			output_file.write(f"{residue_id}\t{residue_names[index]}\t{omega_center:.12f}\t{omega_radius:.12f}\t{phi_center:.12f}\t{phi_radius:.12f}\t{psi_center:.12f}\t{psi_radius:.12f}\n")

	print_ok(f"TALOS-N-like phi/psi intervals saved to: {output_angular_file}")

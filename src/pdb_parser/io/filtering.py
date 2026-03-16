from __future__ import annotations
from pathlib import Path
from typing import Optional, Set

from pdb_parser.utils import ensure_dir

_BACKBONE: Set[str] = {"N", "CA", "C"}
# -----------------------------------------------------------------------------------------------------
def _require_backbone(present: Set[str], *, resname: str, resid: int, mode: str) -> None:
	"""Ensure N, CA, C exist when the mode depends on backbone atoms."""
	missing = [a for a in _BACKBONE if a not in present]
	if missing:
		raise ValueError(f"Residue {resname}_{resid}: missing backbone atoms {missing} for option '{mode}'.")
# -----------------------------------------------------------------------------------------------------
def _require_any_H_variant(present: Set[str], *, resname: str, resid: int, mode: str) -> None:
	"""Require at least one of H/H1/H2/H3 to exist."""
	if not any(h in present for h in ("H", "H1", "H2", "H3")):
		raise ValueError(f"Residue {resname}_{resid}: missing N-H variant; need one of ['H','H1','H2','H3'] for option '{mode}'.")
# -----------------------------------------------------------------------------------------------------
def _allowed_names_for_residue(
	resname: str,
	mode: str,
	present_names: Set[str],
	*,
	resid: int,
) -> Optional[Set[str]]:
	"""Return allowed atom names for a residue under the given mode.
	
	Returns
	-------
	None -> 'full_chain' (keep all)
	Set[str] -> allowed names for other modes; raises on missing required atoms.
	"""
	m = mode.strip().lower()
	r = resname.upper()
	
	if m == "full_chain":
		return None
	
	if m == "backbone":
		_require_backbone(present_names, resname=r, resid=resid, mode=m)
		return set(_BACKBONE)
	
	# For both plus_* modes start by requiring backbone
	_require_backbone(present_names, resname=r, resid=resid, mode=m)
	allowed = set(_BACKBONE)
	
	# Helper: add N-H variants if present (never required for PRO)
	def add_n_h_variants_if_present():
		for h in ("H", "H1", "H2", "H3"):
			if h in present_names:
				allowed.add(h)
	
	if m == "backbone_plus_hydrogens":
		# For non-PRO residues, require at least one N-H variant
		if r != "PRO":
			_require_any_H_variant(present_names, resname=r, resid=resid, mode=m)
		# Add any N-H variants that are present (optional for PRO)
		add_n_h_variants_if_present()
		
		if r == "GLY":
			# Need HA2 or, if missing, HA3
			if "HA2" in present_names:
				allowed.add("HA2")
			elif "HA3" in present_names:
				allowed.add("HA3")
			else:
				raise ValueError(f"Residue GLY_{resid}: missing HA2/HA3 required for option '{m}'.")
		elif r == "PRO":
			# Need HA and one of HD2/HD3 (N-H not required)
			if "HA" not in present_names:
				raise ValueError(f"Residue PRO_{resid}: missing HA required for option '{m}'.")
			allowed.add("HA")
			if "HD2" in present_names:
				allowed.add("HD2")
			elif "HD3" in present_names:
				allowed.add("HD3")
			else:
				raise ValueError(f"Residue PRO_{resid}: missing one of ['HD2','HD3'] required for option '{m}'.")
		else:
			# Other residues: need HA
			if "HA" not in present_names:
				raise ValueError(f"Residue {r}_{resid}: missing HA required for option '{m}'.")
			allowed.add("HA")
		
		return allowed
	
	if m == "backbone_plus_neighbors":
		# For non-PRO residues, require at least one N–H variant
		if r != "PRO":
			_require_any_H_variant(present_names, resname=r, resid=resid, mode=m)
		# Add any N–H variants that are present (optional for PRO)
		add_n_h_variants_if_present()
		
		if r == "GLY":
			required = {"HA2", "HA3", "O"}
			missing = [a for a in required if a not in present_names]
			if missing:
				raise ValueError(f"Residue GLY_{resid}: missing {missing} required for option '{m}'.")
			allowed.update(required)
		elif r == "PRO":
			required = {"HA", "CD", "CB", "O"}
			missing = [a for a in required if a not in present_names]
			if missing:
				raise ValueError(f"Residue PRO_{resid}: missing {missing} required for option '{m}'.")
			allowed.update(required)
		else:
			required = {"HA", "CB", "O"}
			missing = [a for a in required if a not in present_names]
			if missing:
				raise ValueError(f"Residue {r}_{resid}: missing {missing} required for option '{m}'.")
			allowed.update(required)
		
		return allowed
	
	raise ValueError("Invalid mode. Use one of: full_chain, backbone, backbone_plus_hydrogens, backbone_plus_neighbors")
# -----------------------------------------------------------------------------------------------------
def save_filtered_atoms(
	fname: str | Path,
	protein_chain, # MDAnalysis AtomGroup (from ensure_nmr_model_chain_ready)
	mode: str, # 'full_chain' | 'backbone' | 'backbone_plus_hydrogens' | 'backbone_plus_neighbors'
) -> None:
	"""Filter atoms from 'protein_chain' according to 'mode' and save a sorted TSV.
	
	- Validates required atoms per residue; raises ValueError on missing ones.
	- Sorts by original atom.id for determinism.
	- Renumbers atom_id in the output (1..N) to avoid gaps.
	"""
	
	kept = []
	for res in protein_chain.residues:
		resname = res.resname.upper()
		resid = int(res.resid)
		present = {a.name.upper() for a in res.atoms}
		
		allowed = _allowed_names_for_residue(resname, mode, present, resid=resid)
		
		if allowed is None:   # full_chain
			kept.extend(res.atoms)
		else:
			for a in res.atoms:
				if a.name.upper() in allowed:
					kept.append(a)
	
	# Deterministic order by original PDB atom id
	kept.sort(key=lambda a: int(a.id))
	
	# Write TSV with renumbered atom_id (1..N)
	with fname.open("w", encoding="utf-8") as f:
		f.write("atom_id\tatom_name\tresid\tresname\tx\ty\tz\n")
		for new_id, a in enumerate(kept, start=1):
			x, y, z = a.position
			f.write(f"{new_id}\t{a.name.upper()}\t{int(a.resid)}\t{a.resname.upper()}\t{x:.3f}\t{y:.3f}\t{z:.3f}\n")
	
	print(f"[OK] Selected atoms from PDB file saved to: {fname}")
# -----------------------------------------------------------------------------------------------------
def convert_tsv_structure_to_pdb_format(
	input_file: str | Path,
	output_file: str | Path,
) -> None:
	"""
	Read a TSV structure file and write it in fixed-width format.

	Input columns expected:
	atom_id, atom_name, resid, resname, x, y, z

	Output format per line:
	%5d %-4s %4d %3s %8.3f%8.3f%8.3f
	"""

	input_file = Path(input_file)
	output_file = Path(output_file)

	if not input_file.exists():
		raise ValueError(f"Input file not found: {input_file}")

	with input_file.open("r", encoding="utf-8") as fin, \
		 output_file.open("w", encoding="utf-8") as fout:

		# Skip header line
		next(fin)

		for line in fin:
			line = line.strip()
			if not line:
				continue

			parts = line.split()

			atom_id = int(parts[0])
			atom_name = parts[1]
			resid = int(parts[2])
			resname = parts[3]

			x = float(parts[4])
			y = float(parts[5])
			z = float(parts[6])

			fout.write("%5d %-4s %4d %3s %8.3f%8.3f%8.3f\n"	% (atom_id, atom_name, resid, resname, x, y, z))
	
	print(f"[OK] Structure file saved to: {output_file}")
# -----------------------------------------------------------------------------------------------------
def convert_tsv_angular_to_pdb_format(
	input_file: str | Path,
	output_file: str | Path,
) -> None:
	"""
	Convert a TSV angular-constraint file to a fixed-width text format.

	Input columns expected:
	resid, resname, omega_center, omega_radius,
	phi_center, phi_radius, psi_center, psi_radius

	Output format per line:
	%4d %3s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f
	"""
	input_file = Path(input_file)
	output_file = Path(output_file)

	if not input_file.exists():
		raise ValueError(f"Input file not found: {input_file}")

	output_file.parent.mkdir(parents=True, exist_ok=True)

	with input_file.open("r", encoding="utf-8") as input_stream, \
		 output_file.open("w", encoding="utf-8") as output_stream:

		next(input_stream, None)

		for line in input_stream:
			line = line.strip()
			if not line:
				continue

			parts = line.split()

			if len(parts) < 8:
				raise ValueError(f"Invalid angular TSV line: {line}")

			resid = int(parts[0])
			resname = parts[1]
			omega_center = float(parts[2])
			omega_radius = float(parts[3])
			phi_center = float(parts[4])
			phi_radius = float(parts[5])
			psi_center = float(parts[6])
			psi_radius = float(parts[7])

			output_stream.write("%4d %3s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n" % (
				resid, resname, 
				omega_center, omega_radius,
				phi_center, phi_radius,
				psi_center, psi_radius)
			)
	print(f"[OK] Angular constraint file saved to: {output_file}")
# -----------------------------------------------------------------------------------------------------

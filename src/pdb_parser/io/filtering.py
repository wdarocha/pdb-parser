from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from pdb_parser.messages import print_ok, print_warning
from pdb_parser.remapping.ordering import (
	choose_first_existing_atom,
	first_residue_order,
	get_ha_group,
	get_nonterminal_h_group,
)
from pdb_parser.utils import ensure_dir

BACKBONE_ATOMS = {"N", "CA", "C"}
VALID_FILTERING_MODES = (
	"full_chain",
	"backbone",
	"backbone_plus_hydrogens",
	"backbone_plus_neighbors",
)


def _deduplicate_atom_group(protein_chain):
	"""Keep one deterministic conformer for each residue/atom identity.

	The filtered TSV format does not preserve alternate-location identifiers, so
	writing more than one conformer would make the output ambiguous and invalid
	for all downstream topology consumers.
	"""
	chosen: dict[tuple[int, str], tuple[str, int]] = {}
	duplicate_keys: set[tuple[int, str]] = set()
	for index, atom in enumerate(protein_chain):
		key = (int(atom.resid), str(atom.name).strip().upper())
		altloc = ""
		if hasattr(protein_chain, "altLocs"):
			altloc = str(protein_chain.altLocs[index]).strip()
		candidate = (altloc, index)
		previous = chosen.get(key)
		if previous is not None:
			duplicate_keys.add(key)
		if previous is None or (
			altloc == "" and previous[0] != ""
		) or (
			altloc != "" and previous[0] != "" and altloc < previous[0]
		):
			chosen[key] = candidate

	if not duplicate_keys:
		return protein_chain

	keep_indices = sorted(index for _, index in chosen.values())
	print_warning(
		"altloc_policy='prefer_blank': discarded duplicate alternate atom "
		f"records for {len(duplicate_keys)} residue/atom identities before writing "
		"the filtered file."
	)
	return protein_chain[keep_indices]


def _drop_coincident_atoms(protein_chain):
	"""Keep the first atom (file order) of every set of atoms with identical coordinates.

	A deposited model can place two atoms at the same position (1PIJ: H1 and H2 of residue 1).
	The zero distance between them breaks every downstream geometry (torsions divide by it),
	so the later atoms are dropped with a warning.
	"""
	seen: dict[tuple[float, float, float], str] = {}
	keep: list[int] = []
	dropped: list[str] = []
	for index, atom in enumerate(protein_chain):
		key = tuple(np.round(np.asarray(atom.position, dtype=float), 3))
		label = f"{str(atom.resname).strip().upper()}_{int(atom.resid)} {str(atom.name).strip().upper()}"
		if key in seen:
			dropped.append(f"{label} (same coordinates as {seen[key]})")
			continue
		seen[key] = label
		keep.append(index)
	if not dropped:
		return protein_chain
	print_warning(
		f"discarded {len(dropped)} atom(s) whose coordinates repeat an earlier atom: "
		+ "; ".join(dropped)
	)
	return protein_chain[keep]


def _nearest_proline_hd(residue_atoms, previous_residue_atoms) -> str | None:
	"""Return HD2 or HD3, whichever is nearer the alpha hydrogen of the previous residue.

	The remap places h(k) of a proline from ha(k-1) (the HA(k-1)-HD(k) NOE), so the nearer
	delta hydrogen gives the stronger, always available restraint; in cis prolines HD2 can be
	farther than the 5 A NOE cutoff. Ties, a missing HD or a missing previous alpha hydrogen
	return ``None`` (the default first-present choice applies).
	"""
	positions = {str(atom.name).strip().upper(): atom.position for atom in residue_atoms}
	previous = {str(atom.name).strip().upper(): atom.position for atom in previous_residue_atoms}
	ha = next((previous[name] for name in ("HA", "HA2", "HA3") if name in previous), None)
	if ha is None or "HD2" not in positions or "HD3" not in positions:
		return None
	d2 = float(np.linalg.norm(np.asarray(positions["HD2"]) - np.asarray(ha)))
	d3 = float(np.linalg.norm(np.asarray(positions["HD3"]) - np.asarray(ha)))
	return "HD3" if d3 < d2 else "HD2"


def _require_backbone(
	present: set[str],
	*,
	resname: str,
	resid: int,
	mode: str,
) -> None:
	"""Ensure N, CA, and C exist when the mode depends on backbone atoms."""
	missing = [atom_name for atom_name in BACKBONE_ATOMS if atom_name not in present]
	if missing:
		raise ValueError(
			f"Residue {resname}_{resid}: missing backbone atoms {missing} for option '{mode}'."
		)


def _backbone_neighbor_names(
	resname: str,
	present_names: set[str],
	*,
	include_heavy_atoms: bool,
	include_hydrogens: bool,
	include_proline_hydrogens: bool = False,
) -> set[str]:
	"""Return direct backbone-neighbor atoms that are actually present.

	The PDB atom set is authoritative.  This deliberately makes glycine and
	partial experimental structures work without residue-specific requirements.
	"""
	allowed: set[str] = set()
	if include_heavy_atoms:
		allowed.update(name for name in ("O", "CB") if name in present_names)
		if resname == "PRO" and "CD" in present_names:
			allowed.add("CD")
	if include_hydrogens:
		# H/H1... are N neighbours; HA/HA2/HA3 are CA neighbours.  HD2/HD3
		# close proline's N-backbone bond.
		for name in present_names:
			if name in {"H", "HN"} or name.startswith("H1") or name.startswith("H2") or name.startswith("H3"):
				allowed.add(name)
			elif name in {"HA", "HA2", "HA3"}:
				allowed.add(name)
			elif include_proline_hydrogens and resname == "PRO" and name in {"HD2", "HD3"}:
				allowed.add(name)
	return allowed


def _remap_backbone_hydrogens(
	resname: str,
	present_names: set[str],
	*,
	resid: int,
	is_first: bool,
	proline_hd: str | None = None,
) -> set[str]:
	"""Backbone hydrogens kept by ``backbone_plus_hydrogens``: exactly those the remap orders.

	The remap places one atom per position and resolves each position with the first present
	candidate (``remapping/ordering.py``). Keeping a second alternative (HA3 next to HA2 in
	Gly, HD3 next to HD2 in Pro) left atoms in the X/I files without a vertex in the remap
	order, and an N-terminal Gly or Pro with both alternatives no longer matched the 5/6/7-atom
	first-residue patterns. The selection therefore delegates to the remap's own choice:

	* first residue: the hydrogens of ``first_residue_order`` (H3/H2/H1/H, HD2/HD3 in Pro);
	* other residues: the first present of ``get_nonterminal_h_group`` (H; HD2, HD3 in Pro)
	  and of ``get_ha_group`` (HA; HA2, HA3 in Gly). In Pro, ``proline_hd`` (the delta
	  hydrogen nearer the previous alpha hydrogen) takes the place of the first present.

	When the remap cannot order the residue (e.g. no N hydrogen), the previous permissive
	choice is kept so that the parser stage still runs; the remap then rejects the residue
	with its own explicit message.
	"""
	permissive = _backbone_neighbor_names(
		resname,
		present_names,
		include_heavy_atoms=False,
		include_hydrogens=True,
		include_proline_hydrogens=True,
	)
	if is_first:
		candidates = [name for name in sorted(present_names)
			if name in BACKBONE_ATOMS or name in permissive]
		frame = pd.DataFrame([[index, name, resid, resname] for index, name in enumerate(candidates, 1)])
		try:
			return {name for _, name in first_residue_order(frame, resid)} - set(BACKBONE_ATOMS)
		except ValueError:
			return permissive

	h_atom = (
		proline_hd
		if resname == "PRO" and proline_hd in present_names
		else choose_first_existing_atom(present_names, get_nonterminal_h_group(resname))
	)
	if h_atom is None and resname != "PRO" and "HN" in present_names:
		h_atom = "HN"  # kept for parser-only use; the remap requires the name H
	ha_atom = choose_first_existing_atom(present_names, get_ha_group(resname))
	return {name for name in (h_atom, ha_atom) if name is not None}


def _allowed_names_for_residue(
	resname: str,
	mode: str,
	present_names: set[str],
	*,
	resid: int,
	method: str,
	is_first: bool = False,
	proline_hd: str | None = None,
) -> Optional[set[str]]:
	"""Return allowed atom names for one residue under the selected public/internal mode."""
	normalized_mode = mode.strip().lower()
	normalized_resname = resname.upper()

	if normalized_mode == "full_chain":
		return None

	if normalized_mode == "backbone":
		_require_backbone(
			present_names,
			resname=normalized_resname,
			resid=resid,
			mode=normalized_mode,
		)
		return set(BACKBONE_ATOMS)

	_require_backbone(
		present_names,
		resname=normalized_resname,
		resid=resid,
		mode=normalized_mode,
	)
	allowed = set(BACKBONE_ATOMS)

	if normalized_mode == "backbone_plus_hydrogens":
		allowed.update(
			_remap_backbone_hydrogens(
				normalized_resname,
				present_names,
				resid=resid,
				is_first=is_first,
				proline_hd=proline_hd,
			)
		)
		return allowed

	if normalized_mode == "backbone_plus_neighbors":
		allowed.update(
			_backbone_neighbor_names(
				normalized_resname,
				present_names,
				include_heavy_atoms=True,
				include_hydrogens=method.strip().lower() == "nmr",
			)
		)
		return allowed

	raise ValueError(
		f"Invalid mode. Use one of: {', '.join(VALID_FILTERING_MODES)}"
	)


def resolve_effective_atom_selection(method: str, mode: str, *, warn: bool = True) -> str:
	"""Resolve method-specific atom-selection fallbacks in one place."""
	normalized_method = method.strip().lower()
	normalized_mode = mode.strip().lower()
	if normalized_method == "xray" and normalized_mode == "backbone_plus_hydrogens":
		if warn:
			print_warning(
				"atom_selection='backbone_plus_hydrogens' is unavailable for X-ray coordinates; "
				"using 'backbone_plus_neighbors'."
			)
		return "backbone_plus_neighbors"
	return normalized_mode


def save_filtered_atoms(
	fname: str | Path,
	protein_chain,
	mode: str,
	*,
	method: str = "nmr",
) -> None:
	"""Filter atoms from `protein_chain` according to `mode` and save a sorted TSV."""
	output_path = Path(fname)
	ensure_dir(output_path.parent)
	mode = resolve_effective_atom_selection(method, mode)
	protein_chain = _drop_coincident_atoms(_deduplicate_atom_group(protein_chain))
	kept_atoms = []

	atoms_by_residue: dict[tuple[int, str], list] = {}
	for atom in protein_chain:
		key = (int(atom.resid), str(atom.resname).strip().upper())
		atoms_by_residue.setdefault(key, []).append(atom)

	first_resid = min((resid for resid, _ in atoms_by_residue), default=None)
	atoms_by_resid = {resid: atoms for (resid, _), atoms in atoms_by_residue.items()}
	for (resid, resname), residue_atoms in atoms_by_residue.items():
		present_names = {str(atom.name).strip().upper() for atom in residue_atoms}
		proline_hd = (
			_nearest_proline_hd(residue_atoms, atoms_by_resid.get(resid - 1, []))
			if resname == "PRO" and resid != first_resid
			else None
		)
		allowed = _allowed_names_for_residue(
			resname, mode, present_names, resid=resid, method=method, is_first=resid == first_resid,
			proline_hd=proline_hd,
		)

		if allowed is None:
			kept_atoms.extend(residue_atoms)
		else:
			kept_atoms.extend(
				atom for atom in residue_atoms if str(atom.name).strip().upper() in allowed
			)

	kept_atoms.sort(key=lambda atom: int(atom.id))

	with output_path.open("w", encoding="utf-8") as output_stream:
		output_stream.write("atom_id\tatom_name\tresid\tresname\tx\ty\tz\n")
		for new_id, atom in enumerate(kept_atoms, start=1):
			x, y, z = atom.position
			output_stream.write(
				f"{new_id}\t{atom.name.upper()}\t{int(atom.resid)}\t"
				f"{atom.resname.upper()}\t{x:.3f}\t{y:.3f}\t{z:.3f}\n"
			)

	print_ok(f"Selected atoms from PDB file saved to: {output_path}")


# Minimum widths keep columns aligned; the explicit spaces keep them separable at any magnitude.
STRUCTURE_ROW_FORMAT = "%5d %-4s %4d %3s %9.3f %9.3f %9.3f\n"


def convert_tsv_structure_to_pdb_format(
	input_file: str | Path,
	output_file: str | Path,
) -> None:
	"""Read one TSV structure file and rewrite it as the aligned, whitespace-separated X layout.

	Every column is separated by at least one space. ``%9.3f`` fits every coordinate of the
	PDB format (-999.999 to 9999.999) with room to spare, and it is only a minimum width: a
	wider value widens its column instead of touching the previous one, so
	``read_space_separated_file`` always sees seven fields.
	"""
	input_path = Path(input_file)
	output_path = Path(output_file)

	if not input_path.exists():
		raise ValueError(f"Input file not found: {input_path}")

	with input_path.open("r", encoding="utf-8") as input_stream, output_path.open(
		"w",
		encoding="utf-8",
	) as output_stream:
		next(input_stream)

		for line in input_stream:
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

			output_stream.write(
				STRUCTURE_ROW_FORMAT
				% (atom_id, atom_name, resid, resname, x, y, z)
			)

	print_ok(f"Structure file saved to: {output_path}")


def convert_tsv_angular_to_pdb_format(
	input_file: str | Path,
	output_file: str | Path,
) -> None:
	"""Convert one TSV angular file to the public fixed-width parser layout."""
	input_path = Path(input_file)
	output_path = Path(output_file)

	if not input_path.exists():
		raise ValueError(f"Input file not found: {input_path}")

	output_path.parent.mkdir(parents=True, exist_ok=True)

	with input_path.open("r", encoding="utf-8") as input_stream, output_path.open(
		"w",
		encoding="utf-8",
	) as output_stream:
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

			output_stream.write(
				"%4d %3s %17.12f %17.12f %17.12f %17.12f %17.12f %17.12f\n"
				% (
					resid,
					resname,
					omega_center,
					omega_radius,
					phi_center,
					phi_radius,
					psi_center,
					psi_radius,
				)
			)

	print_ok(f"Angular constraint file saved to: {output_path}")

"""PDB I/O and model/chain utilities (MDAnalysis-based)."""

from __future__ import annotations

import os
import re
import tempfile
import time
from pathlib import Path
from typing import TYPE_CHECKING, Iterable

import numpy as np
import requests

from pdb_parser.io.sequence import validate_residue_sequence
from pdb_parser.messages import print_warning
from pdb_parser.utils.fs import ensure_dir

if TYPE_CHECKING:
	import MDAnalysis as mda
	from MDAnalysis.core.groups import AtomGroup


def _load_universe(pdb_path: str | Path) -> mda.Universe:
	"""Open a PDB file with MDAnalysis, imported on first use (the remap stage never needs it)."""
	import MDAnalysis as mda

	return mda.Universe(str(pdb_path), multiframe=True)

_PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")

# -----------------------------------------------------------------------------------------------------
def _normalize_pdb_id(pdb_id: str) -> str:
	"""Return normalized (uppercased) PDB id and validate format."""
	pid = pdb_id.strip().upper()
	if not _PDB_ID_RE.match(pid):
		raise ValueError(f"Invalid PDB id '{pdb_id}'. Expected 4 alphanumeric chars.")
	return pid
# -----------------------------------------------------------------------------------------------------
def _normalize_chain_id(chain_id: str) -> str:
	"""Return a stripped chain identifier and reject empty values."""
	cid = str(chain_id).strip()
	if not cid:
		raise ValueError("chain_id cannot be empty.")
	return cid
# -----------------------------------------------------------------------------------------------------
def _strip_string_array(values: Iterable[object]) -> np.ndarray:
	"""Return a NumPy array of stripped string values preserving length."""
	return np.array([str(value).strip() for value in values], dtype=object)
# -----------------------------------------------------------------------------------------------------
def _get_protein_atom_group(universe: mda.Universe) -> AtomGroup:
	"""Return protein atoms from the current frame, keeping only ATOM records when available."""
	protein_atoms = universe.select_atoms("protein")
	try:
		protein_atoms = protein_atoms.select_atoms("record_type ATOM")
	except (AttributeError, *_mdanalysis_errors()):
		# Some Universes do not expose record_type metadata; keep the broader protein selection.
		return protein_atoms
	return protein_atoms
# -----------------------------------------------------------------------------------------------------
def _mdanalysis_errors() -> tuple[type[Exception], ...]:
	"""Return the MDAnalysis selection errors, imported on first use."""
	from MDAnalysis.exceptions import NoDataError, SelectionError

	return NoDataError, SelectionError
# -----------------------------------------------------------------------------------------------------
def _list_chain_ids_from_atom_group(atom_group: AtomGroup) -> list[str]:
	"""Collect non-empty chain identifiers from both segid and chainID fields."""
	chain_ids: set[str] = set()

	if hasattr(atom_group, "segids"):
		chain_ids.update(text for text in _strip_string_array(atom_group.segids) if text)

	if hasattr(atom_group, "chainIDs"):
		chain_ids.update(text for text in _strip_string_array(atom_group.chainIDs) if text)

	return sorted(chain_ids)
# -----------------------------------------------------------------------------------------------------
def _select_chain_atoms(atom_group: AtomGroup, chain_id: str) -> AtomGroup:
	"""Select a chain from an AtomGroup using segid first and chainID second."""
	cid = _normalize_chain_id(chain_id)
	seg_matches = np.array([], dtype=int)
	chain_matches = np.array([], dtype=int)
	if hasattr(atom_group, "segids"):
		seg_matches = np.flatnonzero(_strip_string_array(atom_group.segids) == cid)
	if hasattr(atom_group, "chainIDs"):
		chain_matches = np.flatnonzero(_strip_string_array(atom_group.chainIDs) == cid)
	if seg_matches.size and chain_matches.size and not np.array_equal(seg_matches, chain_matches):
		raise ValueError(
			f"Chain identifier '{cid}' is ambiguous: it matches different segid and chainID selections."
		)
	match_indices = seg_matches if seg_matches.size else chain_matches
	if match_indices.size > 0:
		return atom_group[match_indices]

	return atom_group[:0]
# -----------------------------------------------------------------------------------------------------
def download_pdb(pdb_id: str, data_dir: str | Path) -> Path:
	"""Download a PDB file (text) from RCSB if not present; return local Path.
	
	Notes
	-----
	- Uses files.rcsb.org canonical URL.
	- Raises on HTTP errors.
	"""
	pid = _normalize_pdb_id(pdb_id)
	out_dir = ensure_dir(data_dir)
	filename = out_dir / f"{pid}.pdb"
	if filename.exists():
		try:
			_validate_downloaded_pdb(filename.read_text(encoding="utf-8", errors="replace"), pid)
			return filename
		except (OSError, ValueError):
			# Keep the old file until a replacement has been completely validated
			# and written.  A truncated cache must not poison future executions.
			pass

	text = _fetch_pdb_text(f"https://files.rcsb.org/download/{pid}.pdb", pid)
	_validate_downloaded_pdb(text, pid)
	fd, temporary_name = tempfile.mkstemp(prefix=f".{filename.name}.", dir=out_dir)
	try:
		with os.fdopen(fd, "w", encoding="utf-8") as stream:
			stream.write(text)
			stream.flush()
			os.fsync(stream.fileno())
		os.replace(temporary_name, filename)
	finally:
		Path(temporary_name).unlink(missing_ok=True)

	return filename


DOWNLOAD_ATTEMPTS = 3
DOWNLOAD_RETRY_DELAYS = (2.0, 5.0)  # seconds before the second and third attempts


def _fetch_pdb_text(url: str, pdb_id: str) -> str:
	"""Download one PDB file, retrying transient failures (timeouts, connection, 5xx, 429).

	A definitive HTTP answer (for example 404, an unknown id) is a data issue
	(``ValueError``); a network failure on every attempt is a ``ConnectionError``.
	"""
	error: Exception | None = None
	for attempt in range(DOWNLOAD_ATTEMPTS):
		try:
			response = requests.get(url, timeout=(5, 30))
			response.raise_for_status()
			return response.text
		except requests.HTTPError as exc:
			status = exc.response.status_code if exc.response is not None else None
			if status is not None and status != 429 and status < 500:
				raise ValueError(f"Failed to download PDB '{pdb_id}' from RCSB: {exc}") from exc
			error = exc
		except requests.RequestException as exc:
			error = exc
		if attempt + 1 < DOWNLOAD_ATTEMPTS:
			time.sleep(DOWNLOAD_RETRY_DELAYS[min(attempt, len(DOWNLOAD_RETRY_DELAYS) - 1)])
	raise ConnectionError(
		f"Failed to download PDB '{pdb_id}' from RCSB after {DOWNLOAD_ATTEMPTS} attempts: {error}"
	) from error


def _validate_downloaded_pdb(content: str, pdb_id: str) -> None:
	"""Reject truncated/non-PDB responses before they reach the parser."""
	if not content.strip() or "EXPDTA" not in content or not re.search(r"^(ATOM  |HETATM)", content, re.MULTILINE):
		raise ValueError(f"Downloaded PDB '{pdb_id}' is invalid or incomplete.")
# -----------------------------------------------------------------------------------------------------
def is_nmr_structure(pdb_path: str | Path, max_lines: int = 400) -> bool:
	"""Backward-compatible NMR predicate based on the shared EXPDTA classifier."""
	try:
		return detect_experimental_method(pdb_path, max_lines=max_lines) == "nmr"
	except ValueError:
		return False


def detect_experimental_method(pdb_path: str | Path, max_lines: int = 400) -> str:
	"""Classify the PDB EXPDTA record as ``nmr`` or ``xray``.

	Only these two methods are supported by the parser.  Keeping this decision
	here avoids silently treating a different experimental technique as X-ray.
	"""
	p = Path(pdb_path)
	detected_methods: set[str] = set()
	with p.open("r", encoding="utf-8", errors="replace") as fh:
		for index, line in enumerate(fh):
			if index > max_lines:
				break
			if not line.startswith("EXPDTA"):
				continue
			expdta = line.upper()
			if "NMR" in expdta:
				detected_methods.add("nmr")
			if "X-RAY" in expdta or "X RAY" in expdta:
				detected_methods.add("xray")
	if len(detected_methods) == 1:
		return detected_methods.pop()
	if len(detected_methods) > 1:
		raise ValueError(
			"Ambiguous EXPDTA method: multiple supported methods were declared "
			f"({', '.join(sorted(detected_methods))})."
		)
	raise ValueError("Unsupported or missing EXPDTA method; only NMR and X-RAY DIFFRACTION are supported.")
# -----------------------------------------------------------------------------------------------------
def list_number_of_models(pdb_path) -> int:
	"""
	Return the number of models available in the PDB file.
	"""
	u = _load_universe(pdb_path)
	return len(u.trajectory)
# -----------------------------------------------------------------------------------------------------
def list_chains_for_model(pdb_path: str | Path, model_number: int) -> list[str]:
	"""Return a sorted list of chain identifiers for the given 1-based model_number."""
	u = _load_universe(pdb_path)
	nmodels = len(u.trajectory)
	if not (1 <= model_number <= nmodels):
		raise ValueError(f"model_number out of range: {model_number} (1..{nmodels})")
	
	# convert to 0-based for MDAnalysis
	u.trajectory[model_number - 1]

	protein_atoms = _get_protein_atom_group(u)
	return _list_chain_ids_from_atom_group(protein_atoms)
# -----------------------------------------------------------------------------------------------------
# Bonded C(i)-N(i+1) distances stay below 1.8 A (1.33 A ideal) in the NMR and X-ray test sets,
# while chain breaks hidden by consecutive numbering start at 2.3 A.
MAX_PEPTIDE_BOND_LENGTH = 2.0


def _require_peptide_bonds(selection: AtomGroup, chain_id: str) -> None:
	"""Reject a chain whose consecutively numbered residues are not bonded.

	Consecutive numbering does not guarantee a continuous chain: residues missing from the
	model, or insertion-coded residues discarded by ``insertion_code_policy``, leave two
	residues with consecutive numbers far apart. Their C(i)-N(i+1) distance shows it.
	"""
	names = [str(name).strip().upper() for name in selection.names]
	resids = [int(resid) for resid in selection.resids]
	positions = np.asarray(selection.positions, dtype=float)
	carbon = {resid: positions[index] for index, (resid, name) in enumerate(zip(resids, names)) if name == "C"}
	nitrogen = {resid: positions[index] for index, (resid, name) in enumerate(zip(resids, names)) if name == "N"}
	breaks = []
	for resid in sorted(carbon):
		if resid + 1 in nitrogen:
			distance = float(np.linalg.norm(carbon[resid] - nitrogen[resid + 1]))
			if distance > MAX_PEPTIDE_BOND_LENGTH:
				breaks.append(f"{resid}->{resid + 1} ({distance:.2f} A)")
	if breaks:
		raise ValueError(
			f"chain {chain_id} is broken between residues {', '.join(breaks)}: the C-N distance exceeds "
			f"{MAX_PEPTIDE_BOND_LENGTH:.1f} A, so these consecutively numbered residues are not bonded "
			"(for example, residues missing from the model, or insertion-coded residues discarded by "
			"insertion_code_policy='prefer_blank')."
		)
# -----------------------------------------------------------------------------------------------------
def extract_model_chain(
	pdb_path: str | Path,
	model_number: int,  # 1-based
	chain_id: str,
	*,
	allow_gaps: bool = False,
	drop_nonpositive: bool = True,
	insertion_code_policy: str = "prefer_blank",
) -> AtomGroup:
	"""
	Extract atoms from a specified model and chain, excluding heteroatoms and water.
	Optionally verify residue continuity.
	
	Raises
	------
	ValueError
		If the chain is not found or gaps exist (and allow_gaps=False).
	"""
	chain_id = _normalize_chain_id(chain_id)
	u = _load_universe(pdb_path)
	nmodels = len(u.trajectory)
	if not (1 <= model_number <= nmodels):
		raise ValueError(f"model_number out of range: {model_number} (1..{nmodels})")
	u.trajectory[model_number - 1]

	protein_atoms = _get_protein_atom_group(u)
	selection = _select_chain_atoms(protein_atoms, chain_id)
	
	if selection.n_atoms == 0:
		available_chains = _list_chain_ids_from_atom_group(protein_atoms)
		raise ValueError(
			f"Chain {chain_id} not found or contains no protein atoms in {pdb_path}. "
			f"Available chain identifiers: {', '.join(available_chains)}"
		)
	
	if insertion_code_policy != "prefer_blank":
		raise ValueError("insertion_code_policy must be 'prefer_blank'.")

	# Positional parser outputs cannot preserve insertion codes. Keep one
	# deterministic residue for each numeric residue number: prefer the blank
	# code, otherwise keep the first code alphabetically.
	icodes = (
		np.asarray(selection.icodes, dtype=str)
		if hasattr(selection, "icodes")
		else np.full(selection.n_atoms, "", dtype=str)
	)
	residue_codes: dict[int, set[str]] = {}
	for resid, icode in zip(selection.resids, icodes, strict=True):
		residue_codes.setdefault(int(resid), set()).add(str(icode).strip())
	chosen_codes = {
		resid: ("" if "" in codes else sorted(codes)[0])
		for resid, codes in residue_codes.items()
	}
	discarded = {
		resid: sorted(code for code in codes if code != chosen_codes[resid])
		for resid, codes in residue_codes.items()
		if any(code != chosen_codes[resid] for code in codes)
	}
	if discarded:
		keep_mask = np.array(
			[
				str(icode).strip() == chosen_codes[int(resid)]
				for resid, icode in zip(selection.resids, icodes, strict=True)
			],
			dtype=bool,
		)
		selection = selection[keep_mask]
		details = ", ".join(
			f"{resid} (discarded: {','.join(code or '<blank>' for code in codes)})"
			for resid, codes in sorted(discarded.items())
		)
		print_warning(
			f"insertion_code_policy='prefer_blank': chain {chain_id}; "
			f"kept one residue per numeric ID and discarded insertion-coded residues: {details}."
		)

	# Deterministic altloc policy: blank is preferred, otherwise the
	# lexicographically first conformer is retained for each residue/atom.
	keep_indices: list[int] = []
	chosen: dict[tuple[int, str], tuple[str, int]] = {}
	for index, atom in enumerate(selection):
		altloc = str(selection.altLocs[index]).strip() if hasattr(selection, "altLocs") else ""
		key = (int(atom.resid), str(atom.name).strip().upper())
		candidate = (altloc, index)
		previous = chosen.get(key)
		if previous is None or (altloc == "" and previous[0] != "") or (altloc != "" and previous[0] != "" and altloc < previous[0]):
			chosen[key] = candidate
	keep_indices = sorted(index for _, index in chosen.values())
	selection = selection[keep_indices]
	resids = np.sort(np.unique(selection.residues.resids))
	if drop_nonpositive:
		nonpositive = resids[resids <= 0]
		if len(nonpositive):
			selection = selection[np.asarray(selection.resids, dtype=int) >= 1]
			resids = np.sort(np.unique(selection.residues.resids))
			if len(resids) == 0:
				raise ValueError(
					f"chain {chain_id} contains no residues with resid >= 1 after "
					"applying drop_nonpositive=true."
				)
			print_warning(
				f"drop_nonpositive=true: removed {len(nonpositive)} "
				f"residue IDs <= 0 from chain {chain_id}; retaining residues "
				f"{int(resids[0])}..{int(resids[-1])}."
			)
	
	if not allow_gaps:
		validate_residue_sequence(resids.tolist(), context=f"chain {chain_id}")
		_require_peptide_bonds(selection, chain_id)
	
	return selection
# -----------------------------------------------------------------------------------------------------

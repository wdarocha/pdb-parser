"""PDB I/O and model/chain utilities (MDAnalysis-based)."""

from __future__ import annotations
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import requests
import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup

from pdb_parser.utils.fs import ensure_dir

_PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")

# -----------------------------------------------------------------------------------------------------
def _normalize_pdb_id(pdb_id: str) -> str:
	"""Return normalized (uppercased) PDB id and validate format."""
	pid = pdb_id.strip().upper()
	if not _PDB_ID_RE.match(pid):
		raise ValueError(f"Invalid PDB id '{pdb_id}'. Expected 4 alphanumeric chars.")
	return pid
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
	filename = out_dir/f"{pid}.pdb"
	if filename.exists():
		return filename
	
	url = f"https://files.rcsb.org/download/{pid}.pdb"
	try:
		resp = requests.get(url, timeout=(5, 30))
		resp.raise_for_status()
	except requests.RequestException as exc:
		raise ValueError(f"Failed to download PDB '{pid}' from RCSB: {exc}") from exc
	
	# Write as text (PDB is ASCII)
	filename.write_text(resp.text, encoding="utf-8")
	
	return filename
# -----------------------------------------------------------------------------------------------------
def is_nmr_structure(pdb_path: str | Path, max_lines: int = 400) -> bool:
	"""Return True if EXPDTA indicates NMR; scans only the first lines for speed."""
	p = Path(pdb_path)
	with p.open("r", encoding="utf-8", errors="replace") as fh:
		for i, line in enumerate(fh):
			if i > max_lines:
				break
			if line.startswith("EXPDTA") and "NMR" in line.upper():
				return True
	return False
# -----------------------------------------------------------------------------------------------------
def list_number_of_models(pdb_path) -> int:
	"""
	Return the number of models available in the PDB file.
	"""
	u = mda.Universe(pdb_path, multiframe=True)
	return len(u.trajectory)
# -----------------------------------------------------------------------------------------------------
def list_chains_for_model(pdb_path: str | Path, model_number: int) -> list[str]:
	"""Return a sorted list of chain identifiers for the given 1-based model_number."""
	u = mda.Universe(str(pdb_path), multiframe=True)
	nmodels = len(u.trajectory)
	if not (1 <= model_number <= nmodels):
		raise ValueError(f"model_number out of range: {model_number} (1..{nmodels})")
	
	# convert to 0-based for MDAnalysis
	u.trajectory[model_number - 1]
	
	if hasattr(u.atoms, "segids"):
		chains = set(u.atoms.segids)
	elif hasattr(u.atoms, "chainIDs"):
		chains = set(u.atoms.chainIDs)
	else:
		chains = set(atom.segid for atom in u.atoms)
	
	return sorted(c for c in chains if c is not None and str(c).strip() != "")
# -----------------------------------------------------------------------------------------------------
def extract_model_chain(
	pdb_path: str | Path,
	model_number: int,  # 1-based
	chain_id: str,
	*,
	allow_gaps: bool = False,
) -> AtomGroup:
	"""
	Extract atoms from a specified model and chain, excluding heteroatoms and water.
	Optionally verify residue continuity.
	
	Raises
	------
	ValueError
		If the chain is not found or gaps exist (and allow_gaps=False).
	"""
	u = mda.Universe(str(pdb_path), multiframe=True)
	u.trajectory[model_number - 1]
	
	# Select only standard protein atoms
	selection = u.select_atoms(f"segid {chain_id} and protein")
	
	if selection.n_atoms == 0:
		raise ValueError(f"Chain {chain_id} not found or contains no protein atoms in {pdb_path}")
	
	resids = np.sort(np.unique(selection.residues.resids))
	
	if not allow_gaps:
		diffs = np.diff(resids)
		if not np.all(diffs == 1):
			gaps = resids[np.where(diffs > 1)[0] + 1].tolist()
			raise ValueError(f"Gaps found in residue sequence at residues: {gaps}")
	
	return selection
# -----------------------------------------------------------------------------------------------------

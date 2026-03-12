"""High-level validators for PDB structures (model/chain readiness).

This module centralizes the common checks you do before processing a structure:
- experimental method must be NMR (based on EXPDTA)
- chosen model must exist (1-based)
- chosen chain must be present in that model
- optional continuity check is delegated to extract_model_chain()

If all checks pass, it returns the extracted AtomGroup.
Otherwise, it raises ValueError with an informative message.
"""

from __future__ import annotations
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup

from pdb_parser.utils import *
from pdb_parser.io.pdb_ops import *
# -----------------------------------------------------------------------------------------------------
def ensure_nmr_model_chain_ready(
	pdb_path: str | Path,
	model_number: int, # 1-based
	chain_id: str,
	*,
	allow_gaps: bool = False,
) -> AtomGroup:
	"""Validate NMR method, model range, and chain availability; then extract AtomGroup.
	
	Parameters
	----------
	pdb_path : str | Path
		Path to the PDB file (text .pdb).
	model_number : int
		1-based model index to use.
	chain_id : str
		Chain (segid) to extract from the selected model.
	allow_gaps : bool
		If False, raise if residue id sequence has gaps. If True, ignore gaps.
	
	Returns
	-------
	AtomGroup
		The extracted selection (protein atoms for the given model/chain).
	
	Raises
	------
	ValueError
		If the structure is not NMR, model_number is out of range, chain is absent,
		or residue gaps are present (when allow_gaps=False).
	"""
	p = Path(pdb_path)
	
	# 1) Experimental method must be NMR
	if not is_nmr_structure(p):
		raise ValueError("structure was not determined by NMR")
	
	# 2) Model range check (1-based)
	num_models = list_number_of_models(p)
	if model_number < 1 or model_number > num_models:
		raise ValueError(f"requested model {model_number} is out of range (valid: 1–{num_models})")
	
	# 3) Chain availability in that model
	chains = list_chains_for_model(p, model_number)
	if chain_id not in chains:
		raise ValueError(
			f"model {model_number} does not contain chain '{chain_id}'. "
			f"available chains: {', '.join(chains)}"
		)
	
	# 4) Extract (this also validates gaps if allow_gaps=False)
	return extract_model_chain(p, model_number, chain_id, allow_gaps=allow_gaps)
# -----------------------------------------------------------------------------------------------------

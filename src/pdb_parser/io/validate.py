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
from typing import TYPE_CHECKING

if TYPE_CHECKING:
	from MDAnalysis.core.groups import AtomGroup

from pdb_parser.io.pdb_ops import (
	detect_experimental_method,
	extract_model_chain,
	list_chains_for_model,
	list_number_of_models,
)
from pdb_parser.messages import print_warning


def ensure_model_chain_ready(
	pdb_path: str | Path,
	model_number: int,
	chain_id: str,
	*,
	method: str,
	allow_gaps: bool = False,
	drop_nonpositive: bool = True,
) -> tuple[AtomGroup, int, str]:
	"""Validate method/model/chain and return atoms, effective model and method."""
	p = Path(pdb_path)
	detected_method = detect_experimental_method(p)
	requested_method = method.strip().lower()
	if requested_method not in {"auto", "nmr", "xray"}:
		raise ValueError("method must be one of: auto, nmr, xray")
	if requested_method != "auto" and requested_method != detected_method:
		raise ValueError(
			f"Requested method '{requested_method}' does not match PDB EXPDTA method '{detected_method}'."
		)

	num_models = list_number_of_models(p)
	effective_model = model_number
	if detected_method == "xray" and num_models == 1 and model_number != 1:
		effective_model = 1
		print_warning(
			f"X-ray structure contains one coordinate model; requested model {model_number} was replaced by model 1."
		)
	elif model_number < 1 or model_number > num_models:
		raise ValueError(f"requested model {model_number} is out of range (valid: 1–{num_models})")

	chains = list_chains_for_model(p, effective_model)
	if chain_id not in chains:
		raise ValueError(
			f"model {effective_model} does not contain chain '{chain_id}'. available chains: {', '.join(chains)}"
		)
	return (
		extract_model_chain(
			p,
			effective_model,
			chain_id,
			allow_gaps=allow_gaps,
			drop_nonpositive=drop_nonpositive,
		),
		effective_model,
		detected_method,
	)


def ensure_nmr_model_chain_ready(
	pdb_path: str | Path,
	model_number: int,  # 1-based
	chain_id: str,
	*,
	allow_gaps: bool = False,
	drop_nonpositive: bool = True,
) -> AtomGroup:
	"""Validate NMR method, model range, and chain availability; then extract AtomGroup.

	Parameters
	----------
	pdb_path : str | Path
		Path to the PDB file (text .pdb).
	model_number : int
		1-based model index to use.
	chain_id : str
		Chain identifier (segid or chainID) to extract from the selected model.
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
	try:
		is_nmr = detect_experimental_method(Path(pdb_path)) == "nmr"
	except ValueError:
		is_nmr = False
	if not is_nmr:
		raise ValueError("structure was not determined by NMR")
	atoms, _, _ = ensure_model_chain_ready(
		pdb_path,
		model_number,
		chain_id,
		method="nmr",
		allow_gaps=allow_gaps,
		drop_nonpositive=drop_nonpositive,
	)
	return atoms

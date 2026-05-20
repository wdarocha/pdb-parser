"""Input readers for pdb-parser: PDB id list and params file."""

from __future__ import annotations
from pathlib import Path
from typing import Dict, List
import re
import sys

_PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")
# -----------------------------------------------------------------------------------------------------
def read_pdb_ids(path: Path) -> List[str]:
	"""Read PDB ids from a text file, validate, and deduplicate preserving order.
	
	Rules:
	- Lines starting with '#' are comments.
	- Inline comments after '#' are removed.
	- Blank lines are ignored.
	- IDs are uppercased and must match ^[A-Za-z0-9]{4}$.
	"""
	if not path.exists():
		raise FileNotFoundError(f"PDB id list not found: {path}")
	
	seen: set[str] = set()
	result: List[str] = []
	
	with path.open("r", encoding="utf-8", errors="replace") as fh:
		for lineno, raw in enumerate(fh, start=1):
			line = raw.strip()
			if not line or line.startswith("#"):
				continue
			if "#" in line:
				line = line.split("#", 1)[0].strip()
			if not line:
				continue
			
			pid = line.upper()
			if not _PDB_ID_RE.match(pid):
				raise ValueError(
					f"Invalid PDB id at {path}:{lineno} -> '{line}'. "
					"Expected exactly 4 alphanumeric characters."
				)
			if pid not in seen:
				seen.add(pid)
				result.append(pid)
		
	if not result:
		raise ValueError(f"No valid PDB ids found in {path}.")
	return result
# -----------------------------------------------------------------------------------------------------
def read_params(path: Path) -> Dict[str, str]:
	"""Read a key: value parameter file.
	
	Rules
	- Lines starting with '#' are comments.
	- Inline comments after '#' are removed.
	- Keys are lowercased and stripped; values are stripped (kept as strings).
	- Empty lines are ignored.
	- Only the ':' separator is accepted.
	"""
	if not path.exists():
		raise FileNotFoundError(f"Params file not found: {path}")
	
	params: Dict[str, str] = {}
	with path.open("r", encoding="utf-8", errors="replace") as fh:
		for lineno, raw in enumerate(fh, start=1):
			line = raw.strip()
			if not line or line.startswith("#"):
				continue
			if "#" in line:
				line = line.split("#", 1)[0].strip()
			if not line:
				continue
			if ":" not in line:
				raise ValueError(
					f"Invalid params line at {path}:{lineno} -> '{raw.rstrip()}'. "
					"Expected 'key: value'."
				)
			key, _, value = line.partition(":")
			k = key.strip().lower()
			v = value.strip()
			if not k:
				raise ValueError(f"Empty key at {path}:{lineno}.")
			params[k] = v
		
	if not params:
		print(f"[WARN] Params file {path} parsed empty.", file=sys.stderr)
	return params
# -----------------------------------------------------------------------------------------------------

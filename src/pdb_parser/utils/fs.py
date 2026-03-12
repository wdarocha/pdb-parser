"""Filesystem utilities."""

from __future__ import annotations
from pathlib import Path

# -----------------------------------------------------------------------------------------------------
def ensure_dir(path: str | Path) -> Path:
	"""Create directory if it does not exist; return Path.
	
	This is idempotent and safe to call repeatedly.
	"""
	p = Path(path)
	p.mkdir(parents=True, exist_ok=True)
	return p
# -----------------------------------------------------------------------------------------------------


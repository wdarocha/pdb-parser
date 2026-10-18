"""Filesystem utilities."""

from __future__ import annotations

from pathlib import Path


def ensure_dir(path: str | Path) -> Path:
	"""Create directory if it does not exist; return Path.

	This is idempotent and safe to call repeatedly.
	"""
	p = Path(path)
	p.mkdir(parents=True, exist_ok=True)
	return p


def ensure_existing_dir(path: str | Path, label: str = "Directory") -> Path:
	"""Require one existing directory path and return it as `Path`."""
	directory = Path(path)
	if not directory.exists():
		raise ValueError(f"{label} does not exist: {directory}")
	if not directory.is_dir():
		raise ValueError(f"{label} is not a directory: {directory}")
	return directory

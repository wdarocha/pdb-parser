"""Shared user-facing message helpers."""

from __future__ import annotations

import sys


def print_ok(message: str) -> None:
	"""Print a user-facing success message to stdout."""
	print(f"[OK] {message}")


def print_warning(message: str) -> None:
	"""Print a user-facing warning message to stderr."""
	print(f"[WARN] {message}", file=sys.stderr)


def print_error(message: str) -> None:
	"""Print a user-facing error message to stderr."""
	print(f"[ERROR] {message}", file=sys.stderr)


def print_skip(pdb_id: str, message: str) -> None:
	"""Print the standard message for a skipped structure."""
	print(f"[{pdb_id}] Skipped: {message}")

"""Shared validation of the residue sequence selected for one execution."""

from __future__ import annotations

from collections.abc import Iterable


def residue_sequence_from_ids(residue_ids: Iterable[int], *, context: str) -> list[int]:
	"""Validate an unordered collection of residue IDs and return sorted IDs."""
	return validate_residue_sequence(sorted({int(value) for value in residue_ids}), context=context)


def are_adjacent_residues(
	residue_ids: Iterable[int], left: int, right: int, *, context: str
) -> bool:
	"""Return whether ``right`` is the immediate successor of ``left``."""
	ids = residue_sequence_from_ids(residue_ids, context=context)
	try:
		return ids[ids.index(int(left)) + 1] == int(right)
	except (ValueError, IndexError):
		return False


def validate_residue_sequence(
	residue_ids: Iterable[int], *, context: str = "selected chain"
) -> list[int]:
	"""Return residue IDs after requiring strict ``n, n+1, ...`` continuity."""
	ids = [int(value) for value in residue_ids]
	if not ids:
		raise ValueError(f"{context} contains no residues.")
	if any(value <= 0 for value in ids):
		raise ValueError(f"{context} contains non-positive residue IDs: {ids!r}.")
	if len(set(ids)) != len(ids):
		raise ValueError(f"{context} contains repeated residue IDs: {ids!r}.")
	gaps = [(left, right) for left, right in zip(ids, ids[1:]) if right != left + 1]
	if gaps:
		details = ", ".join(f"{left}->{right}" for left, right in gaps)
		raise ValueError(
			f"{context} residue sequence must be strictly consecutive; "
			f"gaps are delimited by residues {details}."
		)
	return ids

"""First-residue clique initialization."""

from __future__ import annotations

from .contracts import AtomRefList, CliqueAtomRow
from .ordering import first_residue_name_order

Residue1Row = tuple[CliqueAtomRow, int]


def _row(names: list[str | None], residue_id: int, row_id: int) -> Residue1Row:
    if len(names) != 4:
        raise ValueError(f"atom_names must have length 4, got {len(names)}.")
    return ([(residue_id, name) if name is not None else None for name in names], row_id)


def build_first_residue_pattern(available_atoms: AtomRefList, res1_id: int) -> list[Residue1Row]:
	"""Build the five-, six-, or seven-atom initial clique pattern of residue 1.

	The atoms come from ``first_residue_name_order``, the same order used for the remapped X
	file, so every N-terminal pattern (for example ``H3 h H1`` or ``h H2 H1``) gets matching
	cliques: n1, n2, n3 are the N-terminal hydrogens in that order.
	"""
	names = {name for residue_id, name in available_atoms if residue_id == res1_id}
	if len(names) not in (5, 6, 7):
		raise ValueError(f"Unsupported number of atoms for residue 1: {len(names)}. Expected 5, 6, or 7. Found atoms: {sorted(names)}.")
	order = first_residue_name_order(names)
	if set(order) != names:
		raise ValueError(f"Residue 1 atoms {sorted(names)} do not match the first-residue pattern {order}.")
	if len(order) == 7:
		n1, n2, n3, _, _, ha, _ = order
		rows = [(None, None, None, n1, 0), (None, None, n1, n2, 0), (None, n1, n2, n3, 0), (n1, n2, n3, "N", 1), (n2, n3, "N", "CA", 1), (n3, "N", "CA", ha, 2), ("N", "CA", ha, "C", 1)]
	elif len(order) == 6:
		n1, n2, _, _, ha, _ = order
		rows = [(None, None, None, n1, 0), (None, None, n1, n2, 0), (None, n1, n2, "N", 0), (n1, n2, "N", "CA", 1), (n2, "N", "CA", ha, 2), ("N", "CA", ha, "C", 1)]
	else:
		_, _, _, ha, h = order
		rows = [(None, None, None, "N", 0), (None, None, "N", "CA", 0), (None, "N", "CA", "C", 0), ("N", "CA", "C", ha, 1), (ha, "CA", "N", h, 2)]
	return [_row(list(row[:4]), res1_id, row[4]) for row in rows]

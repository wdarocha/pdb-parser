"""Clique-pattern definitions and atom-name resolution."""

from __future__ import annotations

from .contracts import AtomRefList

CLIQUE_PATTERN_MAP: dict[int, tuple[tuple[int, str], ...]] = {
    1: ((-1, "N"), (-1, "CA"), (-1, "C"), (0, "N")),
    2: ((-1, "CA"), (-1, "C"), (0, "N"), (0, "CA")),
    3: ((-1, "C"), (0, "N"), (0, "CA"), (0, "C")),
    4: ((-1, "C"), (0, "N"), (0, "CA"), (0, "HN")),
    5: ((0, "N"), (0, "CA"), (0, "C"), (0, "HA")),
    6: ((0, "HN"), (0, "N"), (0, "CA"), (0, "HA")),
    7: ((0, "N"), (0, "CA"), (0, "HA"), (0, "C")),
    8: ((-1, "CA"), (-1, "C"), (0, "N"), (0, "HN")),
    9: ((-1, "C"), (0, "N"), (0, "HN"), (0, "CA")),
    10: ((-1, "HA"), (-1, "CA"), (-1, "C"), (0, "HN")),
    11: ((-1, "CA"), (-1, "C"), (0, "HN"), (0, "N")),
    12: ((-1, "CA"), (-1, "C"), (0, "HN"), (0, "CA")),
    13: ((-1, "C"), (0, "HN"), (0, "CA"), (0, "N")),
}
DDGP_ORDER_PATTERN_MAP: dict[int, tuple[int, int, int, int, int]] = {
    1: (1, 2, 3, 4, 5), 2: (1, 2, 3, 5, 4), 3: (1, 2, 4, 6, 7),
    4: (1, 2, 4, 3, 5), 5: (1, 8, 9, 6, 7), 6: (1, 8, 9, 3, 5),
    7: (10, 11, 9, 6, 7), 8: (10, 11, 9, 3, 5),
    9: (10, 12, 13, 6, 7), 10: (10, 12, 13, 3, 5),
}
ATOM_ALIAS_MAP = {"HN": ("H1", "H", "HD2", "HD3"), "HA": ("HA", "HA2", "HA3")}


def build_available_atoms(df_Xremapped) -> AtomRefList:
    return [(int(row[2]), str(row[1])) for row in df_Xremapped.itertuples(index=False)]


def get_atoms_of_residue(residue_id: int, available_atoms: AtomRefList) -> set[str]:
    return {name for current_id, name in available_atoms if current_id == residue_id}


def resolve_atom_name(atom_name: str, residue_id: int, available_atoms: AtomRefList | None = None) -> str:
    if available_atoms is None or atom_name not in ATOM_ALIAS_MAP:
        return atom_name
    residue_atoms = get_atoms_of_residue(residue_id, available_atoms)
    return next((name for name in ATOM_ALIAS_MAP[atom_name] if name in residue_atoms), atom_name)


def get_clique_pattern_ids(ddgp_pattern_id: int) -> tuple[int, int, int, int, int]:
    try:
        return DDGP_ORDER_PATTERN_MAP[ddgp_pattern_id]
    except KeyError as error:
        raise ValueError(f"DDGP order pattern {ddgp_pattern_id} is not implemented. Valid values are 1 to 10.") from error


def map_clique_type(pattern_id: int) -> int:
    if pattern_id in {2, 4, 5, 7, 8, 9, 11, 12, 13}:
        return 1
    if pattern_id in {6, 10}:
        return 2
    if pattern_id in {1, 3}:
        return 3
    raise ValueError(f"Invalid pattern id: {pattern_id}")


def build_clique_pairs(residue_id: int, clique_pattern_id: int, available_atoms: AtomRefList | None = None) -> AtomRefList:
    try:
        pattern = CLIQUE_PATTERN_MAP[clique_pattern_id]
    except KeyError as error:
        raise ValueError(f"Clique pattern {clique_pattern_id} is not implemented. Valid values are 1 to 13.") from error
    return [(residue_id + offset, resolve_atom_name(name, residue_id + offset, available_atoms)) for offset, name in pattern]


def build_ddgp_pattern_entries(ddgp_pattern_id: int, residue_id: int, available_atoms: AtomRefList | None = None) -> list[tuple[AtomRefList, int]]:
    return [(build_clique_pairs(residue_id, pattern_id, available_atoms), map_clique_type(pattern_id)) for pattern_id in get_clique_pattern_ids(ddgp_pattern_id)]

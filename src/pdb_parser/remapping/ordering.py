"""Atom-order construction helpers."""

from __future__ import annotations

import pandas as pd

from .contracts import AtomRef, AtomRefList


def get_nonterminal_h_group(resname: str) -> tuple[str, ...]:
    """
    Return the allowed atom names for the backbone-N hydrogen position
    in non-terminal residues.
    """
    resname = resname.strip().upper()

    if resname == "PRO":
        return ("HD2", "HD3")

    return ("H",)


def get_ha_group(resname: str) -> tuple[str, ...]:
    """
    Return the allowed atom names for the alpha-hydrogen position.
    """
    resname = resname.strip().upper()

    if resname == "GLY":
        return ("HA2", "HA3")

    return ("HA",)


def get_internal_residue_ordering_template(
    ordering_id: int,
    resname: str,
) -> list[tuple[str, ...]]:
    """
    Return the concrete atom ordering template for a non-terminal residue.

    Each position is represented as a tuple of allowed atom names.
    For fixed atoms, the tuple has length 1.
    For ambiguous positions, the tuple contains all accepted alternatives.
    """
    h_group = get_nonterminal_h_group(resname)
    ha_group = get_ha_group(resname)

    ordering_templates: dict[int, list[tuple[str, ...]]] = {
        1:  [("N",), ("CA",), ("C",), h_group, ha_group],
        2:  [("N",), ("CA",), ("C",), ha_group, h_group],
        3:  [("N",), ("CA",), h_group, ha_group, ("C",)],
        4:  [("N",), ("CA",), h_group, ("C",), ha_group],
        5:  [("N",), h_group, ("CA",), ha_group, ("C",)],
        6:  [("N",), h_group, ("CA",), ("C",), ha_group],
        7:  [h_group, ("N",), ("CA",), ha_group, ("C",)],
        8:  [h_group, ("N",), ("CA",), ("C",), ha_group],
        9:  [h_group, ("CA",), ("N",), ha_group, ("C",)],
        10: [h_group, ("CA",), ("N",), ("C",), ha_group],
    }

    if ordering_id not in ordering_templates:
        raise ValueError(
            f"Invalid order_id: {ordering_id}. Expected an integer from 1 to 10."
        )

    return ordering_templates[ordering_id]


def resolve_atom_name_choice(
    atom_name_to_atom_id: dict[str, int],
    allowed_atom_names: tuple[str, ...],
    residue_id: int,
    resname: str,
) -> AtomRef:
    """
    Resolve one atom position by selecting the first allowed atom name
    that is present in the residue.
    """
    for atom_name in allowed_atom_names:
        atom_id = atom_name_to_atom_id.get(atom_name)
        if atom_id is not None:
            return atom_id, atom_name

    raise ValueError(f"For option 'backbone_plus_hydrogens' in residue {resname}_{residue_id} allowed atom names {allowed_atom_names} were not found")


def get_numeric_atom_order_from_named_order(
    df_X_i: pd.DataFrame,
    ordered_atom_name_groups: list[tuple[str, ...]],
) -> AtomRefList:
    """
    Return the ordered list of (atom_id, atom_name) pairs for one residue
    from an ordered list of allowed atom-name groups.

    The input DataFrame must contain exactly one residue and follow the
    column convention:
    - column 0: atom_id
    - column 1: atom_name
    - column 2: residue_id
    - column 3: residue_name
    """
    if df_X_i.empty:
        raise ValueError("The residue DataFrame is empty.")

    residue_id = int(df_X_i.iat[0, 2])
    resname = str(df_X_i.iat[0, 3]).strip().upper()

    atom_name_to_atom_id = {
        str(atom_name).strip().upper(): int(atom_id)
        for atom_id, atom_name in zip(df_X_i[0], df_X_i[1], strict=False)
    }

    resolved_order: AtomRefList = []

    for allowed_atom_names in ordered_atom_name_groups:
        atom_id, atom_name = resolve_atom_name_choice(
            atom_name_to_atom_id,
            allowed_atom_names,
            residue_id,
            resname,
        )
        resolved_order.append((atom_id, atom_name))

    return resolved_order


def get_internal_residue_numeric_order(
    df_X_i: pd.DataFrame,
    ordering_id: int,
) -> AtomRefList:
    """
    Return the ordered list of (atom_id, atom_name) pairs for one
    non-terminal residue.
    """
    if df_X_i.empty:
        raise ValueError("The residue DataFrame is empty.")

    resname = str(df_X_i.iat[0, 3]).strip().upper()

    ordered_atom_name_groups = get_internal_residue_ordering_template(ordering_id, resname)

    return get_numeric_atom_order_from_named_order(df_X_i, ordered_atom_name_groups)


def get_backbone_residue_numeric_order(
    df_X_i: pd.DataFrame,
) -> AtomRefList:
    """Return the fixed backbone-only order `(N, CA, C)` for one residue."""
    return get_numeric_atom_order_from_named_order(
        df_X_i,
        [("N",), ("CA",), ("C",)],
    )


def choose_first_existing_atom(
    atom_names_present: set[str],
    candidates: tuple[str, ...],
    excluded_atoms: set[str] | None = None,
) -> str | None:
    """
    Return the first candidate atom name that exists in the residue and is
    not excluded. Return None if no candidate is available.
    """
    if excluded_atoms is None:
        excluded_atoms = set()

    for atom_name in candidates:
        if atom_name in atom_names_present and atom_name not in excluded_atoms:
            return atom_name

    return None

def get_first_residue_ha_atom(
    atom_names_present: set[str],
) -> str:
    """
    Return the alpha-hydrogen atom to be used in the first residue.
    """
    ha_atom = choose_first_existing_atom(
        atom_names_present,
        ("HA", "HA2", "HA3"),
    )

    if ha_atom is None:
        raise ValueError("The first residue does not contain any HA/HA2/HA3 atom.")

    return ha_atom


def first_residue_name_order(atom_names_present: set[str]) -> list[str]:
    """
    Return the atom names of the first residue in its remapping order.

    Supported patterns
    ------------------
    1. H3, H2, h, N, CA, ha, C
    2. H3, h, H1, N, CA, ha, C
    3. h, H2, H1, N, CA, ha, C
    4. H2, h, N, CA, ha, C
    5. H3, h, N, CA, ha, C
    6. N, CA, C, ha, h

    where h can be H, H1, HD2, or HD3, and ha can be HA, HA2, or HA3.
    """
    required_backbone = {"N", "CA", "C"}
    missing_backbone = sorted(required_backbone - atom_names_present)

    if missing_backbone:
        raise ValueError(f"The first residue is missing required backbone atom(s): {missing_backbone}.")

    ha_atom = get_first_residue_ha_atom(atom_names_present)

    h_candidates = ("H", "H1", "HD2", "HD3")
    patterns: list[tuple[tuple[str, ...], tuple[str | None, ...]]] = [
        (("H3", "H2"), ("H3", "H2", None, "N", "CA", ha_atom, "C")),
        (("H3", "H1"), ("H3", None, "H1", "N", "CA", ha_atom, "C")),
        (("H2", "H1"), (None, "H2", "H1", "N", "CA", ha_atom, "C")),
        (("H2",), ("H2", None, "N", "CA", ha_atom, "C")),
        (("H3",), ("H3", None, "N", "CA", ha_atom, "C")),
    ]

    for required_atoms, atom_layout in patterns:
        if not set(required_atoms).issubset(atom_names_present):
            continue

        h_atom = choose_first_existing_atom(
            atom_names_present,
            h_candidates,
            excluded_atoms=set(required_atoms),
        )
        if h_atom is None:
            continue

        return [h_atom if atom_name is None else atom_name for atom_name in atom_layout]

    h_atom = choose_first_existing_atom(atom_names_present, h_candidates)
    if h_atom is not None:
        return ["N", "CA", "C", ha_atom, h_atom]

    raise ValueError("Could not match any valid ordering pattern for the first residue.")


def first_residue_order(
    df_X_i: pd.DataFrame,
    resnum_1: int,
) -> AtomRefList:
    """
    Return the ordered list of (atom_id, atom_name) pairs for the first residue,
    following ``first_residue_name_order``.
    """
    if df_X_i.empty:
        raise ValueError("The first-residue DataFrame is empty.")

    atom_name_to_atom_id = {
        str(atom_name).strip().upper(): int(atom_id)
        for atom_id, atom_name in zip(df_X_i[0], df_X_i[1], strict=False)
    }
    return [
        (atom_name_to_atom_id[atom_name], atom_name)
        for atom_name in first_residue_name_order(set(atom_name_to_atom_id))
    ]

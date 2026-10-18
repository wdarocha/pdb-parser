"""Validation helpers for supported remapping atom layouts."""

from __future__ import annotations

import pandas as pd

from pdb_parser.messages import print_skip

from .contracts import RemapAtomLayout


def get_missing_backbone_atoms(atom_names: set[str]) -> list[str]:
    """Return the missing atoms required by the backbone-only public layout."""
    return sorted({"N", "CA", "C"} - atom_names)

def get_missing_backbone_plus_hydrogen_atoms(
    atom_names: set[str],
    resname: str,
    resnum: int,
    resnum_1: int,
) -> list[str]:
    """
    Return the list of missing required atoms for one residue.
    """
    resname = resname.strip().upper()
    missing_atoms: list[str] = []

    # --- backbone atoms always required
    for atom_name in ("N", "CA", "C"):
        if atom_name not in atom_names:
            missing_atoms.append(atom_name)

    # --- backbone nitrogen hydrogen requirement
    if resname != "PRO":
        if resnum == resnum_1:
            if not any(h in atom_names for h in ("H", "H1", "H2", "H3")):
                missing_atoms.append("H/H1/H2/H3")
        else:
            if "H" not in atom_names:
                missing_atoms.append("H")

    # --- residue-specific side/backbone hydrogen requirements
    if resname == "GLY":
        if ("HA2" not in atom_names) and ("HA3" not in atom_names):
            missing_atoms.append("HA2/HA3")

    elif resname == "PRO":
        if "HA" not in atom_names:
            missing_atoms.append("HA")

        if ("HD2" not in atom_names) and ("HD3" not in atom_names):
            missing_atoms.append("HD2/HD3")

    else:
        if "HA" not in atom_names:
            missing_atoms.append("HA")

    return missing_atoms


def validate_backbone_plus_hydrogens_residue(
    df_X_i: pd.DataFrame,
    resnum_i: int,
    pdb_id: str,
    chosen_model: int,
    chosen_chain: str,
    resnum_1: int,
) -> int:
    """
    Validate one residue for the 'backbone_plus_hydrogens' atom selection.

    Returns
    -------
    int
        1 if the structure must be skipped, 0 otherwise.
    """
    if df_X_i.empty:
        print_skip(
            pdb_id,
            f"model {chosen_model}, chain {chosen_chain} does not contain residue {resnum_i}",
        )
        return 1

    atom_names_i = set(df_X_i[1].astype(str).str.strip().str.upper().tolist())

    resname_i = str(df_X_i.iat[0, 3]).strip().upper()

    missing_atoms = get_missing_backbone_plus_hydrogen_atoms(
        atom_names_i, resname_i, resnum_i, resnum_1
    )

    if missing_atoms:
        print_skip(
            pdb_id,
            f"(model {chosen_model}, chain {chosen_chain}) "
            f"- option 'backbone_plus_hydrogens': missing atom(s) {missing_atoms} "
            f"in residue {resnum_i} ({resname_i}).",
        )
        return 1

    return 0


def detect_remap_atom_layout(
    df_X: pd.DataFrame,
    res_ids: list[int],
) -> RemapAtomLayout:
    """Detect which public parser atom layout produced the remap input files."""
    if not res_ids:
        raise ValueError(
            "Remapping requires parser outputs generated with a supported public atom_selection, "
            "but no residues were found in the parser structure file."
        )

    first_residue_id = int(res_ids[0])
    is_backbone_only = True
    is_backbone_plus_hydrogens = True

    for residue_id in res_ids:
        df_X_i = df_X[df_X[2] == residue_id]
        if df_X_i.empty:
            raise ValueError(
                "Remapping requires parser outputs generated with a supported public "
                "atom_selection, "
                f"but residue {residue_id} is missing from the parser structure file."
            )

        atom_names = set(df_X_i[1].astype(str).str.strip().str.upper().tolist())
        resname = str(df_X_i.iat[0, 3]).strip().upper()
        missing_backbone_atoms = get_missing_backbone_atoms(atom_names)
        if missing_backbone_atoms:
            raise ValueError(
                "Remapping requires parser outputs generated with a supported public "
                "atom_selection. "
                f"Residue {resname}_{residue_id} is missing backbone atoms "
                f"{missing_backbone_atoms}."
            )

        if atom_names != {"N", "CA", "C"}:
            is_backbone_only = False

        missing_backbone_plus_hydrogens_atoms = get_missing_backbone_plus_hydrogen_atoms(
            atom_names,
            resname,
            int(residue_id),
            first_residue_id,
        )
        if missing_backbone_plus_hydrogens_atoms:
            is_backbone_plus_hydrogens = False

    if is_backbone_only:
        return "backbone"

    if is_backbone_plus_hydrogens:
        return "backbone_plus_hydrogens"

    raise ValueError(
        "Remapping requires parser outputs generated with atom_selection='backbone' "
        "or atom_selection='backbone_plus_hydrogens'. "
        "The parser structure file does not match either supported public layout."
    )


def ensure_backbone_plus_hydrogens_dataframe(
    df_X: pd.DataFrame,
    res_ids: list[int],
) -> None:
    """Require remap input files to match the public backbone_plus_hydrogens layout."""
    atom_layout = detect_remap_atom_layout(df_X, res_ids)
    if atom_layout != "backbone_plus_hydrogens":
        raise ValueError(
            "Remapping requires parser outputs generated with atom_selection="
            "'backbone_plus_hydrogens'."
        )

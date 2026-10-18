"""Dataframe transformations used by the remapping stage."""

from __future__ import annotations

import pandas as pd


def remap_structure_dataframe(
    df: pd.DataFrame,
    new_order: list[int],
) -> pd.DataFrame:
    """Remap structure rows and renumber atom IDs sequentially."""
    df_remapped = df.set_index(0).loc[new_order].reset_index()
    df_remapped[0] = range(1, len(df_remapped) + 1)
    return df_remapped


def remap_distance_dataframe(
    df_D: pd.DataFrame,
    atom_ids_new_order: list[int],
) -> pd.DataFrame:
    """Remap distance-restraint atom IDs to the supplied 1-based order."""
    df_D_new = df_D.copy()
    old_to_new = {
        old_atom_id: new_atom_id
        for new_atom_id, old_atom_id in enumerate(atom_ids_new_order, start=1)
    }
    remapped_columns = df_D_new.loc[:, [0, 1]].apply(
        lambda column: column.map(old_to_new)
    )
    # Validate source references before assigning mapped values or coercing the
    # columns to integers, so the diagnostic identifies the original atom IDs.
    source_ids = df_D.loc[:, [0, 1]].to_numpy().ravel()
    unmapped = pd.notna(source_ids) & remapped_columns.isna().to_numpy().ravel()
    missing_atom_ids = {
        atom_id.item() if hasattr(atom_id, "item") else atom_id
        for atom_id in source_ids[unmapped]
    }
    if missing_atom_ids:
        raise ValueError(
            "Distance constraints reference atom IDs absent from the remapping "
            f"order: {sorted(missing_atom_ids, key=str)}."
        )
    df_D_new.loc[:, [0, 1]] = remapped_columns
    df_D_new[0] = df_D_new[0].astype(int)
    df_D_new[1] = df_D_new[1].astype(int)

    swap_mask = df_D_new[0] < df_D_new[1]
    if swap_mask.any():
        for left, right in ((0, 1), (2, 3), (6, 7), (8, 9)):
            df_D_new.loc[swap_mask, [left, right]] = df_D_new.loc[
                swap_mask, [right, left]
            ].to_numpy()
    return df_D_new.sort_values(by=[0, 1]).reset_index(drop=True)

"""Torsion matrix helpers for the remapping stage."""

from __future__ import annotations

import numpy as np
import pandas as pd

from pdb_parser.geometry import distances_2_abs_torsion_angle, torsion_angle_with_points

from .contracts import AtomLookup, CliqueMatrixEntry, IntervalBounds, IntervalLookup


def build_interval_lookup(df_I: pd.DataFrame) -> IntervalLookup:
    """Return distance bounds indexed by ordered atom pair."""
    pairs = zip(df_I.iloc[:, 0].astype(int).tolist(), df_I.iloc[:, 1].astype(int).tolist())
    bounds = zip(df_I.iloc[:, 4].astype(float).tolist(), df_I.iloc[:, 5].astype(float).tolist())
    return dict(zip(pairs, bounds))


def get_interval_dl_du(
    atom_1_idx: int, atom_4_idx: int, interval_lookup: IntervalLookup
) -> IntervalBounds:
    """Get bounds for an atom pair; a missing pair is a data issue of the parser I file."""
    try:
        return interval_lookup[(atom_1_idx, atom_4_idx)]
    except KeyError as error:
        raise ValueError(
            f"Interval not found for atom pair ({atom_1_idx}, {atom_4_idx}): the parser I "
            "file has no distance for this discretization edge (for example, the pair is "
            "farther apart than max_distance in the input coordinates)."
        ) from error


def update_torsion_matrix(
    T: np.ndarray,
    df_X: pd.DataFrame,
    df_A: pd.DataFrame | np.ndarray,
    df_I: pd.DataFrame,
) -> np.ndarray:
    """Fill torsion fields in ``T`` in place."""
    coordinates = df_X.iloc[:, 4:7].to_numpy()
    angles = df_A.to_numpy() if hasattr(df_A, "to_numpy") else df_A
    names = {int(row[0]): str(row[1]) for row in df_X.itertuples(index=False)}
    residues = {int(row[0]): int(row[2]) for row in df_X.itertuples(index=False)}
    angle_rows = {int(row[0]): index for index, row in enumerate(angles)}
    intervals = build_interval_lookup(df_I)

    for row_index in range(3, T.shape[0]):
        clique_type = int(T[row_index, 4])
        atom_1, atom_2, atom_3, atom_4 = (
            int(T[row_index, column]) for column in range(4)
        )
        if clique_type == 1:
            tau = torsion_angle_with_points(
                coordinates[atom_4 - 1], coordinates[atom_3 - 1],
                coordinates[atom_2 - 1], coordinates[atom_1 - 1],
            )
            T[row_index, 4:6] = (np.sign(tau), np.abs(tau) * 180.0 / np.pi)
        elif clique_type == 2:
            try:
                d14l, d14u = get_interval_dl_du(atom_1, atom_4, intervals)
            except ValueError as error:
                labels = " and ".join(
                    f"{names.get(atom, '?')}_{residues.get(atom, '?')}" for atom in (atom_1, atom_4)
                )
                raise ValueError(f"{error} Atoms: {labels}.") from error
            x4, x3, x2, x1 = coordinates[
                [atom_1 - 1, atom_2 - 1, atom_3 - 1, atom_4 - 1]
            ]
            d12 = np.linalg.norm(x1 - x2)
            d13 = np.linalg.norm(x1 - x3)
            d23 = np.linalg.norm(x2 - x3)
            d24 = np.linalg.norm(x2 - x4)
            d34 = np.linalg.norm(x3 - x4)
            low = distances_2_abs_torsion_angle(d12, d13, d14l, d23, d24, d34)
            high = distances_2_abs_torsion_angle(d12, d13, d14u, d23, d24, d34)
            T[row_index, 4:7] = (
                0.0,
                (high + low) * 90.0 / np.pi,
                (high - low) * 90.0 / np.pi,
            )
        elif clique_type == 3:
            if atom_1 not in names:
                raise ValueError(f"Atom index {atom_1} not found in df_X.")
            residue_id = residues[atom_1]
            if residue_id not in angle_rows:
                raise ValueError(f"Residue {residue_id} not found in df_A.")
            source = angle_rows[residue_id]
            columns = (4, 5) if names[atom_1] == "C" else (6, 7)
            if names[atom_1] != "C":
                source -= 1
            value = angles[source, columns[0]]
            delta = angles[source, columns[1]]
            T[row_index, 4:7] = (np.sign(value), np.abs(value), delta)
        else:
            raise ValueError(
                f"Unsupported clique type {clique_type} at row {row_index}. "
                "Expected 1, 2 or 3."
            )
    return T


def build_atom_lookup(df_X: pd.DataFrame) -> AtomLookup:
    """Map ``(residue_id, atom_name)`` pairs to their atom IDs."""
    return {
        (int(row[2]), str(row[1])): int(row[0])
        for row in df_X.itertuples(index=False)
    }


def build_atom_clique_index_matrix(
    atom_cliques: list[CliqueMatrixEntry],
    df_X: pd.DataFrame,
    df_A: pd.DataFrame,
    df_I: pd.DataFrame,
) -> np.ndarray:
    """Build and populate the remapping clique matrix."""
    lookup = build_atom_lookup(df_X)
    matrix = np.zeros((len(atom_cliques), 7), dtype=np.float64)
    for row_index, (pairs, clique_type) in enumerate(atom_cliques):
        if len(pairs) != 4:
            raise ValueError(
                "Each clique must contain exactly 4 positions, "
                f"got {len(pairs)} at row {row_index}."
            )
        for column, pair in enumerate(reversed(pairs)):
            if pair is None:
                continue
            key = (int(pair[0]), str(pair[1]))
            if key not in lookup:
                raise KeyError(f"Atom pair {key} was not found in df_X.")
            matrix[row_index, column] = lookup[key]
        matrix[row_index, 4] = float(clique_type)
    return update_torsion_matrix(matrix, df_X, df_A, df_I)


def build_backbone_torsion_matrix(
    df_X: pd.DataFrame, df_A: pd.DataFrame
) -> np.ndarray:
    """Build the simplified remap matrix for backbone-only output."""
    records = [
        (int(row[0]), str(row[1]).strip().upper(), int(row[2]))
        for row in df_X.itertuples(index=False)
    ]
    matrix = np.zeros((len(records), 7), dtype=np.float64)
    angle_rows = {int(row[0]): row for row in df_A.itertuples(index=False)}
    fields = {
        "N": (6, 7, -1, "psi"),
        "CA": (2, 3, 0, "omega"),
        "C": (4, 5, 0, "phi"),
    }
    for index, (atom_id, atom_name, residue_id) in enumerate(records):
        matrix[index, 0] = atom_id
        for column in range(1, min(index, 3) + 1):
            matrix[index, column] = records[index - column][0]
        if index < 3:
            continue
        if atom_name not in fields:
            raise ValueError(
                f"Unsupported atom {atom_name!r} in backbone remap structure. "
                "Expected only N, CA, or C."
            )
        value_column, delta_column, offset, label = fields[atom_name]
        source = residue_id + offset
        if source not in angle_rows:
            raise ValueError(
                f"{label} for residue {source} was not found in the parser angular file."
            )
        value = float(angle_rows[source][value_column])
        matrix[index, 4:7] = (
            np.sign(value), np.abs(value), float(angle_rows[source][delta_column])
        )
    return matrix

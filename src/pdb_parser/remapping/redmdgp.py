"""iDMDGP re-orders (``order_type: redmdgp``) for ``backbone_plus_hydrogens`` outputs.

Each DDGP ``order_id`` (1 to 10) is rewritten as a re-order in the sense of Lavor et al.,
*Minimal NMR distance information for rigidity of protein graphs*, Discrete Appl. Math. 256
(2019) 91-104, Def. 3: every vertex ``i >= 4`` is placed from its three immediate predecessors
``i-1, i-2, i-3``. Atoms may be repeated; a repeated atom becomes a new vertex with the same
coordinates (a copy).

Each re-order has one pattern for residue 1, one for residues 2 to n-1, and one for residue n:

- residues 2 to n-1 use ``REDMDGP_INTERNAL_PATTERNS[order_id]``; offsets are relative to the
  current residue (``-1`` = previous residue). The new atoms appear in the order of the DDGP
  ``order_id`` and each one is placed with the same kind of information as its DDGP clique
  (exact distance, NOE, or phi/psi interval);
- residue n uses the same pattern cut after its last new atom (the remaining copies only prepare
  the next residue);
- residue 1 uses ``REDMDGP_FIRST_RESIDUE_PATTERNS[(atom_count, entry)]``: ``atom_count`` is the
  size of the DDGP first-residue layout (7, 6 or 5 atoms, see ``first_residue_order``) and the
  entry is the window where the internal pattern starts (``A`` = N CA C for orders 1 to 6,
  ``B`` = ha C CA for orders 7 to 10).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pdb_parser.geometry import distances_2_abs_torsion_angle, torsion_angle_with_points

from .contracts import AtomRefList
from .ordering import first_residue_order, get_ha_group, get_nonterminal_h_group

PatternEntry = tuple[int, str]
DistanceBounds = tuple[float, float]
DistanceLookup = dict[frozenset[int], DistanceBounds]

# Residues 2 to n-1. "h" = H (HD2/HD3 for PRO), "ha" = HA (HA2/HA3 for GLY).
REDMDGP_INTERNAL_PATTERNS: dict[int, tuple[PatternEntry, ...]] = {
    # N CA C N CA C-1 h CA-1 C-1 N CA C ha N CA C
    1: ((0, "N"), (0, "CA"), (0, "C"), (0, "N"), (0, "CA"), (-1, "C"), (0, "h"), (-1, "CA"), (-1, "C"), (0, "N"), (0, "CA"), (0, "C"), (0, "ha"), (0, "N"), (0, "CA"), (0, "C")),
    # N CA C ha CA C N CA C-1 h CA-1 C-1 N CA C
    2: ((0, "N"), (0, "CA"), (0, "C"), (0, "ha"), (0, "CA"), (0, "C"), (0, "N"), (0, "CA"), (-1, "C"), (0, "h"), (-1, "CA"), (-1, "C"), (0, "N"), (0, "CA"), (0, "C")),
    # N CA h N CA ha N CA C
    3: ((0, "N"), (0, "CA"), (0, "h"), (0, "N"), (0, "CA"), (0, "ha"), (0, "N"), (0, "CA"), (0, "C")),
    # N CA h C-1 N CA C ha N CA C
    4: ((0, "N"), (0, "CA"), (0, "h"), (-1, "C"), (0, "N"), (0, "CA"), (0, "C"), (0, "ha"), (0, "N"), (0, "CA"), (0, "C")),
    # N h CA N ha C N CA C
    5: ((0, "N"), (0, "h"), (0, "CA"), (0, "N"), (0, "ha"), (0, "C"), (0, "N"), (0, "CA"), (0, "C")),
    # N h C-1 N CA C ha N CA C
    6: ((0, "N"), (0, "h"), (-1, "C"), (0, "N"), (0, "CA"), (0, "C"), (0, "ha"), (0, "N"), (0, "CA"), (0, "C")),
    # h N CA ha C CA
    7: ((0, "h"), (0, "N"), (0, "CA"), (0, "ha"), (0, "C"), (0, "CA")),
    # h C-1 N CA C N ha C CA
    8: ((0, "h"), (-1, "C"), (0, "N"), (0, "CA"), (0, "C"), (0, "N"), (0, "ha"), (0, "C"), (0, "CA")),
    # h CA N ha C CA
    9: ((0, "h"), (0, "CA"), (0, "N"), (0, "ha"), (0, "C"), (0, "CA")),
    # h CA C-1 N CA C ha CA C
    10: ((0, "h"), (0, "CA"), (-1, "C"), (0, "N"), (0, "CA"), (0, "C"), (0, "ha"), (0, "CA"), (0, "C")),
}

# Window of residue 1 after which the internal pattern starts: A = (N, CA, C), B = (ha, C, CA).
# Order 10 ends each internal residue on (ha, CA, C) instead of B; h is placed from either window
# by the same ha(k-1)-h(k) NOE, and B keeps residue 1 one position shorter.
REDMDGP_ENTRY = {order_id: ("A" if order_id <= 6 else "B") for order_id in REDMDGP_INTERNAL_PATTERNS}

# Residue 1, by DDGP first-residue layout (``first_residue_order``):
#   7 atoms: n1 n2 n3 N CA ha C   (three N-terminal hydrogens, e.g. H3 H2 H1)
#   6 atoms: n1 n2 N CA ha C      (two N-terminal hydrogens, e.g. H2 H1)
#   5 atoms: N CA C ha h          (one N-terminal hydrogen)
REDMDGP_FIRST_RESIDUE_LAYOUTS: dict[int, tuple[str, ...]] = {
    7: ("n1", "n2", "n3", "N", "CA", "ha", "C"),
    6: ("n1", "n2", "N", "CA", "ha", "C"),
    5: ("N", "CA", "C", "ha", "h"),
}
REDMDGP_FIRST_RESIDUE_PATTERNS: dict[tuple[int, str], tuple[str, ...]] = {
    (7, "A"): ("n1", "n2", "n3", "N", "CA", "ha", "N", "CA", "C"),
    (7, "B"): ("n1", "n2", "n3", "N", "CA", "ha", "C", "CA"),
    (6, "A"): ("n1", "n2", "N", "CA", "ha", "N", "CA", "C"),
    (6, "B"): ("n1", "n2", "N", "CA", "ha", "C", "CA"),
    # One N-terminal hydrogen: its only non-exact information is the h-ha NOE, so leaving the
    # window after h needs one copy of ha placed by that NOE.
    (5, "A"): ("N", "CA", "C", "ha", "N", "CA", "h", "N", "CA", "ha", "N", "CA", "C"),
    (5, "B"): ("N", "CA", "C", "ha", "N", "CA", "h", "N", "CA", "ha", "C", "CA"),
}


def last_residue_pattern(order_id: int) -> tuple[PatternEntry, ...]:
    """Return the internal pattern cut after the last new atom of the residue."""
    pattern = REDMDGP_INTERNAL_PATTERNS[order_id]
    seen: set[PatternEntry] = set()
    last_new = 0
    for index, entry in enumerate(pattern):
        if entry[0] == 0 and entry not in seen:
            last_new = index
        seen.add(entry)
    return pattern[: last_new + 1]


def _residue_atom_ids(df_X_i: pd.DataFrame) -> dict[str, int]:
    return {
        str(atom_name).strip().upper(): int(atom_id)
        for atom_id, atom_name in zip(df_X_i[0], df_X_i[1], strict=False)
    }


def _resolve_internal_atom(df_X_i: pd.DataFrame, role: str) -> int:
    """Resolve N, CA, C, h or ha in one residue, with the DDGP remap aliases."""
    names = _residue_atom_ids(df_X_i)
    resname = str(df_X_i.iat[0, 3]).strip().upper()
    residue_id = int(df_X_i.iat[0, 2])
    candidates = {"h": get_nonterminal_h_group(resname), "ha": get_ha_group(resname)}.get(role, (role,))
    for name in candidates:
        if name in names:
            return names[name]
    raise ValueError(
        f"Re-order role {role!r} not found in residue {resname}_{residue_id}; expected one of {candidates}."
    )


def build_redmdgp_atom_sequence(
    df_X: pd.DataFrame,
    res_ids: list[int],
    order_id: int,
) -> list[int]:
    """Return the parser atom id placed at each re-order position (copies included)."""
    if order_id not in REDMDGP_INTERNAL_PATTERNS:
        raise ValueError(f"Invalid order_id for redmdgp: {order_id}. Expected an integer from 1 to 10.")
    residues = {residue_id: df_X[df_X[2] == residue_id] for residue_id in res_ids}

    first: AtomRefList = first_residue_order(residues[res_ids[0]], res_ids[0])
    layout = REDMDGP_FIRST_RESIDUE_LAYOUTS.get(len(first))
    if layout is None:
        raise ValueError(f"Unsupported first-residue layout with {len(first)} atoms for redmdgp.")
    first_ids = dict(zip(layout, (atom_id for atom_id, _ in first), strict=True))
    sequence = [first_ids[role] for role in REDMDGP_FIRST_RESIDUE_PATTERNS[(len(first), REDMDGP_ENTRY[order_id])]]

    for position in range(1, len(res_ids)):
        pattern = (
            last_residue_pattern(order_id)
            if position == len(res_ids) - 1
            else REDMDGP_INTERNAL_PATTERNS[order_id]
        )
        for offset, role in pattern:
            residue_id = res_ids[position + offset]
            if position + offset == 0:
                # Previous residue is residue 1: only its backbone heavy atoms are referenced.
                sequence.append(first_ids[role])
            else:
                sequence.append(_resolve_internal_atom(residues[residue_id], role))

    missing = sorted(set(df_X[0].astype(int)) - set(sequence))
    if missing:
        raise ValueError(f"Atoms absent from the redmdgp re-order: {missing}.")
    return sequence


def build_distance_lookup(df_I: pd.DataFrame) -> DistanceLookup:
    """Return parser distance bounds indexed by unordered atom pair."""
    pairs = zip(df_I.iloc[:, 0].astype(int).tolist(), df_I.iloc[:, 1].astype(int).tolist())
    bounds = zip(df_I.iloc[:, 4].astype(float).tolist(), df_I.iloc[:, 5].astype(float).tolist())
    return {frozenset(pair): bound for pair, bound in zip(pairs, bounds)}


def _is_exact(bounds: DistanceBounds) -> bool:
    return bounds[0] == bounds[1]


def validate_redmdgp_sequence(
    sequence: list[int],
    df_X: pd.DataFrame,
    distances: DistanceLookup,
) -> None:
    """Check Def. 3 of a re-order on the parser data; raise ``ValueError`` on the first violation."""
    names = {int(row[0]): f"{row[1]}_{int(row[2])}" for row in df_X.itertuples(index=False)}
    coordinates = {int(row[0]): np.array(row[4:7], dtype=float) for row in df_X.itertuples(index=False)}

    def fail(position: int, reason: str) -> None:
        window = " ".join(names[atom] for atom in sequence[max(0, position - 3): position + 1])
        raise ValueError(f"redmdgp re-order is invalid at position {position + 1} ({window}): {reason}.")

    for position, atom in enumerate(sequence):
        previous = sequence[max(0, position - 3): position]
        if position < 3:
            for other in previous:
                bounds = distances.get(frozenset((atom, other)))
                if other == atom or bounds is None or not _is_exact(bounds) or bounds[0] <= 0.0:
                    fail(position, "the initial clique needs exact, positive distances")
            continue
        p3, p2, p1 = previous
        for other, label in ((p1, "i-1"), (p2, "i-2")):
            bounds = distances.get(frozenset((atom, other)))
            if other == atom or bounds is None or not _is_exact(bounds) or bounds[0] <= 0.0:
                fail(position, f"the distance to {label} must be exact and positive")
        if p3 != atom and frozenset((atom, p3)) not in distances:
            fail(position, "no distance to i-3 in the parser I file")
        x3, x2, x1 = (coordinates[other] for other in (p3, p2, p1))
        a, b, c = np.linalg.norm(x3 - x2), np.linalg.norm(x2 - x1), np.linalg.norm(x3 - x1)
        if not (a + b > c and a + c > b and b + c > a):
            fail(position, "the reference triangle is degenerate")


def build_redmdgp_distance_rows(
    sequence: list[int],
    df_I: pd.DataFrame,
    distances: DistanceLookup,
) -> list[tuple[int, int, float, float]]:
    """Return ``(i, j, d_l, d_u)`` rows (1-based positions, ``i > j``) for the re-order I file.

    - every position gets its pairs with i-1, i-2 and i-3 (``[0, 0]`` when i-3 is the same atom);
    - a new atom also gets every parser restraint with atoms already placed, on their first
      occurrence (unless that atom is in the window);
    - a copy whose i-3 distance is an interval also gets ``[0, 0]`` with its first occurrence, so
      the solver intersects that interval down to the angle already chosen.
    Each parser restraint appears on exactly one pair of first occurrences or window positions.
    """
    partners: dict[int, list[int]] = {}
    for a, b in zip(df_I.iloc[:, 0].astype(int).tolist(), df_I.iloc[:, 1].astype(int).tolist()):
        partners.setdefault(a, []).append(b)
        partners.setdefault(b, []).append(a)

    first_position: dict[int, int] = {}
    rows: list[tuple[int, int, float, float]] = []
    for position, atom in enumerate(sequence):
        window = range(max(0, position - 3), position)
        for other_position in window:
            other = sequence[other_position]
            bounds = (0.0, 0.0) if other == atom else distances[frozenset((atom, other))]
            rows.append((position + 1, other_position + 1, *bounds))
        if atom not in first_position:
            window_atoms = {sequence[other_position] for other_position in window}
            for other in partners.get(atom, ()):
                if other in first_position and other not in window_atoms:
                    rows.append((position + 1, first_position[other] + 1, *distances[frozenset((atom, other))]))
            first_position[atom] = position
        elif position >= 3 and sequence[position - 3] != atom:
            if not _is_exact(distances[frozenset((atom, sequence[position - 3]))]):
                rows.append((position + 1, first_position[atom] + 1, 0.0, 0.0))
    return rows


def build_redmdgp_distance_dataframe(
    sequence: list[int],
    df_X: pd.DataFrame,
    df_I: pd.DataFrame,
    distances: DistanceLookup,
) -> pd.DataFrame:
    """Return the re-order I table in the positional I layout (sorted by ``i, j``)."""
    atoms = {int(row[0]): (str(row[1]), int(row[2]), str(row[3])) for row in df_X.itertuples(index=False)}
    records = []
    for i, j, d_l, d_u in build_redmdgp_distance_rows(sequence, df_I, distances):
        name_i, resid_i, resname_i = atoms[sequence[i - 1]]
        name_j, resid_j, resname_j = atoms[sequence[j - 1]]
        records.append([i, j, resid_i, resid_j, d_l, d_u, name_i, name_j, resname_i, resname_j])
    return pd.DataFrame(records).sort_values(by=[0, 1]).reset_index(drop=True)


def _torsion_from_angles(
    window_atoms: tuple[int, int, int, int],
    atoms: dict[int, tuple[str, int]],
    angles: dict[int, tuple],
) -> tuple[float, float, float] | None:
    """Return ``(sign, |tau|, delta)`` from A when the window is the phi or psi dihedral.

    phi(k) is C(k-1)-N(k)-CA(k)-C(k) and psi(k-1) is N(k-1)-CA(k-1)-C(k-1)-N(k); a dihedral read
    backwards has the same value.
    """
    names = tuple(atoms[atom] for atom in window_atoms)
    for candidate in (names, names[::-1]):
        (n3, r3), (n2, r2), (n1, r1), (n0, r0) = candidate
        if (n3, n2, n1, n0) == ("C", "N", "CA", "C") and r3 == r2 - 1 and r2 == r1 == r0:
            value, delta = float(angles[r0][4]), float(angles[r0][5])
            return float(np.sign(value)), abs(value), delta
        if (n3, n2, n1, n0) == ("N", "CA", "C", "N") and r3 == r2 == r1 and r0 == r1 + 1:
            value, delta = float(angles[r1][6]), float(angles[r1][7])
            return float(np.sign(value)), abs(value), delta
    return None


def build_redmdgp_clique_matrix(
    sequence: list[int],
    df_X: pd.DataFrame,
    df_A: pd.DataFrame,
    distances: DistanceLookup,
) -> np.ndarray:
    """Return T rows ``(i, i-1, i-2, i-3, sign, |tau|, delta)`` for the re-order (degrees).

    - i-3 is the same atom: torsion 0;
    - exact i-3 distance: signed torsion of the reference coordinates;
    - phi/psi dihedral window: the A interval, with its sign;
    - any other interval (NOE, or a torsion-derived distance in another window): the |tau| interval
      equivalent to the i-3 distance interval, with unknown sign (0).
    """
    coordinates = {int(row[0]): np.array(row[4:7], dtype=float) for row in df_X.itertuples(index=False)}
    atoms = {int(row[0]): (str(row[1]), int(row[2])) for row in df_X.itertuples(index=False)}
    # A columns: resid, resname, omega, delta_omega, phi, delta_phi, psi, delta_psi.
    angles = {int(row[0]): tuple(row) for row in df_A.itertuples(index=False)}

    matrix = np.zeros((len(sequence), 7), dtype=float)
    for position, atom in enumerate(sequence):
        vertex = position + 1
        matrix[position, 0] = vertex
        for column in range(1, min(4, vertex)):
            matrix[position, column] = vertex - column
        if position < 3:
            continue
        p3, p2, p1 = sequence[position - 3: position]
        if p3 == atom:
            continue
        bounds = distances[frozenset((atom, p3))]
        if _is_exact(bounds):
            tau = torsion_angle_with_points(coordinates[p3], coordinates[p2], coordinates[p1], coordinates[atom])
            matrix[position, 4:7] = (np.sign(tau), abs(tau) * 180.0 / np.pi, 0.0)
            continue
        from_angles = _torsion_from_angles((p3, p2, p1, atom), atoms, angles)
        if from_angles is not None:
            matrix[position, 4:7] = from_angles
            continue
        x1, x2, x3, x4 = (coordinates[other] for other in (p3, p2, p1, atom))
        d12, d13, d23 = np.linalg.norm(x1 - x2), np.linalg.norm(x1 - x3), np.linalg.norm(x2 - x3)
        d24, d34 = np.linalg.norm(x2 - x4), np.linalg.norm(x3 - x4)
        low = distances_2_abs_torsion_angle(d12, d13, bounds[0], d23, d24, d34)
        high = distances_2_abs_torsion_angle(d12, d13, bounds[1], d23, d24, d34)
        low, high = min(low, high), max(low, high)
        matrix[position, 4:7] = (0.0, (high + low) * 90.0 / np.pi, (high - low) * 90.0 / np.pi)
    return matrix

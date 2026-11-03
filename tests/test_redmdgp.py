"""iDMDGP re-orders (``order_type: redmdgp``)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pdb_parser.api import run_remap_job
from pdb_parser.geometry import torsion_angle_with_points
from pdb_parser.io.files import read_space_separated_file
from pdb_parser.output_layout import (
    parser_distance_file,
    remapped_clique_file,
    remapped_distance_file,
    remapped_file_suffix,
    remapped_structure_file,
)
from pdb_parser.pipeline.parser_pipeline import parser
from pdb_parser.public_config import parse_public_order_type
from pdb_parser.remapping.ordering import get_internal_residue_ordering_template
from pdb_parser.remapping.redmdgp import (
    REDMDGP_ENTRY,
    REDMDGP_FIRST_RESIDUE_LAYOUTS,
    REDMDGP_FIRST_RESIDUE_PATTERNS,
    REDMDGP_INTERNAL_PATTERNS,
    last_residue_pattern,
)

FIXTURE = Path(__file__).parent / "fixtures" / "1ng7_model1_residues_1_6.pdb"
ENTRY_WINDOWS = {"A": ("N", "CA", "C"), "B": ("ha", "C", "CA")}
POSITIONS = {1: 16, 2: 15, 3: 9, 4: 11, 5: 9, 6: 10, 7: 6, 8: 9, 9: 6, 10: 9}
LAST_RESIDUE_CUT = {1: 3, 2: 5, 3: 0, 4: 3, 5: 3, 6: 3, 7: 1, 8: 2, 9: 1, 10: 2}
DDGP_ROLE = {("N",): "N", ("CA",): "CA", ("C",): "C", ("H",): "h", ("HA",): "ha"}


@pytest.mark.parametrize("order_id", range(1, 11))
def test_internal_pattern_places_the_new_atoms_in_ddgp_order(order_id: int) -> None:
    pattern = REDMDGP_INTERNAL_PATTERNS[order_id]
    new_atoms = list(dict.fromkeys(role for offset, role in pattern if offset == 0))

    assert new_atoms == [DDGP_ROLE[group] for group in get_internal_residue_ordering_template(order_id, "ALA")]
    assert len(pattern) == POSITIONS[order_id]
    # A copy is never one of the two previous positions (repetitions at least 3 apart).
    assert all(pattern[i] not in pattern[i - 2:i] for i in range(2, len(pattern)))


@pytest.mark.parametrize("order_id", range(1, 11))
def test_last_residue_pattern_is_cut_after_its_last_new_atom(order_id: int) -> None:
    pattern = REDMDGP_INTERNAL_PATTERNS[order_id]
    last = last_residue_pattern(order_id)

    assert pattern[: len(last)] == last
    assert len(pattern) - len(last) == LAST_RESIDUE_CUT[order_id]
    assert {entry for entry in last if entry[0] == 0} == {entry for entry in pattern if entry[0] == 0}


@pytest.mark.parametrize(("atom_count", "entry"), sorted(REDMDGP_FIRST_RESIDUE_PATTERNS))
def test_first_residue_pattern_covers_its_layout_and_ends_on_its_entry(atom_count: int, entry: str) -> None:
    pattern = REDMDGP_FIRST_RESIDUE_PATTERNS[(atom_count, entry)]

    assert set(pattern) == set(REDMDGP_FIRST_RESIDUE_LAYOUTS[atom_count])
    assert pattern[-3:] == ENTRY_WINDOWS[entry]
    assert all(pattern[i] not in pattern[i - 2:i] for i in range(2, len(pattern)))


def test_entries_follow_where_each_order_starts() -> None:
    assert {order_id for order_id, entry in REDMDGP_ENTRY.items() if entry == "A"} == set(range(1, 7))
    assert {order_id for order_id, entry in REDMDGP_ENTRY.items() if entry == "B"} == set(range(7, 11))


def test_order_type_is_optional_and_validated() -> None:
    assert parse_public_order_type({"order_id": "3"}) == "ddgp"
    assert parse_public_order_type({"order_id": "3", "order_type": " REDMDGP "}) == "redmdgp"
    with pytest.raises(ValueError, match="Invalid order_type"):
        parse_public_order_type({"order_id": "3", "order_type": "ddmgp"})


def test_redmdgp_outputs_have_their_own_suffix() -> None:
    assert remapped_file_suffix("1abc", 1, "A", 7, "backbone_plus_hydrogens", "redmdgp") == (
        "1ABC_model1_chainA_redmdgp_hc_order_7.dat"
    )
    # Backbone outputs keep the fixed DMDGP order: order_type is ignored like order_id.
    assert remapped_file_suffix("1abc", 1, "A", 7, "backbone", "redmdgp") == "1ABC_model1_chainA_dmdgp_order.dat"


def _write_fixture(pdb_dir: Path, n_terminal_hydrogens: tuple[str, ...]) -> None:
    """1NG7 residues 1-6 (MET GLY PRO LEU GLN TYR) keeping only the given N-terminal hydrogens."""
    removed = {"H1", "H2", "H3"} - set(n_terminal_hydrogens)
    lines = [
        line for line in FIXTURE.read_text(encoding="utf-8").splitlines()
        if not (line.startswith("ATOM") and line[22:26].strip() == "1" and line[12:16].strip() in removed)
    ]
    pdb_dir.mkdir(parents=True, exist_ok=True)
    (pdb_dir / "TEST.pdb").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _parse(tmp_path: Path, atom_selection: str = "backbone_plus_hydrogens") -> Path:
    out_dir = tmp_path / "outputs"
    params = {
        "method": "nmr",
        "model_number": "1",
        "chain_id": "A",
        "atom_selection": atom_selection,
        "distance_constraints": "interval_experimental",
        "noe_strong": "2.5",
        "noe_medium": "3.5",
        "noe_weak": "5.0",
        "vdwr_hh": "1.8",
        "vdw_constraints": "yes",
        "tightening_upper_bounds": "no",
        "torsion_angle_width": "60.0",
        "percentage_backbone_torsion_angles": "100.0",
    }
    parser(params, tmp_path / "pdb", tmp_path / "seeds", out_dir, "TEST", True)
    return out_dir


@pytest.mark.parametrize(
    ("n_terminal_hydrogens", "first_residue_atoms"),
    [(("H1", "H2", "H3"), 7), (("H1", "H2"), 6), (("H1",), 5)],
)
def test_redmdgp_remap_writes_a_consistent_reorder_for_every_order(
    tmp_path: Path, n_terminal_hydrogens: tuple[str, ...], first_residue_atoms: int
) -> None:
    _write_fixture(tmp_path / "pdb", n_terminal_hydrogens)
    out_dir = _parse(tmp_path)
    parser_pairs = {
        frozenset(((int(row[2]), row[6]), (int(row[3]), row[7])))
        for row in read_space_separated_file(parser_distance_file(out_dir, "TEST", 1, "A")).itertuples(index=False, name=None)
    }

    for order_id in range(1, 11):
        params = {"model_number": "1", "chain_id": "A", "order_id": str(order_id), "order_type": "redmdgp"}
        assert run_remap_job(params, out_dir, "TEST", order_id) == 0
        files = [
            function(out_dir, "TEST", 1, "A", order_id, "backbone_plus_hydrogens", "redmdgp")
            for function in (remapped_structure_file, remapped_distance_file, remapped_clique_file)
        ]
        x_rows, i_rows, t_rows = (read_space_separated_file(path) for path in files)
        n = len(x_rows)
        atom_at = {}
        for i, j, resid_i, resid_j, _, _, name_i, name_j, *_ in i_rows.itertuples(index=False, name=None):
            atom_at[int(i)], atom_at[int(j)] = (int(resid_i), name_i), (int(resid_j), name_j)
        xyz = x_rows.to_numpy(dtype=float)

        assert sorted(atom_at) == list(range(1, n + 1))
        pairs = list(zip(i_rows[0], i_rows[1], strict=True))
        assert len(pairs) == len(set(pairs)) and (i_rows[0] > i_rows[1]).all()
        for i, j, d_l, d_u in zip(i_rows[0], i_rows[1], i_rows[4], i_rows[5], strict=True):
            distance = float(np.linalg.norm(xyz[i - 1] - xyz[j - 1]))
            if d_l == d_u:
                assert distance == pytest.approx(d_l, abs=2e-3)
            if d_u == 0.0:
                assert atom_at[i] == atom_at[j]
        # Every parser restraint is kept, on some pair of positions of the same atoms.
        assert parser_pairs <= {frozenset((atom_at[i], atom_at[j])) for i, j in pairs}
        # T is a DMDGP order: every vertex is placed from i-1, i-2, i-3.
        assert (t_rows[0] == np.arange(1, n + 1)).all()
        for column in (1, 2, 3):
            assert (t_rows[column].to_numpy()[3:] == np.arange(4, n + 1) - column).all()
        first_residue = [atom for position, atom in sorted(atom_at.items()) if atom[0] == 1]
        assert len(set(first_residue)) == first_residue_atoms
        # Every copy has four predecessors with exactly one [0, 0], with the previous occurrence
        # of its atom; they are all exact except i-3 when it is an interval (the zero is then
        # the extra row). A torsion of 0 and a phi/psi window (also with its end atoms swapped)
        # have a known sign, and a known sign contains the reference torsion.
        rows_of = {}
        for i, j, d_l, d_u in zip(i_rows[0], i_rows[1], i_rows[4], i_rows[5], strict=True):
            rows_of.setdefault(i, {})[j] = (d_l, d_u)
        previous = {}
        for k, (sign, abs_tau, delta) in enumerate(t_rows[[4, 5, 6]].itertuples(index=False, name=None), 1):
            zeros = [j for j, (d_l, d_u) in rows_of.get(k, {}).items() if d_u == 0.0]
            assert zeros == ([previous[atom_at[k]]] if atom_at[k] in previous else []), (order_id, k)
            if atom_at[k] in previous:
                rows = rows_of[k]
                assert len(rows) == 4, (order_id, k, rows)
                assert all(d_l == d_u for j, (d_l, d_u) in rows.items() if j != k - 3 or k - 3 in zeros), (order_id, k)
            previous[atom_at[k]] = k
            if k < 4:
                continue
            assert sign != 0 or delta > 0, (order_id, k)
            if sign != 0:
                native = np.degrees(torsion_angle_with_points(*(xyz[v - 1] for v in (k - 3, k - 2, k - 1, k))))
                gap = abs((native - sign * abs_tau + 180.0) % 360.0 - 180.0)
                assert gap <= delta + 1e-6, (order_id, k, native, sign, abs_tau, delta)


def test_redmdgp_is_ignored_for_backbone_outputs(tmp_path: Path) -> None:
    _write_fixture(tmp_path / "pdb", ("H1", "H2", "H3"))
    out_dir = _parse(tmp_path, atom_selection="backbone")

    written = {}
    for order_type in ("ddgp", "redmdgp"):
        params = {"model_number": "1", "chain_id": "A", "order_id": "1", "order_type": order_type}
        assert run_remap_job(params, out_dir, "TEST", 1) == 0
        written[order_type] = remapped_clique_file(out_dir, "TEST", 1, "A", 1, "backbone").read_text(encoding="utf-8")

    assert written["redmdgp"] == written["ddgp"]
    assert not list((out_dir / "TEST" / "remapped").glob("*redmdgp*"))

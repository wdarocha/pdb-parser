"""Files written by earlier releases must stay readable.

Readers split on whitespace, so the column widths and decimal places of a file do not matter.
These tests pin that contract with the layouts written before 0.3.0: parser-stage ``X_*`` with
``%8.3f%8.3f%8.3f`` (rows whose coordinates are not glued together) and ``I_*`` with ``%20.16f``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pdb_parser.api import run_remap_job
from pdb_parser.io.files import (
	DISTANCE_ROW_FORMAT,
	read_validated_positional_file,
	save_distances_from_df_structure,
)
from pdb_parser.output_layout import remapped_clique_file, remapped_distance_file

LEGACY_X_ROW = "%5d %-4s %4d %3s %8.3f%8.3f%8.3f\n"
LEGACY_I_ROW = "%5d %5d %6d %6d %20.16f %20.16f %4.4s %4.4s %s %s\n"
LEGACY_A_ROW = "%4d %3s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n"

ATOMS = [
	(1, "N", 1, "ALA", 0.0, 0.0, 0.0),
	(2, "CA", 1, "ALA", 1.0, 0.0, 0.0),
	(3, "C", 1, "ALA", 1.0, 1.0, 0.0),
	(4, "HA", 1, "ALA", 0.0, 1.0, 0.0),
	(5, "H", 1, "ALA", 0.0, 0.0, 1.0),
	(6, "N", 2, "GLY", 2.0, 1.0, 0.0),
	(7, "CA", 2, "GLY", 2.0, 2.0, 0.0),
	(8, "C", 2, "GLY", 3.0, 2.0, 0.0),
	(9, "HA2", 2, "GLY", 2.0, 2.0, 1.0),
	(10, "H", 2, "GLY", 2.0, 1.0, 1.0),
]
DISTANCES = [
	(5, 4, 1, 1, 1.7320508075688772, 1.7320508075688772, "H", "HA", "ALA", "ALA"),
	(10, 4, 2, 1, 2.4142135623730949, 3.1415926535897931, "H", "HA", "GLY", "ALA"),
	(10, 5, 2, 1, 2.2360679774997898, 2.2360679774997898, "H", "H", "GLY", "ALA"),
	(9, 4, 2, 1, 1.8, 999.0, "HA2", "HA", "GLY", "ALA"),
	(9, 5, 2, 1, 2.4, 2.4, "HA2", "H", "GLY", "ALA"),
	(9, 10, 2, 2, 1.4142135623730951, 1.4142135623730951, "HA2", "H", "GLY", "GLY"),
]


def _write_legacy_parser_outputs(pdb_dir: Path) -> None:
	pdb_dir.mkdir(parents=True)
	(pdb_dir / "X_1ABC_model1_chainA.dat").write_text("".join(LEGACY_X_ROW % atom for atom in ATOMS))
	(pdb_dir / "I_1ABC_model1_chainA.dat").write_text("".join(LEGACY_I_ROW % row for row in DISTANCES))
	(pdb_dir / "A_1ABC_model1_chainA.dat").write_text(
		LEGACY_A_ROW % (1, "ALA", 0, 0, 0, 10, 0, 10) + LEGACY_A_ROW % (2, "GLY", 0, 0, 0, 10, 0, 10)
	)


def test_legacy_parser_outputs_are_read_with_their_full_precision(tmp_path: Path) -> None:
	pdb_dir = tmp_path / "1ABC"
	_write_legacy_parser_outputs(pdb_dir)

	df_x = read_validated_positional_file(pdb_dir / "X_1ABC_model1_chainA.dat", "X")
	df_i = read_validated_positional_file(pdb_dir / "I_1ABC_model1_chainA.dat", "I")
	read_validated_positional_file(pdb_dir / "A_1ABC_model1_chainA.dat", "A")

	assert df_x.iloc[:, 4:7].to_numpy(dtype=float).tolist() == [list(atom[4:]) for atom in ATOMS]
	# %20.16f keeps every digit a double carries, so the old values come back exactly.
	assert df_i.iloc[:, 4:6].to_numpy(dtype=float).tolist() == [list(row[4:6]) for row in DISTANCES]


def test_legacy_parser_outputs_can_still_be_remapped(tmp_path: Path) -> None:
	out_dir = tmp_path / "outputs"
	_write_legacy_parser_outputs(out_dir / "1ABC")

	assert run_remap_job({"model_number": "1", "chain_id": "A", "order_id": "9"}, out_dir, "1ABC", 9) == 0

	remapped_i = read_validated_positional_file(remapped_distance_file(out_dir, "1ABC", 1, "A", 9), "I")
	assert sorted(remapped_i.iloc[:, 5].round(12)) == sorted(round(row[5], 12) for row in DISTANCES)
	assert remapped_clique_file(out_dir, "1ABC", 1, "A", 9).read_text().strip()


def test_twelve_decimal_rewrite_keeps_values_and_interval_order(tmp_path: Path) -> None:
	"""Rewriting an I file at 12 decimals moves each value by at most 5e-13, keeps d_l <= d_u,
	keeps exact (d_l == d_u) constraints exact, and is idempotent."""
	rng = np.random.default_rng(0)
	lower = rng.uniform(0.0, 999.0, 2000)
	upper = lower + rng.choice([0.0, 1e-13, 1e-9, 1.0], lower.size)
	rows = [(i, i + 1, 1, 1, lo, up, "N", "CA", "ALA", "ALA") for i, (lo, up) in enumerate(zip(lower, upper), 1)]
	first = tmp_path / "I_first.dat"
	first.write_text("".join(DISTANCE_ROW_FORMAT % row for row in rows))

	df = read_validated_positional_file(first, "I")
	d_l, d_u = df.iloc[:, 4].to_numpy(dtype=float), df.iloc[:, 5].to_numpy(dtype=float)
	# half of the last printed decimal, plus the spacing of doubles near the value
	assert (np.abs(d_l - lower) <= 5e-13 + np.spacing(lower)).all()
	assert (np.abs(d_u - upper) <= 5e-13 + np.spacing(upper)).all()
	assert (d_l <= d_u).all()
	assert (d_l[lower == upper] == d_u[lower == upper]).all()

	second = tmp_path / "I_second.dat"
	save_distances_from_df_structure(df, second)
	assert second.read_text() == first.read_text()

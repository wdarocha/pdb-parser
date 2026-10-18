from __future__ import annotations

from pathlib import Path

import numpy as np

from pdb_parser.io.distance_constraints import nmr_distance_constraints


def write_hydrogen_tsv(path: Path) -> Path:
	path.write_text(
		(
			"atom_id\tatom_name\tresid\tresname\tx\ty\tz\n"
			"1\tH\t1\tALA\t0.0\t0.0\t0.0\n"
			"2\tH\t2\tALA\t2.0\t0.0\t0.0\n"
			"3\tH\t3\tALA\t4.0\t0.0\t0.0\n"
		),
		encoding="utf-8",
	)
	return path


def test_nmr_distance_constraints_precise_writes_exact_intervals(tmp_path: Path) -> None:
	input_file = write_hydrogen_tsv(tmp_path / "filtered_atoms.dat")
	output_file = tmp_path / "distance_constraints_4.dat"

	nmr_distance_constraints(
		tsv_structure_file=input_file,
		output_distance_file=output_file,
		distance_constraints="precise",
		atom_selection="backbone_plus_hydrogens",
		epsilon_short=0.0,
		epsilon_long=0.0,
		max_distance=3.0,
		noe_strong=0.0,
		noe_medium=0.0,
		noe_weak=0.0,
		vdw_threshold=0.0,
		rng=None,
	)

	lines = output_file.read_text(encoding="utf-8").splitlines()

	assert lines[0] == "atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j"
	assert "2\t1\t2\t1\t2.000000000000\t2.000000000000\tH\tH\tALA\tALA" in lines
	assert "3\t2\t3\t2\t2.000000000000\t2.000000000000\tH\tH\tALA\tALA" in lines
	assert all("3\t1\t3\t1" not in line for line in lines[1:])


def test_nmr_distance_constraints_experimental_assigns_expected_noe_classes(tmp_path: Path) -> None:
	input_file = write_hydrogen_tsv(tmp_path / "filtered_atoms.dat")
	output_file = tmp_path / "distance_constraints_4.dat"

	nmr_distance_constraints(
		tsv_structure_file=input_file,
		output_distance_file=output_file,
		distance_constraints="interval_experimental",
		atom_selection="backbone_plus_hydrogens",
		epsilon_short=0.0,
		epsilon_long=0.0,
		max_distance=0.0,
		noe_strong=2.5,
		noe_medium=3.5,
		noe_weak=5.0,
		vdw_threshold=1.8,
		rng=None,
	)

	lines = output_file.read_text(encoding="utf-8").splitlines()

	assert "2\t1\t2\t1\t1.800000000000\t2.500000000000\tH\tH\tALA\tALA" in lines
	assert "3\t2\t3\t2\t1.800000000000\t2.500000000000\tH\tH\tALA\tALA" in lines
	assert "3\t1\t3\t1\t1.800000000000\t5.000000000000\tH\tH\tALA\tALA" in lines


def test_nmr_distance_constraints_interval_centered_writes_strict_interval(tmp_path: Path) -> None:
	input_file = write_hydrogen_tsv(tmp_path / "filtered_atoms.dat")
	output_file = tmp_path / "distance_constraints_4.dat"

	nmr_distance_constraints(
		tsv_structure_file=input_file,
		output_distance_file=output_file,
		distance_constraints="interval_centered",
		atom_selection="backbone_plus_hydrogens",
		epsilon_short=0.5,
		epsilon_long=1.0,
		max_distance=3.0,
		noe_strong=0.0,
		noe_medium=0.0,
		noe_weak=0.0,
		vdw_threshold=1.8,
		rng=np.random.default_rng(12345),
	)

	lines = output_file.read_text(encoding="utf-8").splitlines()

	assert len(lines) == 3
	assert all("3\t1\t3\t1" not in line for line in lines[1:])

	for line in lines[1:]:
		parts = line.split("\t")
		d_l = float(parts[4])
		d_u = float(parts[5])
		assert d_l < 2.0 < d_u
		assert d_l >= 1.8
		assert d_u <= 3.0

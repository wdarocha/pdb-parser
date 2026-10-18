from __future__ import annotations

from pathlib import Path

import numpy as np

from pdb_parser.io.distance_constraints import (
	build_covalent_and_planar_topology,
	covalent_and_planar_distances,
)
from pdb_parser.io.nmr_constraints import nmr_distance_constraints

SNAPSHOT_DIR = Path(__file__).parent / "snapshots"


def assert_matches_snapshot(snapshot_name: str, content: str) -> None:
	expected = (SNAPSHOT_DIR / snapshot_name).read_text(encoding="utf-8")
	assert content == expected


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


def write_filtered_atoms_file(path: Path) -> Path:
	path.write_text(
		(
			"atom_id\tatom_name\tresid\tresname\tx\ty\tz\n"
			"1\tn\t1\tala\t0.0\t0.0\t0.0\n"
			"2\th\t1\tala\t0.0\t0.0\t1.0\n"
			"3\tca\t1\tala\t1.0\t0.0\t0.0\n"
			"4\tha\t1\tala\t1.0\t0.0\t1.0\n"
			"5\tc\t1\tala\t2.0\t0.0\t0.0\n"
			"6\tn\t2\tala\t3.0\t0.0\t0.0\n"
			"7\th\t2\tala\t3.0\t0.0\t1.0\n"
			"8\tca\t2\tala\t4.0\t0.0\t0.0\n"
			"9\tha\t2\tala\t4.0\t0.0\t1.0\n"
			"10\tc\t2\tala\t5.0\t0.0\t0.0\n"
		),
		encoding="utf-8",
	)
	return path


def test_snapshot_covalent_and_planar_distances_output(tmp_path: Path) -> None:
	input_file = write_filtered_atoms_file(tmp_path / "filtered_atoms.dat")
	topology = build_covalent_and_planar_topology(input_file)
	output_file = tmp_path / "distance_constraints_1.dat"

	covalent_and_planar_distances(input_file, topology, output_file)

	assert_matches_snapshot("covalent_and_planar_distances.dat", output_file.read_text(encoding="utf-8"))


def test_snapshot_nmr_precise_output(tmp_path: Path) -> None:
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

	assert_matches_snapshot("nmr_precise.dat", output_file.read_text(encoding="utf-8"))


def test_snapshot_nmr_experimental_output(tmp_path: Path) -> None:
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

	assert_matches_snapshot("nmr_experimental.dat", output_file.read_text(encoding="utf-8"))


def test_snapshot_nmr_interval_centered_output(tmp_path: Path) -> None:
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

	assert_matches_snapshot("nmr_interval_centered.dat", output_file.read_text(encoding="utf-8"))

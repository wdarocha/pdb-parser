from __future__ import annotations

from pathlib import Path

import pytest

from pdb_parser.io.distance_constraints import (
	merge_distance_constraint_files,
	tightening_upper_bounds,
)
from pdb_parser.io.files import read_space_separated_file


def write_distance_instance(path: Path) -> Path:
	path.write_text(
		"\n".join(
			[
				"    1     2      1      1   1.0000000000000000   3.0000000000000000    N   CA ALA ALA",
				"    2     3      1      1   1.0000000000000000   2.0000000000000000   CA    C ALA ALA",
				"    1     3      1      1   1.0000000000000000  10.0000000000000000    N    C ALA ALA",
			]
		)
		+ "\n",
		encoding="utf-8",
	)
	return path


def write_constraint_table(path: Path, rows: list[str]) -> Path:
	path.write_text(
		"\n".join(
			[
				"atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j",
				*rows,
			]
		)
		+ "\n",
		encoding="utf-8",
	)
	return path


def test_tightening_upper_bounds_overwrites_original_file(tmp_path: Path) -> None:
	instance_file = write_distance_instance(tmp_path / "I_test.dat")
	original_df = read_space_separated_file(instance_file)
	assert float(original_df.iloc[2, 5]) == 10.0

	tightening_upper_bounds(instance_file)

	output_df = read_space_separated_file(instance_file)

	assert float(output_df.iloc[2, 5]) == 5.0


def test_merge_distance_constraint_files_prefers_smallest_precise_interval(tmp_path: Path) -> None:
	files_dir = tmp_path / "tmp"
	files_dir.mkdir()

	write_constraint_table(
		files_dir / "distance_constraints_1.dat",
		["1\t2\t1\t1\t1.0\t4.0\tN\tCA\tALA\tALA"],
	)
	write_constraint_table(
		files_dir / "distance_constraints_2.dat",
		["1\t2\t1\t1\t2.5\t2.5\tN\tCA\tALA\tALA"],
	)
	write_constraint_table(
		files_dir / "distance_constraints_3.dat",
		["1\t2\t1\t1\t3.0\t3.0\tN\tCA\tALA\tALA"],
	)

	instance_file = tmp_path / "I_merged.dat"
	merge_distance_constraint_files(files_dir, instance_file)

	merged_df = read_space_separated_file(instance_file)
	assert float(merged_df.iloc[0, 4]) == 2.5
	assert float(merged_df.iloc[0, 5]) == 2.5


def test_merge_distance_constraint_files_rejects_inconsistent_intersection(tmp_path: Path) -> None:
	files_dir = tmp_path / "tmp"
	files_dir.mkdir()

	write_constraint_table(
		files_dir / "distance_constraints_1.dat",
		["1\t2\t1\t1\t4.0\t5.0\tN\tCA\tALA\tALA"],
	)
	write_constraint_table(
		files_dir / "distance_constraints_2.dat",
		["1\t2\t1\t1\t6.0\t7.0\tN\tCA\tALA\tALA"],
	)

	with pytest.raises(ValueError, match="Inconsistent merged interval"):
		merge_distance_constraint_files(files_dir, tmp_path / "I_bad.dat")


def test_merge_distance_constraint_files_accepts_string_directory(tmp_path: Path) -> None:
	files_dir = tmp_path / "constraints"
	files_dir.mkdir()
	write_constraint_table(
		files_dir / "distance_constraints_1.dat",
		["1\t2\t1\t1\t1.0\t2.0\tN\tCA\tALA\tALA"],
	)
	output = tmp_path / "I.dat"
	merge_distance_constraint_files(str(files_dir), output)
	assert output.exists()


def test_merge_tightens_before_rounding_so_bounds_are_rounded_once(tmp_path: Path) -> None:
	"""Tightening the written file would sum already rounded bounds (0.333333333333 twice gives
	0.666666666666); tightening during the merge rounds the exact path sum once."""
	third = f"{1 / 3:.17g}"
	files_dir = tmp_path / "tmp"
	files_dir.mkdir()
	write_constraint_table(
		files_dir / "distance_constraints_1.dat",
		[
			f"1\t2\t1\t1\t0.1\t{third}\tN\tCA\tALA\tALA",
			f"2\t3\t1\t1\t0.1\t{third}\tCA\tC\tALA\tALA",
			"1\t3\t1\t1\t0.1\t999.0\tN\tC\tALA\tALA",
		],
	)

	instance_file = tmp_path / "I_merged.dat"
	merge_distance_constraint_files(files_dir, instance_file, tighten_upper_bounds=True)

	assert instance_file.read_text(encoding="utf-8").splitlines()[1].split()[5] == "0.666666666667"


def test_tightening_rejects_a_pair_without_a_finite_path(tmp_path: Path) -> None:
	"""A 999 bound is replaced by its shortest path; with no finite path it is an error."""
	files_dir = tmp_path / "tmp"
	files_dir.mkdir()
	write_constraint_table(
		files_dir / "distance_constraints_1.dat",
		[
			"1\t2\t1\t1\t1.0\t1.0\tN\tCA\tALA\tALA",
			"3\t4\t2\t2\t1.0\t1.0\tN\tCA\tGLY\tGLY",
			"3\t1\t2\t1\t2.0\t999.0\tN\tN\tGLY\tALA",
		],
	)

	with pytest.raises(ValueError, match=r"No finite shortest-path upper bound found for pair \(3, 1\)"):
		merge_distance_constraint_files(files_dir, tmp_path / "I_merged.dat", tighten_upper_bounds=True)


def test_in_memory_merge_rounds_tightened_sums_once(tmp_path: Path) -> None:
	"""The parser merges its full-precision tables: a path sum is rounded once, not per edge."""
	from pdb_parser.io.constraint_files import merge_distance_constraint_tables
	from pdb_parser.io.files import distance_table

	value = 0.1234567890124  # rounds to 0.123456789012 in a 12-decimal tmp/ file
	table = distance_table([
		(1, 2, 1, 1, 0.1, value, "N", "CA", "ALA", "ALA"),
		(2, 3, 1, 1, 0.1, value, "CA", "C", "ALA", "ALA"),
		(1, 3, 1, 1, 0.1, 999.0, "N", "C", "ALA", "ALA"),
	])
	instance_file = tmp_path / "I_merged.dat"
	merge_distance_constraint_tables([table], instance_file, tighten_upper_bounds=True)

	# 2 * 0.1234567890124 = 0.2469135780248 -> 0.246913578025 (0.246913578024 from rounded edges)
	assert instance_file.read_text(encoding="utf-8").splitlines()[1].split()[5] == "0.246913578025"

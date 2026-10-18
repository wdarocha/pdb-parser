from __future__ import annotations

from pathlib import Path

import pytest

from pdb_parser.io.distance_constraints import planar_peptide_distances


def write_planar_input_file(path: Path) -> Path:
	path.write_text(
		(
			"atom_id\tatom_name\tresid\tresname\tx\ty\tz\n"
			"1\tCA\t1\tALA\t0.0\t0.0\t0.0\n"
			"2\tO\t1\tALA\t0.0\t1.0\t0.0\n"
			"3\tH\t2\tALA\t1.0\t1.0\t0.0\n"
			"4\tCA\t2\tALA\t1.0\t0.0\t0.0\n"
		),
		encoding="utf-8",
	)
	return path


@pytest.mark.parametrize("as_strings", [False, True])
def test_planar_peptide_distances_accepts_path_and_string_inputs(
	tmp_path: Path, as_strings: bool
) -> None:
	input_file = write_planar_input_file(tmp_path / "filtered_atoms.dat")
	output_file = tmp_path / "distance_constraints_3.dat"

	arguments = (str(input_file), str(output_file)) if as_strings else (input_file, output_file)
	planar_peptide_distances(*arguments, "backbone_plus_hydrogens")

	lines = output_file.read_text(encoding="utf-8").splitlines()

	assert lines[0] == "atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j"
	assert "4\t2\t2\t1\t1.414213562373\t1.414213562373\tCA\tO\tALA\tALA" in lines
	assert "3\t1\t2\t1\t1.414213562373\t1.414213562373\tH\tCA\tALA\tALA" in lines
	assert "3\t2\t2\t1\t1.000000000000\t1.000000000000\tH\tO\tALA\tALA" in lines

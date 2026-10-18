from __future__ import annotations

from pathlib import Path

from pdb_parser.io.files import read_validated_positional_file
from pdb_parser.io.filtering import (
	convert_tsv_angular_to_pdb_format,
	convert_tsv_structure_to_pdb_format,
)


def write_text_file(path: Path, content: str) -> Path:
	path.write_text(content, encoding="utf-8")
	return path


def test_convert_tsv_structure_to_pdb_format_keeps_aligned_columns(tmp_path: Path) -> None:
	input_file = write_text_file(
		tmp_path / "filtered_atoms.dat",
		(
			"atom_id\tatom_name\tresid\tresname\tx\ty\tz\n"
			"1\tN\t10\tALA\t1.2349\t2.3451\t3.4567\n"
			"2\tCA\t10\tALA\t-1.0\t0.0\t4.0\n"
		),
	)

	output_file = tmp_path / "X_test.dat"
	convert_tsv_structure_to_pdb_format(input_file, output_file)

	assert output_file.read_text(encoding="utf-8") == (
		"    1 N      10 ALA     1.235     2.345     3.457\n"
		"    2 CA     10 ALA    -1.000     0.000     4.000\n"
	)


def test_structure_layout_stays_parseable_at_any_coordinate_magnitude(tmp_path: Path) -> None:
	"""Regression: "%8.3f%8.3f%8.3f" glued a coordinate <= -100 Å to the previous one
	("35.883-108.528"), so the remap reader saw 6 fields on some rows and 7 on others."""
	values = [
		(35.883, -108.528, 74.601),     # the 2N62 case
		(-999.999, -999.999, -999.999), # most negative value of the PDB format
		(9999.999, 9999.999, 9999.999), # largest value of the PDB format
		(-12345.678, 0.0, 123456.789),  # beyond the PDB format (e.g. mmCIF input)
	]
	rows = "".join(f"{i}\tCA\t{i}\tGLY\t{x}\t{y}\t{z}\n" for i, (x, y, z) in enumerate(values, 1))
	input_file = write_text_file(tmp_path / "filtered_atoms.dat", "atom_id\tatom_name\tresid\tresname\tx\ty\tz\n" + rows)
	output_file = tmp_path / "X_test.dat"
	convert_tsv_structure_to_pdb_format(input_file, output_file)

	df = read_validated_positional_file(output_file, "X")
	assert len(df) == len(values)
	for row, (x, y, z) in zip(df.itertuples(index=False), values):
		assert (row[4], row[5], row[6]) == (round(x, 3), round(y, 3), round(z, 3))
	assert all(len(line.split()) == 7 for line in output_file.read_text().splitlines())


def test_convert_tsv_angular_to_pdb_format_preserves_fixed_width_layout(tmp_path: Path) -> None:
	input_file = write_text_file(
		tmp_path / "angular_constraints.dat",
		(
			"resid\tresname\tomega_center\tomega_radius\tphi_center\tphi_radius\tpsi_center\tpsi_radius\n"
			"10\tALA\t180.0\t5.0\t-60.25\t10.5\t135.75\t20.0\n"
		),
	)

	output_file = tmp_path / "A_test.dat"
	convert_tsv_angular_to_pdb_format(input_file, output_file)

	assert output_file.read_text(encoding="utf-8") == (
		"  10 ALA  180.000000000000    5.000000000000  -60.250000000000   10.500000000000"
		"  135.750000000000   20.000000000000\n"
	)

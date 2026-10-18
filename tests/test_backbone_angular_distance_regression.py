from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from pdb_parser.geometry.distance_geometry import torsion_angle_with_points
from pdb_parser.io.distance_constraints import backbone_angular_interval_to_distance_interval


def write_backbone_structure_file(path: Path) -> Path:
	path.write_text(
		(
			"atom_id\tatom_name\tresid\tresname\tx\ty\tz\n"
			"1\tN\t1\tALA\t0.0\t0.0\t0.0\n"
			"2\tCA\t1\tALA\t1.0\t0.0\t0.0\n"
			"3\tC\t1\tALA\t1.0\t1.0\t0.0\n"
			"4\tN\t2\tALA\t2.0\t1.0\t1.0\n"
			"5\tCA\t2\tALA\t3.0\t1.0\t1.0\n"
			"6\tC\t2\tALA\t3.0\t2.0\t1.0\n"
		),
		encoding="utf-8",
	)
	return path


def write_angular_file(path: Path) -> Path:
	n1 = np.array([0.0, 0.0, 0.0])
	ca1 = np.array([1.0, 0.0, 0.0])
	c1 = np.array([1.0, 1.0, 0.0])
	n2 = np.array([2.0, 1.0, 1.0])
	ca2 = np.array([3.0, 1.0, 1.0])
	c2 = np.array([3.0, 2.0, 1.0])

	psi1 = np.rad2deg(torsion_angle_with_points(n1, ca1, c1, n2))
	omega2 = np.rad2deg(torsion_angle_with_points(ca1, c1, n2, ca2))
	phi2 = np.rad2deg(torsion_angle_with_points(c1, n2, ca2, c2))

	path.write_text(
		(
			"resid\tresname\tomega_center\tomega_radius\tphi_center\tphi_radius\tpsi_center\tpsi_radius\n"
			f"1\tALA\t0.000000\t0.000000\t0.000000\t0.000000\t{psi1:.6f}\t0.000000\n"
			f"2\tALA\t{omega2:.6f}\t0.000000\t{phi2:.6f}\t0.000000\t0.000000\t0.000000\n"
		),
		encoding="utf-8",
	)
	return path


@pytest.mark.parametrize("as_strings", [False, True])
def test_backbone_angular_interval_to_distance_interval_accepts_path_and_string_inputs(
	tmp_path: Path, as_strings: bool
) -> None:
	structure_file = write_backbone_structure_file(tmp_path / "filtered_atoms.dat")
	angular_file = write_angular_file(tmp_path / "angular_constraints_5.dat")
	output_file = tmp_path / "distance_constraints_5.dat"

	arguments = (str(structure_file), str(angular_file), str(output_file)) if as_strings else (
		structure_file, angular_file, output_file
	)
	backbone_angular_interval_to_distance_interval(*arguments)

	lines = output_file.read_text(encoding="utf-8").splitlines()

	assert lines[0] == "atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j"
	assert "4\t1\t2\t1\t2.449489742783\t2.449489742783\tN\tN\tALA\tALA" in lines
	assert "5\t2\t2\t1\t2.449489742783\t2.449489742783\tCA\tCA\tALA\tALA" in lines

	c_c_line = next(line for line in lines[1:] if line.startswith("6\t3\t2\t1\t"))
	parts = c_c_line.split("\t")
	assert math.isclose(float(parts[4]), 2.449489742783178, rel_tol=0.0, abs_tol=1e-12)
	assert math.isclose(float(parts[5]), 2.449489742783178, rel_tol=0.0, abs_tol=1e-12)
	assert parts[6:] == ["C", "C", "ALA", "ALA"]

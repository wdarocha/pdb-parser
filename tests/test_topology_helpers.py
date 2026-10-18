from __future__ import annotations

from pathlib import Path

import pytest

from pdb_parser.io.distance_constraints import (
	build_covalent_and_planar_topology,
	covalent_and_planar_distances,
	load_filtered_atoms_table,
)
from pdb_parser.io.topology import RESIDUE_BONDS


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


def test_load_filtered_atoms_table_preserves_rows_and_normalizes_headers(tmp_path: Path) -> None:
	input_file = write_filtered_atoms_file(tmp_path / "filtered_atoms.dat")

	df = load_filtered_atoms_table(input_file)

	assert list(df.columns) == ["atom_id", "atom_name", "resid", "resname", "x", "y", "z"]
	assert df.shape == (10, 7)
	assert df.iloc[0]["atom_name"] == "n"


def test_build_covalent_and_planar_topology_builds_expected_core_pairs(tmp_path: Path) -> None:
	input_file = write_filtered_atoms_file(tmp_path / "filtered_atoms.dat")

	topology = build_covalent_and_planar_topology(input_file)
	atom_index = topology["atom_index"]

	assert atom_index[(1, "N")] == 0
	assert atom_index[(1, "CA")] == 2
	assert atom_index[(1, "C")] == 4
	assert atom_index[(2, "N")] == 5

	bonded_pairs = topology["bonded_pairs"]
	assert (2, 0) in bonded_pairs
	assert (4, 2) in bonded_pairs
	assert (5, 4) in bonded_pairs

	angle_pairs = topology["angle_pairs"]
	assert (4, 0) in angle_pairs

	excluded_pairs = topology["excluded_pairs"]
	assert bonded_pairs.issubset(excluded_pairs)
	assert angle_pairs.issubset(excluded_pairs)


def test_covalent_and_planar_distances_writes_expected_pairs(tmp_path: Path) -> None:
	input_file = write_filtered_atoms_file(tmp_path / "filtered_atoms.dat")
	topology = build_covalent_and_planar_topology(input_file)

	output_file = tmp_path / "distance_constraints_1.dat"
	covalent_and_planar_distances(input_file, topology, output_file)

	lines = output_file.read_text(encoding="utf-8").splitlines()

	assert lines[0] == "atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j"
	assert "3\t1\t1\t1\t1.000000000000\t1.000000000000\tCA\tN\tALA\tALA" in lines
	assert "6\t5\t2\t1\t1.000000000000\t1.000000000000\tN\tC\tALA\tALA" in lines
	assert "5\t1\t1\t1\t2.000000000000\t2.000000000000\tC\tN\tALA\tALA" in lines


@pytest.mark.parametrize("resname", sorted(RESIDUE_BONDS))
def test_residue_bond_template_is_one_connected_molecule(resname: str) -> None:
	bonds = RESIDUE_BONDS[resname]
	atoms = {atom for bond in bonds for atom in bond}
	reached, frontier = {"N"}, ["N"]
	while frontier:
		atom = frontier.pop()
		for bond in bonds:
			if atom in bond:
				other = bond[1] if bond[0] == atom else bond[0]
				if other not in reached:
					reached.add(other)
					frontier.append(other)

	assert reached == atoms, f"{resname}: atoms not bonded to the backbone: {sorted(atoms - reached)}"


@pytest.mark.parametrize(
	("resname", "ring"),
	[
		("PHE", ("CG", "CD1", "CE1", "CZ", "CE2", "CD2")),
		("TYR", ("CG", "CD1", "CE1", "CZ", "CE2", "CD2")),
		("HIS", ("CG", "ND1", "CE1", "NE2", "CD2")),
		("TRP", ("CG", "CD1", "NE1", "CE2", "CD2")),
		("TRP", ("CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3")),
	],
)
def test_aromatic_rings_are_closed(resname: str, ring: tuple[str, ...]) -> None:
	bonds = {frozenset(bond) for bond in RESIDUE_BONDS[resname]}
	for atom_a, atom_b in zip(ring, ring[1:] + ring[:1]):
		assert frozenset((atom_a, atom_b)) in bonds, f"{resname}: missing {atom_a}-{atom_b}"


def test_tryptophan_indole_hydrogen_is_on_ne1() -> None:
	bonds = {frozenset(bond) for bond in RESIDUE_BONDS["TRP"]}
	assert frozenset(("NE1", "HE1")) in bonds
	assert frozenset(("CE2", "HE2")) not in bonds

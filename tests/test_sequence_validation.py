from pathlib import Path

import pandas as pd
import pytest

from pdb_parser.io import (
	backbone_angular_interval_to_distance_interval,
	build_covalent_and_planar_topology,
	ensure_model_chain_ready,
	extract_model_chain,
	load_filtered_atoms_table,
	nmr_distance_constraints,
	planar_peptide_distances,
	talos_n_like,
	validate_positional_cross_references,
	validate_residue_sequence,
	vdw_distances,
)
from pdb_parser.pipeline.parser_pipeline import parser


def test_residue_sequence_allows_non_one_start_and_rejects_gaps() -> None:
	assert validate_residue_sequence([10, 11, 12]) == [10, 11, 12]
	with pytest.raises(ValueError, match=r"10->12"):
		validate_residue_sequence([10, 12], context="chain A")


def test_filtered_atom_table_rejects_duplicate_atom_identity(tmp_path: Path) -> None:
	path = tmp_path / "filtered_atoms.dat"
	path.write_text(
		"atom_id\tatom_name\tresid\tresname\tx\ty\tz\n"
		"1\tCA\t10\tALA\t0\t0\t0\n"
		"2\tCA\t10\tALA\t1\t0\t0\n",
		encoding="utf-8",
	)
	with pytest.raises(ValueError, match="duplicate atom"):
		load_filtered_atoms_table(path)


def test_filtered_atom_table_rejects_missing_columns_and_nonfinite_coordinates(tmp_path: Path) -> None:
	missing = tmp_path / "missing.dat"
	missing.write_text("atom_id\tatom_name\tresid\tresname\tx\ty\n1\tN\t1\tALA\t0\t0\n", encoding="utf-8")
	with pytest.raises(ValueError, match="missing required columns.*z"):
		load_filtered_atoms_table(missing)

	nonfinite = tmp_path / "nonfinite.dat"
	nonfinite.write_text(
		"atom_id\tatom_name\tresid\tresname\tx\ty\tz\n1\tN\t1\tALA\t0\tinf\t0\n",
		encoding="utf-8",
	)
	with pytest.raises(ValueError, match="finite numeric"):
		load_filtered_atoms_table(nonfinite)


def test_topology_rejects_residue_without_a_covalent_template(tmp_path: Path) -> None:
	path = tmp_path / "unknown.dat"
	path.write_text(
		"atom_id\tatom_name\tresid\tresname\tx\ty\tz\n"
		"1\tN\t1\tMSE\t0\t0\t0\n"
		"2\tCA\t1\tMSE\t1\t0\t0\n",
		encoding="utf-8",
	)
	with pytest.raises(ValueError, match="without covalent topology templates.*MSE"):
		build_covalent_and_planar_topology(path)


def test_all_filtered_table_constraint_stages_reject_gaps(tmp_path: Path) -> None:
	path = tmp_path / "gap.dat"
	path.write_text(
		"atom_id\tatom_name\tresid\tresname\tx\ty\tz\n"
		"1\tN\t1\tALA\t0\t0\t0\n"
		"2\tCA\t1\tALA\t1\t0\t0\n"
		"3\tC\t1\tALA\t2\t0\t0\n"
		"4\tN\t3\tALA\t3\t0\t0\n"
		"5\tCA\t3\tALA\t4\t0\t0\n"
		"6\tC\t3\tALA\t5\t0\t0\n",
		encoding="utf-8",
	)
	with pytest.raises(ValueError, match=r"gaps are delimited by residues 1->3"):
		build_covalent_and_planar_topology(path)
	for writer, suffix, args in (
		(planar_peptide_distances, "planar", ("backbone",)),
		(nmr_distance_constraints, "nmr", ("precise", "backbone", 0, 0, 5, 0, 0, 0, 0)),
	):
		with pytest.raises(ValueError, match=r"1->3"):
			writer(path, tmp_path / f"{suffix}.dat", *args)
	with pytest.raises(ValueError, match=r"1->3"):
		backbone_angular_interval_to_distance_interval(
			path, tmp_path / "missing-angles.dat", tmp_path / "angular.dat"
		)
	with pytest.raises(ValueError, match=r"1->3"):
		vdw_distances(
			path,
			{"atom_index": {}, "residue_ids": (), "logical_graph": {},
			 "bonded_pairs": set(), "angle_pairs": set(), "excluded_pairs": set()},
			tmp_path / "vdw.dat",
		)


def test_unknown_residue_fixture_fails_pipeline_before_publishing_outputs(tmp_path: Path) -> None:
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	(pdb_dir / "UNKN.pdb").write_bytes(
		(Path(__file__).parent / "fixtures" / "unknown_residue.pdb").read_bytes()
	)
	params = {
		"method": "nmr", "output_mode": "parser_outputs", "model_number": "1",
		"chain_id": "A", "atom_selection": "full_chain", "distance_constraints": "precise",
		"max_distance": "5.0", "vdw_constraints": "no", "torsion_angle_width": "0.0",
		"percentage_backbone_torsion_angles": "100.0",
	}
	with pytest.raises(ValueError, match="without covalent topology templates.*MSE"):
		parser(params, pdb_dir, tmp_path / "seeds", tmp_path / "outputs", "UNKN")
	assert not list((tmp_path / "outputs" / "UNKN").glob("X_*.dat"))
	assert not list((tmp_path / "outputs" / "UNKN").glob("I_*.dat"))


def test_pdb_fixtures_cover_gap_insertion_altloc_and_explicit_chain_selection() -> None:
	fixture_dir = Path(__file__).parent / "fixtures"
	with pytest.raises(ValueError, match=r"10->12"):
		extract_model_chain(fixture_dir / "gap.pdb", 1, "A")
	inserted = extract_model_chain(fixture_dir / "insertion_code.pdb", 1, "A")
	assert list(inserted.residues.resids) == [10, 11]

	altloc = extract_model_chain(fixture_dir / "altloc.pdb", 1, "A")
	assert len(altloc) == 3
	assert float(altloc.select_atoms("name CA")[0].position[0]) == 1.0

	chain_b, model, method = ensure_model_chain_ready(
		fixture_dir / "multiple_chains.pdb", 1, "B", method="nmr"
	)
	assert model == 1 and method == "nmr"
	assert set(chain_b.residues.resnames) == {"GLY"}


def test_save_filtered_atoms_never_writes_duplicate_altloc_atoms(tmp_path: Path) -> None:
	import MDAnalysis as mda

	from pdb_parser.io.filtering import save_filtered_atoms

	u = mda.Universe(Path(__file__).parent / "fixtures" / "altloc.pdb")
	output = tmp_path / "filtered_atoms.dat"
	save_filtered_atoms(output, u.select_atoms("protein"), "full_chain", method="nmr")

	filtered = load_filtered_atoms_table(output)
	assert len(filtered) == 3
	assert float(filtered.loc[filtered["atom_name"] == "CA", "x"].iloc[0]) == 1.0


def test_parser_rejects_gap_before_creating_output_directory(tmp_path: Path) -> None:
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	(pdb_dir / "GAP1.pdb").write_bytes((Path(__file__).parent / "fixtures" / "gap.pdb").read_bytes())
	params = {
		"method": "nmr",
		"output_mode": "parser_outputs",
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone",
		"distance_constraints": "precise",
		"max_distance": "5.0",
		"vdw_constraints": "no",
		"torsion_angle_width": "0.0",
		"percentage_backbone_torsion_angles": "100.0",
	}
	with pytest.raises(ValueError, match=r"10->12"):
		parser(params, pdb_dir, tmp_path / "seeds", tmp_path / "outputs", "GAP1")
	assert not (tmp_path / "outputs").exists()


def test_talos_rejects_gap_before_writing_angular_artifact(tmp_path: Path) -> None:
	output = tmp_path / "A_gap.dat"
	with pytest.raises(ValueError, match=r"10->12"):
		talos_n_like(
			Path(__file__).parent / "fixtures" / "gap.pdb",
			1,
			"A",
			output,
			0.0,
			10.0,
			10.0,
			100.0,
		)
	assert not output.exists()


def test_positional_tables_require_consistent_x_a_i_references() -> None:
	df_x = pd.DataFrame([
		[1, "N", 1, "ALA", 0.0, 0.0, 0.0],
		[2, "CA", 1, "ALA", 1.0, 0.0, 0.0],
	], columns=range(7))
	df_a = pd.DataFrame([[1, "ALA", 0, 0, 0, 0, 0, 0]], columns=range(8))
	df_i = pd.DataFrame([[1, 99, 1, 1, 1.0, 1.0, "N", "CA", "ALA", "ALA"]], columns=range(10))
	with pytest.raises(ValueError, match="atom_id 99.*absent from X"):
		validate_positional_cross_references(df_x, df_a, df_i)

	df_i.iloc[0, 1] = 2
	df_i.iloc[0, 3] = 2
	df_i.iloc[0, 9] = "GLY"
	with pytest.raises(ValueError, match="inconsistent residue/atom metadata"):
		validate_positional_cross_references(df_x, df_a, df_i)


def test_fixture_multiple_chain_selection_reaches_pipeline_output(tmp_path: Path) -> None:
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	(pdb_dir / "MULT.pdb").write_bytes(
		(Path(__file__).parent / "fixtures" / "multiple_chains.pdb").read_bytes()
	)
	params = {
		"method": "nmr", "output_mode": "filtered_atoms", "model_number": "1",
		"chain_id": "B", "atom_selection": "backbone", "distance_constraints": "precise",
		"max_distance": "5.0", "vdw_constraints": "no", "torsion_angle_width": "0.0",
		"percentage_backbone_torsion_angles": "100.0",
	}
	parser(params, pdb_dir, tmp_path / "seeds", tmp_path / "outputs", "MULT")
	filtered = next((tmp_path / "outputs" / "MULT").glob("filtered_atoms_*.dat"))
	assert "GLY" in filtered.read_text(encoding="utf-8")


@pytest.mark.parametrize(("fixture", "message"), [
	("insertion_code.pdb", "insertion codes"),
	("altloc.pdb", "alternative"),
])
def test_parser_fixture_policies_are_enforced_before_outputs(
	tmp_path: Path, fixture: str, message: str
) -> None:
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	(pdb_dir / "CASE.pdb").write_bytes((Path(__file__).parent / "fixtures" / fixture).read_bytes())
	params = {
		"method": "nmr", "output_mode": "filtered_atoms", "model_number": "1",
		"chain_id": "A", "atom_selection": "backbone", "distance_constraints": "precise",
		"max_distance": "5.0", "vdw_constraints": "no", "torsion_angle_width": "0.0",
		"percentage_backbone_torsion_angles": "100.0",
	}
	if fixture == "altloc.pdb":
		parser(params, pdb_dir, tmp_path / "seeds", tmp_path / "outputs", "CASE")
		assert list((tmp_path / "outputs" / "CASE").glob("filtered_atoms_*.dat"))
	else:
		parser(params, pdb_dir, tmp_path / "seeds", tmp_path / "outputs", "CASE")
		assert list((tmp_path / "outputs" / "CASE").glob("filtered_atoms_*.dat"))


def _backbone_pdb(path: Path, residues: list[tuple[str, float]]) -> Path:
	"""Write N/CA/C for each (residue label, x of N) with 1.33 A bonds inside each residue."""
	lines = ["HEADER    PEPTIDE BOND FIXTURE", "EXPDTA    X-RAY DIFFRACTION"]
	serial = 1
	for label, x in residues:
		number, code = (label[:-1], label[-1]) if label[-1].isalpha() else (label, " ")
		for name, dx in (("N", 0.0), ("CA", 1.46), ("C", 2.0)):
			lines.append(
				f"ATOM  {serial:5d}  {name:<3s} ALA A{int(number):4d}{code}   "
				f"{x + dx:8.3f}{0.0:8.3f}{0.0:8.3f}{1.00:6.2f}{20.00:6.2f}           {name[0]}"
			)
			serial += 1
	lines.extend(("TER", "END", ""))
	path.write_text("\n".join(lines), encoding="utf-8")
	return path


def test_bonded_chain_passes_the_peptide_bond_check(tmp_path: Path) -> None:
	# C at x + 2.0 and the next N at x + 3.33: C-N = 1.33 A.
	pdb = _backbone_pdb(tmp_path / "ok.pdb", [("1", 0.0), ("2", 3.33), ("3", 6.66)])
	assert list(extract_model_chain(pdb, 1, "A").residues.resids) == [1, 2, 3]


def test_consecutive_numbers_without_a_peptide_bond_are_a_chain_break(tmp_path: Path) -> None:
	pdb = _backbone_pdb(tmp_path / "break.pdb", [("1", 0.0), ("2", 3.33), ("3", 10.0)])
	with pytest.raises(ValueError, match=r"chain A is broken between residues 2->3 \(4\.67 A\)"):
		extract_model_chain(pdb, 1, "A")
	assert list(extract_model_chain(pdb, 1, "A", allow_gaps=True).residues.resids) == [1, 2, 3]


def test_discarding_an_inserted_residue_that_breaks_the_chain_is_rejected(tmp_path: Path) -> None:
	# 1 - 1A - 2 is a continuous chain; keeping 1 and 2 only leaves 1 and 2 unbonded (1ORC-like).
	pdb = _backbone_pdb(tmp_path / "insertion.pdb", [("1", 0.0), ("1A", 3.33), ("2", 6.66)])
	with pytest.raises(ValueError, match=r"broken between residues 1->2.*insertion_code_policy"):
		extract_model_chain(pdb, 1, "A")

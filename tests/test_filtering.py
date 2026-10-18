from __future__ import annotations

from pdb_parser.io.filtering import _allowed_names_for_residue


def test_proline_hydrogen_selection_is_limited_to_backbone_plus_hydrogens() -> None:
	present = {"N", "CA", "C", "O", "CB", "CD", "H", "HA", "HD2", "HD3"}

	hydrogen_mode = _allowed_names_for_residue(
		"PRO",
		"backbone_plus_hydrogens",
		present,
		resid=1,
		method="nmr",
	)
	neighbor_mode = _allowed_names_for_residue(
		"PRO",
		"backbone_plus_neighbors",
		present,
		resid=1,
		method="nmr",
	)

	# one HN surrogate only, the one the remap orders (HD2 before HD3)
	assert "HD2" in hydrogen_mode and "HD3" not in hydrogen_mode
	assert "CD" not in hydrogen_mode
	assert "CD" in neighbor_mode
	assert "HD2" not in neighbor_mode
	assert "HD3" not in neighbor_mode


def _bph(resname, present, *, is_first=False):
	return _allowed_names_for_residue(
		resname, "backbone_plus_hydrogens", present, resid=1, method="nmr", is_first=is_first
	)


def test_backbone_plus_hydrogens_keeps_exactly_the_remapped_alternative() -> None:
	"""Regression: HA3 (Gly) and HD3 (Pro) were selected next to HA2/HD2, but the remap orders one
	atom per position, so the remap skipped the instance ("atom IDs absent from the remapping order")."""
	assert _bph("GLY", {"N", "CA", "C", "O", "H", "HA2", "HA3"}) == {"N", "CA", "C", "H", "HA2"}
	assert _bph("GLY", {"N", "CA", "C", "H", "HA3"}) == {"N", "CA", "C", "H", "HA3"}
	assert _bph("PRO", {"N", "CA", "C", "HA", "HD2", "HD3", "HB2"}) == {"N", "CA", "C", "HA", "HD2"}
	assert _bph("PRO", {"N", "CA", "C", "HA", "HD3"}) == {"N", "CA", "C", "HA", "HD3"}
	assert _bph("LYS", {"N", "CA", "C", "H", "HA", "HD2", "HD3"}) == {"N", "CA", "C", "H", "HA"}


def test_first_residue_follows_first_residue_order() -> None:
	from pdb_parser.remapping.first_residue import build_first_residue_pattern

	cases = {
		("GLY", frozenset({"N", "CA", "C", "H1", "H2", "H3", "HA2", "HA3"})): {"N", "CA", "C", "H1", "H2", "H3", "HA2"},
		("PRO", frozenset({"N", "CA", "C", "H2", "H3", "HA", "HD2", "HD3"})): {"N", "CA", "C", "H2", "H3", "HD2", "HA"},
		("MET", frozenset({"N", "CA", "C", "H", "HA"})): {"N", "CA", "C", "H", "HA"},
		("SER", frozenset({"N", "CA", "C", "H1", "H2", "HA"})): {"N", "CA", "C", "H1", "H2", "HA"},
	}
	for (resname, present), expected in cases.items():
		selected = _bph(resname, set(present), is_first=True)
		assert selected == expected, resname
		# the 5/6/7-atom clique pattern accepts the selection
		build_first_residue_pattern([(1, name) for name in selected], 1)


def test_nearest_proline_hd_follows_the_previous_alpha_hydrogen() -> None:
	from types import SimpleNamespace

	from pdb_parser.io.filtering import _nearest_proline_hd

	def atom(name, xyz):
		return SimpleNamespace(name=name, position=xyz)

	pro = [atom("HD2", (0.0, 0.0, 5.0)), atom("HD3", (0.0, 0.0, 3.0))]
	assert _nearest_proline_hd(pro, [atom("HA", (0.0, 0.0, 0.0))]) == "HD3"
	assert _nearest_proline_hd(pro, [atom("HA2", (0.0, 0.0, 8.0))]) == "HD2"  # GLY before PRO
	assert _nearest_proline_hd(pro, [atom("HA", (0.0, 0.0, 4.0))]) == "HD2"  # tie keeps HD2
	assert _nearest_proline_hd(pro[:1], [atom("HA", (0.0, 0.0, 0.0))]) is None
	assert _nearest_proline_hd(pro, [atom("N", (0.0, 0.0, 0.0))]) is None


def _parse_fixture(tmp_path, edit=None):
	from pathlib import Path

	from pdb_parser.pipeline.parser_pipeline import parser

	lines = (Path(__file__).parent / "fixtures" / "1ng7_model1_residues_1_6.pdb").read_text(encoding="utf-8").splitlines()
	(tmp_path / "pdb").mkdir()
	(tmp_path / "pdb" / "TEST.pdb").write_text("\n".join(edit(lines) if edit else lines) + "\n", encoding="utf-8")
	params = {
		"method": "nmr", "model_number": "1", "chain_id": "A", "atom_selection": "backbone_plus_hydrogens",
		"distance_constraints": "interval_experimental", "noe_strong": "2.5", "noe_medium": "3.5",
		"noe_weak": "5.0", "vdwr_hh": "1.8", "vdw_constraints": "no", "torsion_angle_width": "60.0",
		"percentage_backbone_torsion_angles": "100.0",
	}
	parser(params, tmp_path / "pdb", tmp_path / "seeds", tmp_path / "out", "TEST", True)
	x_file = tmp_path / "out" / "TEST" / "X_TEST_model1_chainA.dat"
	return [line.split() for line in x_file.read_text(encoding="utf-8").splitlines()]


def test_backbone_plus_hydrogens_keeps_the_proline_hd_nearest_the_previous_ha(tmp_path) -> None:
	"""1NG7 PRO 3: HD3 is 2.47 A from GLY 2 HA2 and HD2 2.85 A, so HD3 is kept."""
	rows = _parse_fixture(tmp_path)
	proline = {row[1] for row in rows if row[2] == "3"}

	assert "HD3" in proline and "HD2" not in proline


def test_atoms_with_repeated_coordinates_keep_the_first_one(tmp_path, capsys) -> None:
	"""1PIJ deposits H1 and H2 of residue 1 at the same position; H2 is dropped."""
	from pdb_parser.api import run_remap_job

	def copy_h1_onto_h2(lines):
		h1 = next(line for line in lines if line.startswith("ATOM") and line[22:26].strip() == "1" and line[12:16].strip() == "H1")
		return [
			line[:30] + h1[30:54] + line[54:]
			if line.startswith("ATOM") and line[22:26].strip() == "1" and line[12:16].strip() == "H2"
			else line
			for line in lines
		]

	rows = _parse_fixture(tmp_path, copy_h1_onto_h2)

	assert sorted(row[1] for row in rows if row[2] == "1") == ["C", "CA", "H1", "H3", "HA", "N"]
	assert "MET_1 H2 (same coordinates as MET_1 H1)" in capsys.readouterr().err
	for order_type in ("ddgp", "redmdgp"):
		params = {"model_number": "1", "chain_id": "A", "order_id": "1", "order_type": order_type}
		assert run_remap_job(params, tmp_path / "out", "TEST", 1) == 0

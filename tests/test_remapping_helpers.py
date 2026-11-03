from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pdb_parser.io.files import read_space_separated_file
from pdb_parser.remapping.dataframes import remap_distance_dataframe
from pdb_parser.remapping.vertex_ordering import (
	build_backbone_torsion_matrix,
	build_clique_pairs,
	build_ddgp_pattern_entries,
	build_first_residue_pattern,
	choose_first_existing_atom,
	detect_remap_atom_layout,
	get_backbone_residue_numeric_order,
	get_clique_pattern_ids,
	get_ha_group,
	get_internal_residue_ordering_template,
	get_nonterminal_h_group,
	map_clique_type,
	resolve_atom_name_choice,
)


def test_get_nonterminal_h_group_and_ha_group_handle_special_residues() -> None:
	assert get_nonterminal_h_group("PRO") == ("HD2", "HD3")
	assert get_nonterminal_h_group("ALA") == ("H",)
	assert get_ha_group("GLY") == ("HA2", "HA3")
	assert get_ha_group("SER") == ("HA",)


def test_get_internal_residue_ordering_template_uses_special_atom_aliases() -> None:
	template = get_internal_residue_ordering_template(9, "GLY")
	assert template == [("H",), ("CA",), ("N",), ("HA2", "HA3"), ("C",)]

	pro_template = get_internal_residue_ordering_template(1, "PRO")
	assert pro_template == [("N",), ("CA",), ("C",), ("HD2", "HD3"), ("HA",)]


def test_resolve_atom_name_choice_prefers_first_available_candidate() -> None:
	atom_name_to_atom_id = {"HA3": 7, "CA": 2}

	assert resolve_atom_name_choice(atom_name_to_atom_id, ("HA2", "HA3"), 4, "GLY") == (7, "HA3")


def test_choose_first_existing_atom_honors_exclusions() -> None:
	atom_names_present = {"H3", "H2", "H1"}

	assert choose_first_existing_atom(atom_names_present, ("H3", "H2", "H1"), {"H3", "H2"}) == "H1"
	assert choose_first_existing_atom(atom_names_present, ("HD2", "HD3")) is None


def test_get_clique_pattern_ids_and_map_clique_type_cover_known_patterns() -> None:
	assert get_clique_pattern_ids(1) == (1, 2, 3, 4, 5)
	assert map_clique_type(1) == 3
	assert map_clique_type(6) == 2
	assert map_clique_type(13) == 1


def test_build_clique_pairs_resolves_aliases_from_available_atoms() -> None:
	available_atoms = [
		(4, "C"),
		(5, "N"),
		(5, "CA"),
		(5, "HA2"),
		(5, "HD2"),
	]

	pairs = build_clique_pairs(residue_id=5, clique_pattern_id=13, available_atoms=available_atoms)

	assert pairs == [(4, "C"), (5, "HD2"), (5, "CA"), (5, "N")]


def test_build_ddgp_pattern_entries_returns_five_entries() -> None:
	available_atoms = [
		(4, "CA"),
		(4, "C"),
		(5, "N"),
		(5, "CA"),
		(5, "C"),
		(5, "HA"),
		(5, "H"),
	]

	entries = build_ddgp_pattern_entries(ddgp_pattern_id=1, residue_id=5, available_atoms=available_atoms)

	assert len(entries) == 5
	assert entries[0][1] == 3
	assert entries[-1][1] == 1


def test_build_first_residue_pattern_handles_five_atom_case() -> None:
	available_atoms = [
		(1, "N"),
		(1, "CA"),
		(1, "C"),
		(1, "HA"),
		(1, "H"),
	]

	rows = build_first_residue_pattern(available_atoms, 1)

	assert len(rows) == 5
	assert rows[-1] == ([(1, "HA"), (1, "CA"), (1, "N"), (1, "H")], 2)


@pytest.mark.parametrize(
	("n_terminal", "expected_order"),
	[
		(("H3", "H2", "H1"), ("H3", "H2", "H1")),
		(("H3", "H", "H1"), ("H3", "H", "H1")),
		(("H", "H2", "H1"), ("H", "H2", "H1")),
	],
)
def test_build_first_residue_pattern_follows_the_x_order_of_every_seven_atom_pattern(
	n_terminal: tuple[str, ...], expected_order: tuple[str, ...]
) -> None:
	"""Patterns 2 (H3 h H1) and 3 (h H2 H1) used to look for H2/H3 by name and fail."""
	available_atoms = [(1, name) for name in ("N", "CA", "C", "HA", *n_terminal)]

	rows = build_first_residue_pattern(available_atoms, 1)
	h1, h2, h3 = expected_order

	assert [row[0][-1] for row in rows] == [(1, h1), (1, h2), (1, h3), (1, "N"), (1, "CA"), (1, "HA"), (1, "C")]
	assert rows[5] == ([(1, h3), (1, "N"), (1, "CA"), (1, "HA")], 2)


@pytest.mark.parametrize("rename", [{"H2": "H"}, {"H3": "H"}])
def test_ddgp_remap_accepts_seven_atom_n_termini_without_h2_or_h3(tmp_path, rename: dict[str, str]) -> None:
	from pathlib import Path

	import numpy as np

	from pdb_parser.api import run_remap_job
	from pdb_parser.io.files import read_space_separated_file
	from pdb_parser.output_layout import remapped_clique_file
	from pdb_parser.pipeline.parser_pipeline import parser

	lines = []
	for line in (Path(__file__).parent / "fixtures" / "1ng7_model1_residues_1_6.pdb").read_text(encoding="utf-8").splitlines():
		name = line[12:16].strip()
		if line.startswith("ATOM") and line[22:26].strip() == "1" and name in rename:
			line = line[:12] + f"{rename[name]:^4s}" + line[16:]
		lines.append(line)
	(tmp_path / "pdb").mkdir()
	(tmp_path / "pdb" / "TEST.pdb").write_text("\n".join(lines) + "\n", encoding="utf-8")
	params = {
		"method": "nmr", "model_number": "1", "chain_id": "A", "atom_selection": "backbone_plus_hydrogens",
		"distance_constraints": "interval_experimental", "noe_strong": "2.5", "noe_medium": "3.5",
		"noe_weak": "5.0", "vdwr_hh": "1.8", "vdw_constraints": "no", "torsion_angle_width": "60.0",
		"percentage_backbone_torsion_angles": "100.0",
	}
	parser(params, tmp_path / "pdb", tmp_path / "seeds", tmp_path / "out", "TEST", True)

	for order_id in (1, 9):
		assert run_remap_job({"model_number": "1", "chain_id": "A", "order_id": str(order_id)}, tmp_path / "out", "TEST", order_id) == 0
		T = read_space_separated_file(remapped_clique_file(tmp_path / "out", "TEST", 1, "A", order_id))
		assert (T[0] == np.arange(1, len(T) + 1)).all()


def test_build_first_residue_pattern_rejects_unsupported_atom_count() -> None:
	with pytest.raises(ValueError, match="Unsupported number of atoms for residue 1"):
		build_first_residue_pattern([(1, "N"), (1, "CA"), (1, "C"), (1, "HA")], 1)


def test_remap_distance_dataframe_reports_unmapped_atom_ids() -> None:
	import pandas as pd

	distances = pd.DataFrame([[3, 1, 2, 1, 1.0, 2.0, "CA", "N", "GLY", "ALA"]])
	with pytest.raises(ValueError, match=r"\[3\]"):
		remap_distance_dataframe(distances, [1, 2])


def test_detect_remap_atom_layout_accepts_backbone_only() -> None:
	df_X = read_space_separated_file(
		Path(__file__).parent / "fixtures" / "backbone_structure.dat"
	)

	assert detect_remap_atom_layout(df_X, [1, 2, 3]) == "backbone"


def test_get_backbone_residue_numeric_order_returns_n_ca_c_order() -> None:
	df_X = read_space_separated_file(
		Path(__file__).parent / "fixtures" / "backbone_structure.dat"
	)
	df_X_i = df_X[df_X[2] == 2]

	assert get_backbone_residue_numeric_order(df_X_i) == [
		(4, "N"),
		(5, "CA"),
		(6, "C"),
	]


def test_build_backbone_torsion_matrix_uses_previous_three_vertices_and_parser_angles() -> None:
	df_X = read_space_separated_file(
		Path(__file__).parent / "fixtures" / "backbone_structure.dat"
	)
	df_A = read_space_separated_file(
		Path(__file__).parent / "fixtures" / "backbone_angles.dat"
	)

	T = build_backbone_torsion_matrix(df_X, df_A)

	assert np.array_equal(T[0], np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
	assert np.array_equal(T[1], np.array([2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
	assert np.array_equal(T[2], np.array([3.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0]))
	assert np.array_equal(T[3], np.array([4.0, 3.0, 2.0, 1.0, -1.0, 70.0, 3.0]))
	assert np.array_equal(T[4], np.array([5.0, 4.0, 3.0, 2.0, -1.0, 120.0, 4.0]))
	assert np.array_equal(T[5], np.array([6.0, 5.0, 4.0, 3.0, 1.0, 80.0, 5.0]))
	assert np.array_equal(T[6], np.array([7.0, 6.0, 5.0, 4.0, 1.0, 90.0, 6.0]))


def test_missing_discretization_interval_is_a_data_issue_naming_the_atoms() -> None:
	import numpy as np
	import pandas as pd

	from pdb_parser.remapping.torsions import update_torsion_matrix

	df_X = pd.DataFrame([
		[1, "HA", 65, "ASN", 0.0, 0.0, 0.0],
		[2, "CA", 65, "ASN", 1.0, 0.0, 0.0],
		[3, "C", 65, "ASN", 1.0, 1.0, 0.0],
		[4, "H", 66, "GLN", 2.0, 1.0, 1.0],
	])
	df_I = pd.DataFrame([[1, 2, 65, 65, 1.0, 1.0], [2, 3, 65, 65, 1.0, 1.0]])
	T = np.zeros((4, 7))
	T[3, :5] = (4, 3, 2, 1, 2)  # interval clique: needs the HA_65--H_66 distance

	with pytest.raises(ValueError, match=r"Interval not found for atom pair \(4, 1\).*Atoms: H_66 and HA_65"):
		update_torsion_matrix(T, df_X, np.zeros((0, 8)), df_I)

from __future__ import annotations

from pathlib import Path

import pytest

from pdb_parser.output_layout import (
	parser_angular_file,
	parser_distance_file,
	parser_structure_file,
	pdb_output_dir,
	remapped_clique_file,
	remapped_distance_file,
	remapped_output_dir,
	remapped_structure_file,
)


def test_parser_output_layout_matches_current_contract(tmp_path: Path) -> None:
	out_dir = tmp_path / "outputs"

	assert pdb_output_dir(out_dir, "1ABC") == out_dir / "1ABC"
	assert parser_structure_file(out_dir, "1ABC", 1, "A") == out_dir / "1ABC" / "X_1ABC_model1_chainA.dat"
	assert parser_angular_file(out_dir, "1ABC", 1, "A") == out_dir / "1ABC" / "A_1ABC_model1_chainA.dat"
	assert parser_distance_file(out_dir, "1ABC", 1, "A") == out_dir / "1ABC" / "I_1ABC_model1_chainA.dat"


def test_remapped_output_layout_matches_current_contract(tmp_path: Path) -> None:
	out_dir = tmp_path / "outputs"

	assert remapped_output_dir(out_dir, "1ABC") == out_dir / "1ABC" / "remapped"
	assert remapped_structure_file(out_dir, "1ABC", 1, "A", 9) == (
		out_dir / "1ABC" / "remapped" / "X_1ABC_model1_chainA_ddgp_order_9.dat"
	)
	assert remapped_distance_file(out_dir, "1ABC", 1, "A", 9) == (
		out_dir / "1ABC" / "remapped" / "I_1ABC_model1_chainA_ddgp_order_9.dat"
	)
	assert remapped_clique_file(out_dir, "1ABC", 1, "A", 9) == (
		out_dir / "1ABC" / "remapped" / "T_1ABC_model1_chainA_ddgp_order_9.dat"
	)


def test_remapped_output_layout_rejects_non_numeric_order(tmp_path: Path) -> None:
	with pytest.raises(ValueError, match="order_id"):
		remapped_structure_file(tmp_path / "outputs", "1ABC", 1, "A", "random")


def test_remapped_output_layout_uses_dmdgp_suffix_for_backbone(tmp_path: Path) -> None:
	out_dir = tmp_path / "outputs"

	assert remapped_structure_file(out_dir, "1ABC", 1, "A", 9, "backbone") == (
		out_dir / "1ABC" / "remapped" / "X_1ABC_model1_chainA_dmdgp_order.dat"
	)
	assert remapped_distance_file(out_dir, "1ABC", 1, "A", 9, "backbone") == (
		out_dir / "1ABC" / "remapped" / "I_1ABC_model1_chainA_dmdgp_order.dat"
	)
	assert remapped_clique_file(out_dir, "1ABC", 1, "A", 1, "backbone") == (
		out_dir / "1ABC" / "remapped" / "T_1ABC_model1_chainA_dmdgp_order.dat"
	)


@pytest.mark.parametrize(
	("pdb_id", "model", "chain", "order"),
	[
		("../ABCD", 1, "A", 1),
		("1ABC", 0, "A", 1),
		("1ABC", 1, "../A", 1),
		("1ABC", 1, "A", 0),
	],
)
def test_output_layout_rejects_ambiguous_components(
	tmp_path: Path, pdb_id: str, model: int, chain: str, order: int
) -> None:
	with pytest.raises(ValueError):
		remapped_structure_file(tmp_path, pdb_id, model, chain, order)

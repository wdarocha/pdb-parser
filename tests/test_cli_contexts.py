from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from pdb_parser.instance_remap import load_remap_cli_context
from pdb_parser.pdb_parser import load_parser_cli_context


def write_text_file(path: Path, content: str) -> Path:
	path.write_text(content, encoding="utf-8")
	return path


def test_load_parser_cli_context_reads_ids_params_and_creates_dirs(tmp_path: Path) -> None:
	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1abc\n")
	params_file = write_text_file(
		tmp_path / "params.cfg",
		"""
		model_number: 1
		chain_id: A
		atom_selection: backbone_plus_hydrogens
		distance_constraints: precise
		max_distance: 5
		vdw_constraints: yes
		torsion_angle_width: 0
		percentage_backbone_torsion_angles: 100
		""",
	)

	args = argparse.Namespace(
		pdb_ids_file=pdb_ids_file,
		params_file=params_file,
		pdb_data_dir=tmp_path / "pdb",
		seed_dir=tmp_path / "seeds",
		out_dir=tmp_path / "outputs",
	)

	context = load_parser_cli_context(args)

	assert context["pdb_ids"] == ["1ABC"]
	assert context["params"]["chain_id"] == "A"
	assert Path(context["pdb_data_dir"]).exists()
	assert Path(context["seed_dir"]).exists()
	assert Path(context["out_dir"]).exists()


def test_load_parser_cli_context_accepts_backbone_public_config(tmp_path: Path) -> None:
	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1abc\n")
	params_file = write_text_file(
		tmp_path / "params.cfg",
		"""
		model_number: 1
		chain_id: A
		atom_selection: backbone
		distance_constraints: precise
		max_distance: 5
		vdw_constraints: yes
		torsion_angle_width: 0
		percentage_backbone_torsion_angles: 100
		""",
	)

	args = argparse.Namespace(
		pdb_ids_file=pdb_ids_file,
		params_file=params_file,
		pdb_data_dir=tmp_path / "pdb",
		seed_dir=tmp_path / "seeds",
		out_dir=tmp_path / "outputs",
	)

	context = load_parser_cli_context(args)

	assert context["params"]["atom_selection"] == "backbone"


def test_load_parser_cli_context_accepts_full_chain_atom_selection(tmp_path: Path) -> None:
	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1abc\n")
	params_file = write_text_file(
		tmp_path / "params.cfg",
		"""
		model_number: 1
		chain_id: A
		atom_selection: full_chain
		distance_constraints: precise
		max_distance: 5
		vdw_constraints: yes
		torsion_angle_width: 0
		percentage_backbone_torsion_angles: 100
		""",
	)

	args = argparse.Namespace(
		pdb_ids_file=pdb_ids_file,
		params_file=params_file,
		pdb_data_dir=tmp_path / "pdb",
		seed_dir=tmp_path / "seeds",
		out_dir=tmp_path / "outputs",
	)

	context = load_parser_cli_context(args)

	assert context["params"]["atom_selection"] == "full_chain"


def test_load_remap_cli_context_reads_public_inputs(tmp_path: Path) -> None:
	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1abc\n")
	remap_params_file = write_text_file(
		tmp_path / "instance_remap.cfg",
		"""
		model_number: 1
		chain_id: A
		order_id: 9
		""",
	)
	out_dir = tmp_path / "outputs"
	out_dir.mkdir()

	args = argparse.Namespace(
		pdb_ids_file=pdb_ids_file,
		remap_params_file=remap_params_file,
		out_dir=out_dir,
	)

	context = load_remap_cli_context(args)

	assert context["pdb_ids"] == ["1ABC"]
	assert context["remap_params"]["order_id"] == "9"
	assert context["order_choices"] == [9]
	assert context["order_types"] == ["ddgp"]
	assert context["out_dir"] == out_dir


def test_load_remap_cli_context_rejects_non_numeric_order(tmp_path: Path) -> None:
	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1abc\n")
	remap_params_file = write_text_file(
		tmp_path / "instance_remap.cfg",
		"""
		model_number: 1
		chain_id: A
		order_id: random
		""",
	)
	out_dir = tmp_path / "outputs"
	out_dir.mkdir()

	args = argparse.Namespace(
		pdb_ids_file=pdb_ids_file,
		remap_params_file=remap_params_file,
		out_dir=out_dir,
	)

	with pytest.raises(ValueError, match="order_id"):
		load_remap_cli_context(args)


def test_load_remap_cli_context_rejects_missing_output_dir(tmp_path: Path) -> None:
	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1abc\n")
	remap_params_file = write_text_file(
		tmp_path / "instance_remap.cfg",
		"""
		model_number: 1
		chain_id: A
		order_id: 9
		""",
	)

	args = argparse.Namespace(
		pdb_ids_file=pdb_ids_file,
		remap_params_file=remap_params_file,
		out_dir=tmp_path / "missing_outputs",
	)

	with pytest.raises(ValueError, match="Output directory does not exist"):
		load_remap_cli_context(args)

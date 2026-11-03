"""Several remap orders in one call: ``order_id: 1-10`` and ``order_type: ddgp, redmdgp``."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from test_redmdgp import _parse, _write_fixture

import pdb_parser.api as api
from pdb_parser.api import run_remap_job, run_remap_jobs
from pdb_parser.public_config import (
	parse_public_order_choices,
	parse_public_order_types,
	validate_public_remap_params,
)

BASE = {"model_number": "1", "chain_id": "A"}


@pytest.mark.parametrize(
	("raw", "expected"),
	[("9", [9]), (" 1-10 ", list(range(1, 11))), ("1,3,9", [1, 3, 9]), ("1-4, 9, 3", [1, 2, 3, 4, 9])],
)
def test_order_id_accepts_one_value_lists_and_ranges(raw: str, expected: list[int]) -> None:
	assert parse_public_order_choices(raw) == expected


@pytest.mark.parametrize(
	("raw", "message"),
	[
		("11", "Invalid order_id"),
		("0-3", "Invalid order_id"),
		("3-1", "the start must not exceed the end"),
		("1,,2", "empty item"),
		("1-x", "must be an integer"),
	],
)
def test_order_id_lists_are_validated(raw: str, message: str) -> None:
	with pytest.raises(ValueError, match=message):
		parse_public_order_choices(raw)


def test_order_type_accepts_a_comma_separated_list() -> None:
	assert parse_public_order_types({}) == ["ddgp"]
	assert parse_public_order_types({"order_type": "ddgp, REDMDGP, ddgp"}) == ["ddgp", "redmdgp"]
	with pytest.raises(ValueError, match="Invalid order_type 'ddmgp'"):
		parse_public_order_types({"order_type": "ddgp, ddmgp"})
	with pytest.raises(ValueError, match="dmdgp cannot be combined"):
		parse_public_order_types({"order_type": "ddgp, dmdgp"})
	validate_public_remap_params({**BASE, "order_id": "1-10", "order_type": "ddgp, redmdgp"}, Path("x.cfg"))


def _remapped_files(out_dir: Path) -> dict[str, bytes]:
	return {path.name: path.read_bytes() for path in sorted((out_dir / "TEST" / "remapped").glob("*.dat"))}


def test_one_call_writes_the_same_files_as_one_call_per_order(tmp_path: Path) -> None:
	_write_fixture(tmp_path / "pdb", ("H1", "H2", "H3"))
	out_dir = _parse(tmp_path)

	for order_type in ("ddgp", "redmdgp"):
		for order_id in range(1, 11):
			params = {**BASE, "order_id": str(order_id), "order_type": order_type}
			assert run_remap_job(params, out_dir, "TEST", order_id) == 0
	separate = _remapped_files(out_dir)
	shutil.rmtree(out_dir / "TEST" / "remapped")

	run_remap_jobs({**BASE, "order_id": "1-10", "order_type": "ddgp, redmdgp"}, out_dir, "TEST",
		list(range(1, 11)), ["ddgp", "redmdgp"])

	assert len(separate) == 60
	assert _remapped_files(out_dir) == separate


def test_a_failed_order_does_not_stop_the_others(
	tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
	_write_fixture(tmp_path / "pdb", ("H1", "H2", "H3"))
	out_dir = _parse(tmp_path)
	real_run_remap_job = api.run_remap_job

	def failing_order_3(params, out, pdb_id, order_id, **kwargs):
		if order_id == 3:
			raise ValueError("order 3 is broken")
		return real_run_remap_job(params, out, pdb_id, order_id, **kwargs)

	monkeypatch.setattr(api, "run_remap_job", failing_order_3)
	with pytest.raises(ValueError, match=r"1 of 4 remap orders failed \(ddgp 3\)"):
		run_remap_jobs({**BASE, "order_id": "1-4"}, out_dir, "TEST", [1, 2, 3, 4], ["ddgp"])

	assert "TEST: order_type ddgp, order_id 3 failed: order 3 is broken" in capsys.readouterr().err
	assert sorted(name.split("_order_")[1] for name in _remapped_files(out_dir) if name.startswith("T_")) == [
		"1.dat", "2.dat", "4.dat",
	]


def test_backbone_outputs_are_written_once_for_any_number_of_orders(
	tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
	_write_fixture(tmp_path / "pdb", ("H1", "H2", "H3"))
	out_dir = _parse(tmp_path, atom_selection="backbone")

	run_remap_jobs({**BASE, "order_id": "1-10", "order_type": "ddgp, redmdgp"}, out_dir, "TEST",
		list(range(1, 11)), ["ddgp", "redmdgp"])

	assert sorted(_remapped_files(out_dir)) == [
		"I_TEST_model1_chainA_dmdgp_order.dat",
		"T_TEST_model1_chainA_dmdgp_order.dat",
		"X_TEST_model1_chainA_dmdgp_order.dat",
	]
	assert capsys.readouterr().err.count("order_id and order_type are ignored") == 1


def test_order_type_dmdgp_makes_order_id_optional() -> None:
	validate_public_remap_params({**BASE, "order_type": "dmdgp"}, Path("x.cfg"))
	validate_public_remap_params({**BASE, "order_type": "dmdgp", "order_id": "9"}, Path("x.cfg"))
	with pytest.raises(ValueError, match="Invalid order_id"):
		validate_public_remap_params({**BASE, "order_type": "dmdgp", "order_id": "11"}, Path("x.cfg"))
	with pytest.raises(ValueError, match="order_id"):
		validate_public_remap_params({**BASE, "order_type": "ddgp"}, Path("x.cfg"))


def _run_remap_cli(tmp_path: Path, out_dir: Path, config: str) -> int:
	from pdb_parser import instance_remap

	ids = tmp_path / "ids.txt"
	ids.write_text("TEST\n", encoding="utf-8")
	cfg = tmp_path / "remap.cfg"
	cfg.write_text(config, encoding="utf-8")
	return instance_remap.main([str(ids), str(cfg), str(out_dir)])


def test_order_type_dmdgp_writes_the_backbone_order_without_warning(
	tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
	_write_fixture(tmp_path / "pdb", ("H1", "H2", "H3"))
	out_dir = _parse(tmp_path, atom_selection="backbone")

	assert _run_remap_cli(tmp_path, out_dir, "model_number: 1\nchain_id: A\norder_type: dmdgp\n") == 0
	explicit = _remapped_files(out_dir)
	assert "order_id and order_type are ignored" not in capsys.readouterr().err

	shutil.rmtree(out_dir / "TEST" / "remapped")
	assert _run_remap_cli(tmp_path, out_dir, "model_number: 1\nchain_id: A\norder_id: 9\n") == 0
	assert "use order_type: dmdgp" in capsys.readouterr().err
	assert sorted(explicit) == [
		"I_TEST_model1_chainA_dmdgp_order.dat",
		"T_TEST_model1_chainA_dmdgp_order.dat",
		"X_TEST_model1_chainA_dmdgp_order.dat",
	]
	assert _remapped_files(out_dir) == explicit


def test_order_type_dmdgp_is_rejected_for_backbone_plus_hydrogens_outputs(
	tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
	_write_fixture(tmp_path / "pdb", ("H1", "H2", "H3"))
	out_dir = _parse(tmp_path)

	assert _run_remap_cli(tmp_path, out_dir, "model_number: 1\nchain_id: A\norder_type: dmdgp\n") == 1
	assert "order_type: dmdgp requires parser outputs with atom_selection: backbone" in capsys.readouterr().out
	assert not (out_dir / "TEST" / "remapped").exists() or not list((out_dir / "TEST" / "remapped").iterdir())

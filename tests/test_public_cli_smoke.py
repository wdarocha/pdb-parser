from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from pdb_parser import instance_remap as instance_remap_module
from pdb_parser.instance_remap import build_parser as build_remap_parser
from pdb_parser.pdb_parser import build_parser as build_parser_parser


def test_parser_cli_help_contract() -> None:
	help_text = build_parser_parser().format_help()

	assert "pdb-parser" in help_text
	assert "pdb_ids_file" in help_text
	assert "params_file" in help_text
	assert "seed_dir" in help_text
	assert "out_dir" in help_text



def test_instance_remap_cli_help_contract() -> None:
	help_text = build_remap_parser().format_help()

	assert "instance-remap" in help_text
	assert "remap_params_file" in help_text
	assert "out_dir" in help_text
	assert "Remap previously generated parser output files." in help_text



def test_run_pipeline_help_contract() -> None:
	repo_root = Path(__file__).resolve().parents[1]
	script_path = repo_root / "run_pipeline.sh"
	env = os.environ.copy()
	env["PYTHONPATH"] = f"{repo_root / 'src'}{os.pathsep}{env.get('PYTHONPATH', '')}".rstrip(os.pathsep)

	result = subprocess.run(
		[str(script_path), "--help"],
		check=True,
		cwd=repo_root,
		env=env,
		capture_output=True,
		text=True,
	)

	help_text = result.stderr
	assert "instance_remap.cfg" in help_text
	assert "remap-all" in help_text
	assert "--mode <mode>" in help_text
	assert "Execution mode: full, parser, remap, or remap-all" in help_text


def test_instance_remap_cli_warns_when_backbone_order_is_ignored(
	tmp_path: Path,
	capsys,
) -> None:
	pdb_ids_file = tmp_path / "pdb_ids.txt"
	pdb_ids_file.write_text("1ABC\n", encoding="utf-8")
	remap_params_file = tmp_path / "instance_remap.cfg"
	remap_params_file.write_text(
		"model_number: 1\nchain_id: A\norder_id: 9\n",
		encoding="utf-8",
	)
	out_dir = tmp_path / "outputs"
	pdb_dir = out_dir / "1ABC"
	pdb_dir.mkdir(parents=True)
	(pdb_dir / "X_1ABC_model1_chainA.dat").write_text(
		"1 N 1 ALA 0.0 0.0 0.0\n"
		"2 CA 1 ALA 1.0 0.0 0.0\n"
		"3 C 1 ALA 1.0 1.0 0.0\n"
		"4 N 2 GLY 2.0 1.0 0.0\n"
		"5 CA 2 GLY 2.0 2.0 0.0\n"
		"6 C 2 GLY 3.0 2.0 0.0\n"
		"7 N 3 SER 4.0 2.0 0.0\n"
		"8 CA 3 SER 4.0 3.0 0.0\n"
		"9 C 3 SER 5.0 3.0 0.0\n",
		encoding="utf-8",
	)
	(pdb_dir / "A_1ABC_model1_chainA.dat").write_text(
		"1 ALA 170.0 1.0 60.0 2.0 -70.0 3.0\n"
		"2 GLY -120.0 4.0 80.0 5.0 90.0 6.0\n"
		"3 SER 110.0 7.0 -40.0 8.0 130.0 9.0\n",
		encoding="utf-8",
	)
	(pdb_dir / "I_1ABC_model1_chainA.dat").write_text(
		"9 1 3 1 1.0 1.0 C N SER ALA\n",
		encoding="utf-8",
	)

	result = instance_remap_module.main(
		[
			str(pdb_ids_file),
			str(remap_params_file),
			str(out_dir),
		]
	)

	captured = capsys.readouterr()

	assert result == 0
	assert "configured order_id and order_type are ignored for this layout" in captured.err


def test_instance_remap_cli_returns_one_when_a_structure_is_skipped(
	tmp_path: Path,
	monkeypatch: pytest.MonkeyPatch,
	capsys,
) -> None:
	pdb_ids_file = tmp_path / "pdb_ids.txt"
	pdb_ids_file.write_text("1ABC\n2ABC\n", encoding="utf-8")
	remap_params_file = tmp_path / "instance_remap.cfg"
	remap_params_file.write_text("model_number: 1\nchain_id: A\norder_id: 9\n", encoding="utf-8")
	out_dir = tmp_path / "outputs"
	out_dir.mkdir()
	# 1ABC is skipped (missing hydrogens, for example); 2ABC is remapped.
	monkeypatch.setattr(
		"pdb_parser.api.run_remap_job",
		lambda params, out, pdb_id, order: 1 if pdb_id == "1ABC" else 0,
	)

	result = instance_remap_module.main([str(pdb_ids_file), str(remap_params_file), str(out_dir)])

	assert result == 1
	assert "[SKIP] 1ABC: Data issue: structure skipped by the remap stage" in capsys.readouterr().out


def test_parser_cli_exposes_keep_tmp_flag() -> None:
	assert "--keep-tmp" in build_parser_parser().format_help()


@pytest.mark.parametrize("keep_tmp", [True, False])
def test_parser_cli_keeps_tmp_dir_only_with_keep_tmp(
	tmp_path: Path, keep_tmp: bool
) -> None:
	from test_pipeline_integration import _write_complete_nmr_pdb

	from pdb_parser import pdb_parser as parser_module

	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	_write_complete_nmr_pdb(pdb_dir / "TEST.pdb")
	(tmp_path / "pdb_ids.txt").write_text("TEST\n", encoding="utf-8")
	(tmp_path / "params.cfg").write_text(
		"method: nmr\nmodel_number: 1\nchain_id: A\natom_selection: backbone_plus_hydrogens\n"
		"distance_constraints: precise\nmax_distance: 5.0\nvdw_constraints: no\n"
		"torsion_angle_width: 0.0\npercentage_backbone_torsion_angles: 100.0\n",
		encoding="utf-8",
	)
	# No network: download_pdb reuses the valid fixture already cached in pdb_dir.
	out_dir = tmp_path / "outputs"

	argv = [*(["--keep-tmp"] if keep_tmp else []), str(tmp_path / "pdb_ids.txt"), str(tmp_path / "params.cfg"),
		str(pdb_dir), str(tmp_path / "seeds"), str(out_dir)]
	assert parser_module.main(argv) == 0

	assert (out_dir / "TEST" / "I_TEST_model1_chainA.dat").is_file()
	assert (out_dir / "TEST" / "tmp").is_dir() is keep_tmp

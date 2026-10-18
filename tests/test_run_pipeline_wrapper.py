from __future__ import annotations

import os
import stat
import subprocess
from pathlib import Path

import pytest


def write_text_file(path: Path, content: str) -> Path:
	path.write_text(content, encoding="utf-8")
	return path


def make_python_stub(tmp_path: Path) -> tuple[Path, Path]:
	stub_dir = tmp_path / "stub_calls"
	stub_dir.mkdir()
	stub_path = tmp_path / "python_stub.sh"
	stub_path.write_text(
		"""#!/usr/bin/env bash
set -eu
counter_file="$STUB_DIR/counter"
count=0
if [ -f "$counter_file" ]; then
	count="$(cat "$counter_file")"
fi
count=$((count + 1))
printf '%s' "$count" > "$counter_file"
printf '%s\n' "$@" > "$STUB_DIR/call${count}.args"
cat "$4" > "$STUB_DIR/call${count}.cfg"
""",
		encoding="utf-8",
	)
	stub_path.chmod(stub_path.stat().st_mode | stat.S_IEXEC)
	return stub_path, stub_dir


def test_run_pipeline_remap_all_runs_every_listed_order_in_one_call(tmp_path: Path) -> None:
	repo_root = Path(__file__).resolve().parents[1]
	script_path = repo_root / "run_pipeline.sh"
	stub_path, stub_dir = make_python_stub(tmp_path)

	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1ABC\n")
	parser_params_file = write_text_file(tmp_path / "params.cfg", "model_number: 1\n")
	remap_params_file = write_text_file(
		tmp_path / "instance_remap.cfg",
		"model_number: 1\nchain_id: A\norder_id: 9\n",
	)
	pdb_data_dir = tmp_path / "pdb"
	seed_dir = tmp_path / "seeds"
	out_dir = tmp_path / "outputs"
	pdb_data_dir.mkdir()
	seed_dir.mkdir()
	out_dir.mkdir()

	env = os.environ.copy()
	env["PYTHON"] = str(stub_path)
	env["STUB_DIR"] = str(stub_dir)

	subprocess.run(
		[
			str(script_path),
			"--mode",
			"remap-all",
			"--orders",
			"1,7",
			str(pdb_ids_file),
			str(parser_params_file),
			str(remap_params_file),
			str(pdb_data_dir),
			str(seed_dir),
			str(out_dir),
		],
		check=True,
		cwd=repo_root,
		env=env,
	)

	# One remap call receives every listed order.
	assert (stub_dir / "counter").read_text(encoding="utf-8") == "1"
	call_1_cfg = (stub_dir / "call1.cfg").read_text(encoding="utf-8")
	call_1_args = (stub_dir / "call1.args").read_text(encoding="utf-8")

	assert "order_id: 1,7" in call_1_cfg
	assert "order_id: 9" not in call_1_cfg
	assert "pdb_parser.instance_remap" in call_1_args


def test_run_pipeline_rejects_invalid_mode(tmp_path: Path) -> None:
	repo_root = Path(__file__).resolve().parents[1]
	script_path = repo_root / "run_pipeline.sh"

	result = subprocess.run(
		[str(script_path), "--mode", "bad-mode"],
		cwd=repo_root,
		capture_output=True,
		text=True,
	)

	assert result.returncode == 2
	assert "Invalid mode: bad-mode" in result.stderr


def test_run_pipeline_rejects_remap_when_filtered_atoms_was_requested(tmp_path: Path) -> None:
	repo_root = Path(__file__).resolve().parents[1]
	script_path = repo_root / "run_pipeline.sh"
	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1ABC\n")
	parser_params_file = write_text_file(tmp_path / "params.cfg", "output_mode: filtered_atoms\n")
	remap_params_file = write_text_file(tmp_path / "instance_remap.cfg", "order_id: 1\n")

	result = subprocess.run(
		[
			str(script_path), "--mode", "full", str(pdb_ids_file), str(parser_params_file),
			str(remap_params_file), str(tmp_path / "pdb"), str(tmp_path / "seeds"), str(tmp_path / "outputs"),
		],
		cwd=repo_root,
		capture_output=True,
		text=True,
	)

	assert result.returncode == 2
	assert "requires parser output_mode: parser_outputs" in result.stderr


def test_run_pipeline_rejects_invalid_orders_spec(tmp_path: Path) -> None:
	repo_root = Path(__file__).resolve().parents[1]
	script_path = repo_root / "run_pipeline.sh"
	stub_path, stub_dir = make_python_stub(tmp_path)

	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1ABC\n")
	parser_params_file = write_text_file(tmp_path / "params.cfg", "model_number: 1\n")
	remap_params_file = write_text_file(
		tmp_path / "instance_remap.cfg",
		"model_number: 1\nchain_id: A\norder_id: 9\n",
	)
	pdb_data_dir = tmp_path / "pdb"
	seed_dir = tmp_path / "seeds"
	out_dir = tmp_path / "outputs"
	pdb_data_dir.mkdir()
	seed_dir.mkdir()
	out_dir.mkdir()

	env = os.environ.copy()
	env["PYTHON"] = str(stub_path)
	env["STUB_DIR"] = str(stub_dir)

	result = subprocess.run(
		[
			str(script_path),
			"--mode",
			"remap-all",
			"--orders",
			"1,abc",
			str(pdb_ids_file),
			str(parser_params_file),
			str(remap_params_file),
			str(pdb_data_dir),
			str(seed_dir),
			str(out_dir),
		],
		cwd=repo_root,
		capture_output=True,
		text=True,
		env=env,
	)

	assert result.returncode == 2
	assert "Invalid order specification: abc" in result.stderr


def test_run_pipeline_rejects_missing_required_arguments(tmp_path: Path) -> None:
	repo_root = Path(__file__).resolve().parents[1]
	script_path = repo_root / "run_pipeline.sh"

	result = subprocess.run(
		[str(script_path), "--mode", "parser"],
		cwd=repo_root,
		capture_output=True,
		text=True,
	)

	assert result.returncode == 2
	assert "Usage:" in result.stderr


@pytest.mark.parametrize(
	("arguments", "message"),
	[
		(["--mode", "invalid"], "Invalid mode: invalid"),
		(["--mode", "remap-all", "--orders", "3-1"], "Invalid order range: 3-1"),
		(["--mode", "remap-all", "--orders", "x"], "Invalid order specification: x"),
	]
)
def test_run_pipeline_invalid_input_matrix_is_safe_with_existing_output_dir(
	tmp_path: Path,
	arguments: list[str],
	message: str,
) -> None:
	repo_root = Path(__file__).resolve().parents[1]
	script_path = repo_root / "run_pipeline.sh"
	out_dir = tmp_path / "existing-output"
	out_dir.mkdir()
	marker = out_dir / "keep.txt"
	marker.write_text("do not remove", encoding="utf-8")
	(tmp_path / "ids.txt").write_text("1ABC\n", encoding="utf-8")
	(tmp_path / "params.cfg").write_text("output_mode: parser_outputs\n", encoding="utf-8")
	(tmp_path / "remap.cfg").write_text("order_id: 1\n", encoding="utf-8")
	base = [
		str(script_path),
		*arguments,
		str(tmp_path / "ids.txt"),
		str(tmp_path / "params.cfg"),
		str(tmp_path / "remap.cfg"),
		str(tmp_path / "pdb"),
		str(tmp_path / "seeds"),
		str(out_dir),
	]
	result = subprocess.run(base, cwd=repo_root, capture_output=True, text=True)

	assert result.returncode == 2
	assert message in result.stderr
	assert marker.read_text(encoding="utf-8") == "do not remove"


def test_run_pipeline_remap_mode_warns_and_writes_dmdgp_outputs_for_backbone(tmp_path: Path) -> None:
	repo_root = Path(__file__).resolve().parents[1]
	script_path = repo_root / "run_pipeline.sh"

	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1ABC\n")
	parser_params_file = write_text_file(tmp_path / "params.cfg", "model_number: 1\n")
	remap_params_file = write_text_file(
		tmp_path / "instance_remap.cfg",
		"model_number: 1\nchain_id: A\norder_id: 9\n",
	)
	pdb_data_dir = tmp_path / "pdb"
	seed_dir = tmp_path / "seeds"
	out_dir = tmp_path / "outputs"
	pdb_data_dir.mkdir()
	seed_dir.mkdir()
	pdb_dir = out_dir / "1ABC"
	pdb_dir.mkdir(parents=True)

	write_text_file(
		pdb_dir / "X_1ABC_model1_chainA.dat",
		(
			"1 N 1 ALA 0.0 0.0 0.0\n"
			"2 CA 1 ALA 1.0 0.0 0.0\n"
			"3 C 1 ALA 1.0 1.0 0.0\n"
			"4 N 2 GLY 2.0 1.0 0.0\n"
			"5 CA 2 GLY 2.0 2.0 0.0\n"
			"6 C 2 GLY 3.0 2.0 0.0\n"
			"7 N 3 SER 4.0 2.0 0.0\n"
			"8 CA 3 SER 4.0 3.0 0.0\n"
			"9 C 3 SER 5.0 3.0 0.0\n"
		),
	)
	write_text_file(
		pdb_dir / "A_1ABC_model1_chainA.dat",
		(
			"1 ALA 170.0 1.0 60.0 2.0 -70.0 3.0\n"
			"2 GLY -120.0 4.0 80.0 5.0 90.0 6.0\n"
			"3 SER 110.0 7.0 -40.0 8.0 130.0 9.0\n"
		),
	)
	write_text_file(
		pdb_dir / "I_1ABC_model1_chainA.dat",
		"9 1 3 1 1.0 1.0 C N SER ALA\n",
	)

	env = os.environ.copy()
	env["PYTHONPATH"] = f"{repo_root / 'src'}{os.pathsep}{env.get('PYTHONPATH', '')}".rstrip(os.pathsep)

	result = subprocess.run(
		[
			str(script_path),
			"--mode",
			"remap",
			str(pdb_ids_file),
			str(parser_params_file),
			str(remap_params_file),
			str(pdb_data_dir),
			str(seed_dir),
			str(out_dir),
		],
		check=True,
		cwd=repo_root,
		env=env,
		capture_output=True,
		text=True,
	)

	remapped_dir = out_dir / "1ABC" / "remapped"

	assert "configured order_id and order_type are ignored for this layout" in result.stderr
	assert (remapped_dir / "X_1ABC_model1_chainA_dmdgp_order.dat").exists()
	assert (remapped_dir / "I_1ABC_model1_chainA_dmdgp_order.dat").exists()
	assert (remapped_dir / "T_1ABC_model1_chainA_dmdgp_order.dat").exists()


def test_run_pipeline_defaults_to_full_mode(tmp_path: Path) -> None:
	repo_root = Path(__file__).resolve().parents[1]
	script_path = repo_root / "run_pipeline.sh"
	stub_path, stub_dir = make_python_stub(tmp_path)
	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1ABC\n")
	parser_params_file = write_text_file(tmp_path / "params.cfg", "model_number: 1\n")
	remap_params_file = write_text_file(tmp_path / "instance_remap.cfg", "order_id: 9\n")
	env = os.environ.copy()
	env["PYTHON"] = str(stub_path)
	env["STUB_DIR"] = str(stub_dir)

	subprocess.run(
		[
			str(script_path), str(pdb_ids_file), str(parser_params_file), str(remap_params_file),
			str(tmp_path / "pdb"), str(tmp_path / "seeds"), str(tmp_path / "outputs"),
		],
		check=True,
		cwd=repo_root,
		env=env,
	)

	assert "pdb_parser.pdb_parser" in (stub_dir / "call1.args").read_text(encoding="utf-8")
	assert "pdb_parser.instance_remap" in (stub_dir / "call2.args").read_text(encoding="utf-8")


@pytest.mark.parametrize("keep_tmp", [True, False])
def test_run_pipeline_passes_keep_tmp_only_to_the_parser(tmp_path: Path, keep_tmp: bool) -> None:
	repo_root = Path(__file__).resolve().parents[1]
	script_path = repo_root / "run_pipeline.sh"
	stub_path, stub_dir = make_python_stub(tmp_path)
	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "1ABC\n")
	parser_params_file = write_text_file(tmp_path / "params.cfg", "model_number: 1\n")
	remap_params_file = write_text_file(tmp_path / "instance_remap.cfg", "order_id: 9\n")
	env = os.environ.copy()
	env["PYTHON"] = str(stub_path)
	env["STUB_DIR"] = str(stub_dir)

	subprocess.run(
		[
			str(script_path), *(["--keep-tmp"] if keep_tmp else []),
			str(pdb_ids_file), str(parser_params_file), str(remap_params_file),
			str(tmp_path / "pdb"), str(tmp_path / "seeds"), str(tmp_path / "outputs"),
		],
		check=True,
		cwd=repo_root,
		env=env,
	)

	parser_args = (stub_dir / "call1.args").read_text(encoding="utf-8").splitlines()
	remap_args = (stub_dir / "call2.args").read_text(encoding="utf-8").splitlines()
	assert ("--keep-tmp" in parser_args) is keep_tmp
	assert "--keep-tmp" not in remap_args

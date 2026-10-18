from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pdb_parser.api import (
	BACKBONE_REMAP_ORDER_IGNORED_WARNING,
	build_ddgp_order_vector,
	run_parser_job,
	run_remap_job,
)
from pdb_parser.cli_support import run_cli_jobs
from pdb_parser.output_layout import (
	parser_angular_file,
	parser_distance_file,
	parser_structure_file,
	remapped_clique_file,
	remapped_distance_file,
	remapped_structure_file,
)
from pdb_parser.pipeline.remapping import write_remapped_artifacts


def test_run_parser_job_calls_download_then_parser(tmp_path: Path) -> None:
	calls: list[tuple] = []

	def fake_download(pdb_id: str, pdb_data_dir: Path) -> None:
		calls.append(("download", pdb_id, Path(pdb_data_dir)))

	def fake_parser(
		params: dict[str, str],
		pdb_data_dir: Path,
		seed_dir: Path,
		out_dir: Path,
		pdb_id: str,
		remove_tmp_dir: bool,
	) -> None:
		calls.append(
			("parser", params, Path(pdb_data_dir), Path(seed_dir), Path(out_dir), pdb_id, remove_tmp_dir)
		)

	params = {"model_number": "1"}
	run_parser_job(
		params=params,
		pdb_data_dir=tmp_path / "pdb",
		seed_dir=tmp_path / "seeds",
		out_dir=tmp_path / "outputs",
		pdb_id="1ABC",
		parser_fn=fake_parser,
		download_pdb_fn=fake_download,
		remove_tmp_dir=True,
	)

	assert calls[0] == ("download", "1ABC", tmp_path / "pdb")
	assert calls[1] == (
		"parser",
		params,
		tmp_path / "pdb",
		tmp_path / "seeds",
		tmp_path / "outputs",
		"1ABC",
		True,
	)


def test_run_parser_job_normalizes_pdb_id_and_preserves_download_failure(
	tmp_path: Path,
) -> None:
	calls: list[str] = []

	def fake_download(pdb_id: str, data_dir: Path) -> None:
		calls.append(pdb_id)
		raise ValueError("download failed")

	with pytest.raises(ValueError, match="download failed"):
		run_parser_job(
		params={}, pdb_data_dir=tmp_path / "pdb", seed_dir=tmp_path / "seeds",
		out_dir=tmp_path / "outputs", pdb_id="1abc", parser_fn=lambda *args: None,
		 download_pdb_fn=fake_download,
	)
	assert calls == ["1ABC"]


def test_run_parser_job_rejects_path_separator_in_pdb_id(tmp_path: Path) -> None:
	with pytest.raises(ValueError, match="Invalid PDB id"):
		run_parser_job(
			params={}, pdb_data_dir=tmp_path, seed_dir=tmp_path, out_dir=tmp_path,
			pdb_id="../1ABC", parser_fn=lambda *args: None,
		)


def test_run_remap_job_builds_vector_and_calls_remap_instance(tmp_path: Path) -> None:
	out_dir = tmp_path / "outputs"
	pdb_dir = out_dir / "1ABC"
	pdb_dir.mkdir(parents=True)
	(pdb_dir / "X_1ABC_model1_chainA.dat").write_text(
		"1 N 1 ALA 0.0 0.0 0.0\n"
		"2 CA 1 ALA 1.0 0.0 0.0\n"
		"3 C 1 ALA 1.0 1.0 0.0\n"
		"4 H 1 ALA 0.0 0.0 1.0\n"
		"5 HA 1 ALA 0.0 1.0 0.0\n"
		"6 N 2 GLY 2.0 1.0 0.0\n"
		"7 CA 2 GLY 2.0 2.0 0.0\n"
		"8 C 2 GLY 3.0 2.0 0.0\n"
		"9 H 2 GLY 2.0 1.0 1.0\n"
		"10 HA2 2 GLY 2.0 2.0 1.0\n"
		"11 N 3 SER 4.0 2.0 0.0\n"
		"12 CA 3 SER 4.0 3.0 0.0\n"
		"13 C 3 SER 5.0 3.0 0.0\n"
		"14 H 3 SER 4.0 2.0 1.0\n"
		"15 HA 3 SER 4.0 3.0 1.0\n",
		encoding="utf-8",
	)

	calls: list[tuple] = []

	def fake_build_vector(
		remap_params: dict[str, str],
		out_dir: Path,
		pdb_id: str,
		order_choice: int,
	) -> np.ndarray:
		calls.append(("build_vector", remap_params, Path(out_dir), pdb_id, order_choice))
		return np.array([9, 9, 9], dtype=int)

	def fake_remap_instance(
		remap_params: dict[str, str],
		out_dir: Path,
		pdb_id: str,
		ddgp_order_vec: np.ndarray,
	) -> int:
		calls.append(("remap", remap_params, Path(out_dir), pdb_id, ddgp_order_vec.copy()))
		return 0

	remap_params = {"model_number": "1", "chain_id": "A", "order_id": "9"}
	result = run_remap_job(
		remap_params=remap_params,
		out_dir=out_dir,
		pdb_id="1ABC",
		order_choice=9,
		remap_instance_fn=fake_remap_instance,
		build_ddgp_order_vector_fn=fake_build_vector,
	)

	assert result == 0
	assert calls[0] == ("build_vector", remap_params, out_dir, "1ABC", 9)
	assert calls[1][0:4] == ("remap", remap_params, out_dir, "1ABC")
	assert np.array_equal(calls[1][4], np.array([9, 9, 9], dtype=int))


def test_run_remap_job_rejects_x_gap_before_reading_dependencies_and_preserves_outputs(
	tmp_path: Path,
) -> None:
	out_dir = tmp_path / "outputs"
	pdb_dir = out_dir / "1ABC"
	pdb_dir.mkdir(parents=True)
	(pdb_dir / "X_1ABC_model1_chainA.dat").write_text(
		"1 N 1 ALA 0 0 0\n2 CA 1 ALA 1 0 0\n3 C 1 ALA 2 0 0\n"
		"4 N 3 ALA 3 0 0\n5 CA 3 ALA 4 0 0\n6 C 3 ALA 5 0 0\n",
		encoding="utf-8",
	)
	(pdb_dir / "A_1ABC_model1_chainA.dat").write_text(
		"1 ALA 0 0 0 0 0 0\n3 ALA 0 0 0 0 0 0\n", encoding="utf-8"
	)
	(pdb_dir / "I_1ABC_model1_chainA.dat").write_text(
		"3 4 1 3 1 1 C N ALA ALA\n", encoding="utf-8"
	)
	remapped_dir = pdb_dir / "remapped"
	remapped_dir.mkdir()
	existing = {
		remapped_structure_file(out_dir, "1ABC", 1, "A", 1, "backbone"),
		remapped_distance_file(out_dir, "1ABC", 1, "A", 1, "backbone"),
		remapped_clique_file(out_dir, "1ABC", 1, "A", 1, "backbone"),
	}
	for path in existing:
		path.write_text("previous generation\n", encoding="utf-8")

	with pytest.raises(ValueError, match="Remap cannot proceed because residue gaps"):
		run_remap_job({"model_number": "1", "chain_id": "A", "order_id": "1"}, out_dir, "1ABC", 1)
	assert all(path.read_text(encoding="utf-8") == "previous generation\n" for path in existing)


def test_remap_gap_does_not_abort_the_following_pdb_job(tmp_path: Path) -> None:
	out_dir = tmp_path / "outputs"
	pdb_dir = out_dir / "1ABC"
	pdb_dir.mkdir(parents=True)
	(pdb_dir / "X_1ABC_model1_chainA.dat").write_text(
		"1 N 1 ALA 0 0 0\n2 CA 1 ALA 1 0 0\n3 C 1 ALA 2 0 0\n"
		"4 N 3 ALA 3 0 0\n5 CA 3 ALA 4 0 0\n6 C 3 ALA 5 0 0\n",
		encoding="utf-8",
	)
	completed: list[str] = []

	def run_job(pdb_id: str) -> None:
		if pdb_id == "1ABC":
			run_remap_job(
				{"model_number": "1", "chain_id": "A", "order_id": "1"},
				out_dir,
				pdb_id,
				1,
			)
		completed.append(pdb_id)

	assert run_cli_jobs(["1ABC", "2ABC"], run_job) == 1
	assert completed == ["2ABC"]


def test_run_remap_job_skips_order_vector_and_warns_for_backbone_layout(
	tmp_path: Path,
	capsys: pytest.CaptureFixture[str],
) -> None:
	out_dir = tmp_path / "outputs"
	pdb_dir = out_dir / "1ABC"
	pdb_dir.mkdir(parents=True)
	(pdb_dir / "X_1ABC_model1_chainA.dat").write_text(
		"1 N 1 ALA 0.0 0.0 0.0\n"
		"2 CA 1 ALA 1.0 0.0 0.0\n"
		"3 C 1 ALA 1.0 1.0 0.0\n"
		"4 N 2 GLY 2.0 1.0 0.0\n"
		"5 CA 2 GLY 2.0 2.0 0.0\n"
		"6 C 2 GLY 3.0 2.0 0.0\n",
		encoding="utf-8",
	)

	calls: list[tuple] = []

	def fake_build_vector(
		remap_params: dict[str, str],
		out_dir: Path,
		pdb_id: str,
		order_choice: int | str,
	) -> np.ndarray:
		calls.append(("build_vector", remap_params, Path(out_dir), pdb_id, order_choice))
		return np.array([9, 9], dtype=int)

	def fake_remap_instance(
		remap_params: dict[str, str],
		out_dir: Path,
		pdb_id: str,
		ddgp_order_vec: np.ndarray,
	) -> int:
		calls.append(("remap", remap_params, Path(out_dir), pdb_id, ddgp_order_vec.copy()))
		return 0

	result = run_remap_job(
		remap_params={"model_number": "1", "chain_id": "A", "order_id": "9"},
		out_dir=out_dir,
		pdb_id="1ABC",
		order_choice=9,
		remap_instance_fn=fake_remap_instance,
		build_ddgp_order_vector_fn=fake_build_vector,
	)

	captured = capsys.readouterr()

	assert result == 0
	assert len(calls) == 1
	assert calls[0][0:4] == ("remap", {"model_number": "1", "chain_id": "A", "order_id": "9"}, out_dir, "1ABC")
	assert np.array_equal(calls[0][4], np.zeros(2, dtype=int))
	assert BACKBONE_REMAP_ORDER_IGNORED_WARNING in captured.err






def test_build_ddgp_order_vector_reports_missing_structure_file_for_remapping(tmp_path: Path) -> None:
	out_dir = tmp_path / "outputs"
	(out_dir / "1ABC").mkdir(parents=True)

	try:
		build_ddgp_order_vector(
			remap_params={"model_number": "1", "chain_id": "A", "order_id": "9"},
			out_dir=out_dir,
			pdb_id="1ABC",
			order_id=9,
		)
	except ValueError as exc:
		assert "Parser structure file required by the remapping stage was not found" in str(exc)
	else:
		raise AssertionError("Expected ValueError for missing X_* file")




def test_run_remap_job_writes_outputs_for_fixed_order(tmp_path: Path) -> None:
	out_dir = tmp_path / "outputs"
	pdb_dir = out_dir / "1ABC"
	pdb_dir.mkdir(parents=True)

	(pdb_dir / "X_1ABC_model1_chainA.dat").write_text(
		"1 N 1 ALA 0.0 0.0 0.0\n"
		"2 CA 1 ALA 1.0 0.0 0.0\n"
		"3 C 1 ALA 1.0 1.0 0.0\n"
		"4 HA 1 ALA 0.0 1.0 0.0\n"
		"5 H 1 ALA 0.0 0.0 1.0\n"
		"6 N 2 GLY 2.0 1.0 0.0\n"
		"7 CA 2 GLY 2.0 2.0 0.0\n"
		"8 C 2 GLY 3.0 2.0 0.0\n"
		"9 HA2 2 GLY 2.0 2.0 1.0\n"
		"10 H 2 GLY 2.0 1.0 1.0\n",
		encoding="utf-8",
	)
	(pdb_dir / "A_1ABC_model1_chainA.dat").write_text(
		"1 ALA 0 0 0 10 0 10\n"
		"2 GLY 0 0 0 10 0 10\n",
		encoding="utf-8",
	)
	(pdb_dir / "I_1ABC_model1_chainA.dat").write_text(
		"5 4 1 1 2.4 2.4 H HA ALA ALA\n"
		"10 4 2 1 2.4 2.4 H HA GLY ALA\n"
		"10 5 2 1 2.4 2.4 H H GLY ALA\n"
		"9 4 2 1 2.4 2.4 HA2 HA GLY ALA\n"
		"9 5 2 1 2.4 2.4 HA2 H GLY ALA\n"
		"9 10 2 2 2.4 2.4 HA2 H GLY GLY\n",
		encoding="utf-8",
	)

	remap_params = {"model_number": "1", "chain_id": "A", "order_id": "9"}
	result = run_remap_job(
		remap_params=remap_params,
		out_dir=out_dir,
		pdb_id="1ABC",
		order_choice=9,
	)

	xfile = remapped_structure_file(out_dir, "1ABC", 1, "A", 9)
	ifile = remapped_distance_file(out_dir, "1ABC", 1, "A", 9)
	tfile = remapped_clique_file(out_dir, "1ABC", 1, "A", 9)

	assert result == 0
	assert xfile.exists()
	assert ifile.exists()
	assert tfile.exists()
	assert xfile.read_text(encoding="utf-8").strip()
	assert ifile.read_text(encoding="utf-8").strip()
	assert tfile.read_text(encoding="utf-8").strip()
	# Exercise the extracted remapping stage directly for the
	# backbone_plus_hydrogens layout as well as through the public API.
	assert write_remapped_artifacts(
		remap_params,
		out_dir,
		"1ABC",
		np.array([9, 9], dtype=int),
	) == 0


def test_run_remap_job_writes_outputs_for_backbone_parser_layout(
	tmp_path: Path,
) -> None:
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

	result = run_remap_job(
		remap_params={"model_number": "1", "chain_id": "A", "order_id": "9"},
		out_dir=out_dir,
		pdb_id="1ABC",
		order_choice=9,
	)

	tfile = remapped_clique_file(out_dir, "1ABC", 1, "A", 9, "backbone")

	assert result == 0
	assert tfile.exists()
	assert tfile.read_text(encoding="utf-8").splitlines()[3] == "    4     3     2     1 -1  70.000000000000   3.000000000000"

	# The backbone path intentionally ignores the order vector, but must remain
	# runnable through the extracted remapping stage.
	assert write_remapped_artifacts(
		{"model_number": "1", "chain_id": "A", "order_id": "9"},
		out_dir,
		"1ABC",
		np.array([9, 9, 9], dtype=int),
	) == 0


def test_remap_gap_fails_before_reading_dependents_and_preserves_existing_artifact(tmp_path: Path) -> None:
	out_dir = tmp_path / "outputs"
	pdb_dir = out_dir / "1ABC"
	pdb_dir.mkdir(parents=True)
	(parser_structure_file(out_dir, "1ABC", 1, "A")).write_text(
		"1 N 10 ALA 0 0 0\n2 CA 10 ALA 1 0 0\n3 C 12 ALA 2 0 0\n",
		encoding="utf-8",
	)
	parser_angular_file(out_dir, "1ABC", 1, "A").write_text(
		"10 ALA 0 0 0 0 0 0\n12 ALA 0 0 0 0 0 0\n", encoding="utf-8"
	)
	parser_distance_file(out_dir, "1ABC", 1, "A").write_text(
	"3 1 12 10 1 1 C N ALA ALA\n", encoding="utf-8"
	)
	existing = remapped_structure_file(out_dir, "1ABC", 1, "A", 1, "backbone")
	existing.parent.mkdir(parents=True, exist_ok=True)
	existing.write_text("sentinel\n", encoding="utf-8")

	with pytest.raises(ValueError, match="Remap cannot proceed because residue gaps"):
		run_remap_job(
			{"model_number": "1", "chain_id": "A", "order_id": "1"},
			out_dir,
			"1ABC",
			1,
		)

	assert existing.read_text(encoding="utf-8") == "sentinel\n"

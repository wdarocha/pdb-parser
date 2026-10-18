"""End-to-end contracts for the extracted parser and remapping stages."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pdb_parser.api import run_remap_job
from pdb_parser.io.files import read_space_separated_file
from pdb_parser.output_layout import (
    parser_angular_file,
    parser_distance_file,
    parser_metadata_file,
    parser_structure_file,
    remapped_clique_file,
    remapped_distance_file,
    remapped_structure_file,
)
from pdb_parser.pipeline import artifacts as artifacts_pipeline
from pdb_parser.pipeline import remapping as remapping_pipeline
from pdb_parser.pipeline.parser_pipeline import parser


def _write_complete_nmr_pdb(path: Path) -> None:
    """Create a small complete chain that exercises both public atom layouts."""
    lines = ["HEADER    INTEGRATION FIXTURE", "EXPDTA    SOLUTION NMR"]
    serial = 1
    for residue_id in range(1, 4):
        x = float(residue_id * 3)
        atoms = (
            ("N", x, 0.0, 0.0, "N"),
            ("CA", x + 1.0, 1.0, 0.0, "C"),
            ("C", x + 2.0, 0.0, 1.0, "C"),
            ("H", x, -1.0, 0.0, "H"),
            ("HA", x + 1.0, 2.0, 0.0, "H"),
        )
        for atom_name, atom_x, atom_y, atom_z, element in atoms:
            lines.append(
                f"ATOM  {serial:5d} {atom_name:^4s} ALA A{residue_id:4d}    "
                f"{atom_x:8.3f}{atom_y:8.3f}{atom_z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}"
            )
            serial += 1
    lines.extend(("TER", "END", ""))
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_complete_xray_pdb(path: Path) -> None:
	_write_complete_nmr_pdb(path)
	path.write_text(
		path.read_text(encoding="utf-8").replace("EXPDTA    SOLUTION NMR", "EXPDTA    X-RAY DIFFRACTION"),
		encoding="utf-8",
	)


@pytest.mark.parametrize(
    ("atom_selection", "order_id"),
    [("backbone", 1), ("backbone_plus_hydrogens", 5)],
)
def test_real_parser_and_remap_preserve_xiat_contracts(
    tmp_path: Path,
    atom_selection: str,
    order_id: int,
) -> None:
    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()
    _write_complete_nmr_pdb(pdb_dir / "TEST.pdb")
    out_dir = tmp_path / "outputs"
    params = {
        "method": "nmr",
        "output_mode": "parser_outputs",
        "model_number": "1",
        "chain_id": "A",
        "atom_selection": atom_selection,
        "distance_constraints": "precise",
        "max_distance": "5.0",
        "vdw_constraints": "no",
        "torsion_angle_width": "0.0",
        "percentage_backbone_torsion_angles": "100.0",
    }

    parser(params, pdb_dir, tmp_path / "seeds", out_dir, "TEST", True)
    remap_params = {"model_number": "1", "chain_id": "A", "order_id": str(order_id)}
    assert run_remap_job(remap_params, out_dir, "TEST", order_id) == 0

    xfile = parser_structure_file(out_dir, "TEST", 1, "A")
    afile = parser_angular_file(out_dir, "TEST", 1, "A")
    ifile = parser_distance_file(out_dir, "TEST", 1, "A")
    x_frame = read_space_separated_file(xfile)
    i_frame = read_space_separated_file(ifile)
    a_frame = read_space_separated_file(afile)
    metadata = json.loads(
        parser_metadata_file(out_dir, "TEST", 1, "A").read_text(encoding="utf-8")
    )

    assert metadata["atom_selection"] == atom_selection
    assert len(x_frame) == (9 if atom_selection == "backbone" else 15)
    assert {0, 1, 2, 3, 4, 5, 6}.issubset(x_frame.columns)
    assert {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}.issubset(i_frame.columns)
    assert len(a_frame) == 3

    layout = atom_selection
    remapped = (
        remapped_structure_file(out_dir, "TEST", 1, "A", order_id, layout),
        remapped_distance_file(out_dir, "TEST", 1, "A", order_id, layout),
        remapped_clique_file(out_dir, "TEST", 1, "A", order_id, layout),
    )
    assert all(path.exists() and path.read_text(encoding="utf-8").strip() for path in remapped)
    assert np.loadtxt(remapped[0], ndmin=2).shape[0] == len(x_frame)


@pytest.mark.parametrize("atom_selection", ["full_chain", "backbone_plus_neighbors"])
def test_real_parser_accepts_all_other_public_atom_selections(
	tmp_path: Path,
	atom_selection: str,
) -> None:
	"""The public parser modes must produce normal parser-stage artifacts."""
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	_write_complete_nmr_pdb(pdb_dir / "TEST.pdb")
	out_dir = tmp_path / "outputs"
	params = {
		"method": "nmr",
		"output_mode": "parser_outputs",
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": atom_selection,
		"distance_constraints": "precise",
		"max_distance": "5.0",
		"vdw_constraints": "no",
		"torsion_angle_width": "0.0",
		"percentage_backbone_torsion_angles": "100.0",
	}

	parser(params, pdb_dir, tmp_path / "seeds", out_dir, "TEST", True)
	metadata = json.loads(
		parser_metadata_file(out_dir, "TEST", 1, "A").read_text(encoding="utf-8")
	)
	x_frame = read_space_separated_file(parser_structure_file(out_dir, "TEST", 1, "A"))

	assert metadata["method"] == "nmr"
	assert metadata["atom_selection"] == atom_selection
	assert len(x_frame) == 15
	assert parser_angular_file(out_dir, "TEST", 1, "A").exists()
	assert parser_distance_file(out_dir, "TEST", 1, "A").exists()


def test_xray_parser_falls_back_from_hydrogens_and_preserves_output_contract(
	tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	_write_complete_xray_pdb(pdb_dir / "XRAY.pdb")
	out_dir = tmp_path / "outputs"
	params = {
		"method": "xray",
		"output_mode": "parser_outputs",
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone_plus_hydrogens",
		"distance_constraints": "precise",
		"max_distance": "20.0",
		"vdw_constraints": "no",
		"torsion_angle_width": "0.0",
		"percentage_backbone_torsion_angles": "100.0",
	}

	parser(params, pdb_dir, tmp_path / "seeds", out_dir, "XRAY", True)
	metadata = json.loads(parser_metadata_file(out_dir, "XRAY", 1, "A").read_text())
	assert metadata["method"] == "xray"
	assert metadata["atom_selection"] == "backbone_plus_neighbors"
	assert parser_structure_file(out_dir, "XRAY", 1, "A").exists()
	assert parser_distance_file(out_dir, "XRAY", 1, "A").exists()
	assert "unavailable for X-ray" in capsys.readouterr().err


def test_xray_backbone_parser_output_can_be_remapped(tmp_path: Path) -> None:
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	_write_complete_xray_pdb(pdb_dir / "XRAY.pdb")
	out_dir = tmp_path / "outputs"
	params = {
		"method": "xray",
		"output_mode": "parser_outputs",
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone",
		"distance_constraints": "precise",
		"max_distance": "20.0",
		"vdw_constraints": "no",
		"torsion_angle_width": "0.0",
		"percentage_backbone_torsion_angles": "100.0",
	}

	parser(params, pdb_dir, tmp_path / "seeds", out_dir, "XRAY", True)

	assert run_remap_job(
		{"model_number": "1", "chain_id": "A", "order_id": "1"},
		out_dir,
		"XRAY",
		1,
	) == 0


@pytest.mark.parametrize("tightening", ["yes", "no"])
def test_real_parser_tightening_follows_its_own_option(tmp_path: Path, tightening: str) -> None:
	"""VDW bounds are generated with vdw_constraints; their 999 upper bound is tightened only
	when tightening_upper_bounds is yes."""
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	_write_complete_nmr_pdb(pdb_dir / "TEST.pdb")
	out_dir = tmp_path / "outputs"
	params = {
		"method": "nmr",
		"output_mode": "parser_outputs",
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone_plus_hydrogens",
		"distance_constraints": "interval_centered",
		"epsilon_short": "1.0",
		"epsilon_long": "2.0",
		"max_distance": "20.0",
		"vdwr_hh": "1.8",
		"vdw_constraints": "yes",
		"tightening_upper_bounds": tightening,
		"torsion_angle_width": "0.0",
		"percentage_backbone_torsion_angles": "100.0",
	}

	parser(params, pdb_dir, tmp_path / "seeds", out_dir, "TEST", True)
	instance = read_space_separated_file(parser_distance_file(out_dir, "TEST", 1, "A"))

	assert not instance.empty
	assert (instance.iloc[:, 4] >= 0.0).all()
	if tightening == "yes":
		assert (instance.iloc[:, 5] < 999.0).all()
	else:
		assert (instance.iloc[:, 5] >= 999.0).any()


def test_filtered_atoms_run_removes_stale_parser_and_remap_artifacts(
	tmp_path: Path,
) -> None:
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	_write_complete_nmr_pdb(pdb_dir / "TEST.pdb")
	out_dir = tmp_path / "outputs"
	params = {
		"method": "nmr", "output_mode": "parser_outputs", "model_number": "1", "chain_id": "A",
		"atom_selection": "backbone", "distance_constraints": "precise", "max_distance": "5.0",
		"vdw_constraints": "no", "torsion_angle_width": "0.0", "percentage_backbone_torsion_angles": "100.0",
	}
	parser(params, pdb_dir, tmp_path / "seeds", out_dir, "TEST", True)
	remapped_dir = out_dir / "TEST" / "remapped"
	remapped_dir.mkdir()
	(remapped_dir / "stale.dat").write_text("old", encoding="utf-8")

	filtered_params = {**params, "output_mode": "filtered_atoms"}
	parser(filtered_params, pdb_dir, tmp_path / "seeds", out_dir, "TEST", True)
	assert parser_structure_file(out_dir, "TEST", 1, "A").exists() is False
	assert parser_angular_file(out_dir, "TEST", 1, "A").exists() is False
	assert parser_distance_file(out_dir, "TEST", 1, "A").exists() is False
	assert not remapped_dir.exists()


def test_filtered_atoms_cleanup_failure_restores_complete_previous_generation(
	tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
	output_dir = tmp_path / "outputs" / "TEST"
	output_dir.mkdir(parents=True)
	remapped_dir = output_dir / "remapped"
	remapped_dir.mkdir()
	old_paths = {
		output_dir / name: f"old-{name}"
		for name in ("X.dat", "A.dat", "I.dat")
	}
	for path, content in old_paths.items():
		path.write_text(content, encoding="utf-8")
	(remapped_dir / "old.dat").write_text("old-remap", encoding="utf-8")

	original_rmtree = artifacts_pipeline.shutil.rmtree
	state = {"calls": 0}

	def fail_second_rmtree(path: object, *args: object, **kwargs: object) -> None:
		state["calls"] += 1
		if state["calls"] == 1:
			raise OSError("simulated cleanup failure")
		original_rmtree(path, *args, **kwargs)

	monkeypatch.setattr(artifacts_pipeline.shutil, "rmtree", fail_second_rmtree)
	context = {
		"output_files": {
			"filtered_atoms": output_dir / "X_filtered.dat",
			"metadata": output_dir / "metadata.json",
			"x": output_dir / "X.dat",
			"a": output_dir / "A.dat",
			"i": output_dir / "I.dat",
			"tsv_structure": tmp_path / "tmp" / "filtered_atoms.dat",
		},
		"output_mode": "filtered_atoms",
		"method": "nmr",
		"atom_selection": "backbone",
		"protein_chain": None,
		"tsv_structure_file": tmp_path / "tmp" / "filtered_atoms.dat",
		"tmp_dir": tmp_path / "tmp",
	}
	context["tmp_dir"].mkdir()
	monkeypatch.setattr(artifacts_pipeline, "_prepare_filtered_structure_file", lambda *args: None)
	monkeypatch.setattr(artifacts_pipeline, "_write_parser_metadata", lambda path, *args: None)

	with pytest.raises(OSError, match="simulated cleanup failure"):
		artifacts_pipeline.write_parser_artifacts(context)  # type: ignore[arg-type]

	assert {path: path.read_text(encoding="utf-8") for path in old_paths} == old_paths
	assert (remapped_dir / "old.dat").read_text(encoding="utf-8") == "old-remap"


def test_filtered_atoms_staging_failure_restores_already_moved_artifacts(
	tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
	output_dir = tmp_path / "outputs" / "TEST"
	output_dir.mkdir(parents=True)
	remapped_dir = output_dir / "remapped"
	remapped_dir.mkdir()
	old_paths = {
		output_dir / name: f"old-{name}"
		for name in ("X.dat", "A.dat", "I.dat")
	}
	for path, content in old_paths.items():
		path.write_text(content, encoding="utf-8")
	(remapped_dir / "old.dat").write_text("old-remap", encoding="utf-8")

	original_copy2 = artifacts_pipeline.shutil.copy2
	state = {"calls": 0}

	def fail_second_copy(source: object, destination: object, *args: object, **kwargs: object) -> object:
		state["calls"] += 1
		if state["calls"] == 2:
			raise OSError("simulated staging failure")
		return original_copy2(source, destination, *args, **kwargs)

	monkeypatch.setattr(artifacts_pipeline.shutil, "copy2", fail_second_copy)
	context = {
		"output_files": {
			"filtered_atoms": output_dir / "X_filtered.dat",
			"metadata": output_dir / "metadata.json",
			"x": output_dir / "X.dat",
			"a": output_dir / "A.dat",
			"i": output_dir / "I.dat",
			"tsv_structure": tmp_path / "tmp" / "filtered_atoms.dat",
		},
		"output_mode": "filtered_atoms",
		"method": "nmr",
		"atom_selection": "backbone",
		"protein_chain": None,
		"tsv_structure_file": tmp_path / "tmp" / "filtered_atoms.dat",
		"tmp_dir": tmp_path / "tmp",
	}
	context["tmp_dir"].mkdir()
	monkeypatch.setattr(artifacts_pipeline, "_prepare_filtered_structure_file", lambda *args: None)
	monkeypatch.setattr(artifacts_pipeline, "_write_parser_metadata", lambda path, *args: None)

	with pytest.raises(OSError, match="simulated staging failure"):
		artifacts_pipeline.write_parser_artifacts(context)  # type: ignore[arg-type]

	assert {path: path.read_text(encoding="utf-8") for path in old_paths} == old_paths
	assert (remapped_dir / "old.dat").read_text(encoding="utf-8") == "old-remap"


def test_parser_outputs_remove_persistent_filtered_atoms_artifact(
	tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
	output_dir = tmp_path / "outputs" / "TEST"
	tmp_dir = output_dir / "tmp"
	tmp_dir.mkdir(parents=True)
	filtered = output_dir / "X_filtered.dat"
	filtered.write_text("stale-filtered", encoding="utf-8")
	output_files = {
		"filtered_atoms": filtered,
		"metadata": output_dir / "metadata.json",
		"x": output_dir / "X.dat",
		"a": output_dir / "A.dat",
		"i": output_dir / "I.dat",
		"a_tmp": tmp_dir / "A.tmp",
		"dc_1": tmp_dir / "dc1",
		"dc_2": tmp_dir / "dc2",
		"dc_3": tmp_dir / "dc3",
		"dc_4": tmp_dir / "dc4",
		"dc_5": tmp_dir / "dc5",
		"tsv_structure": tmp_dir / "filtered_atoms.dat",
	}
	monkeypatch.setattr(artifacts_pipeline, "_prepare_filtered_structure_file", lambda *args: None)
	monkeypatch.setattr(artifacts_pipeline, "build_covalent_and_planar_topology", lambda path: {})
	monkeypatch.setattr(artifacts_pipeline, "_write_parser_outputs", lambda context, topology: None)
	context = {
		"output_files": output_files, "output_mode": "parser_outputs", "method": "nmr",
		"atom_selection": "backbone", "protein_chain": None,
		"tsv_structure_file": output_files["tsv_structure"], "tmp_dir": tmp_dir,
	}
	artifacts_pipeline.write_parser_artifacts(context)  # type: ignore[arg-type]
	assert not filtered.exists()


def test_remap_requires_parser_outputs_metadata_when_sidecar_exists(tmp_path: Path) -> None:
	out_dir = tmp_path / "outputs"
	metadata = parser_metadata_file(out_dir, "TEST", 1, "A")
	metadata.parent.mkdir(parents=True)
	metadata.write_text(
		json.dumps({"method": "nmr", "atom_selection": "backbone", "output_mode": "filtered_atoms"}),
		encoding="utf-8",
	)
	with pytest.raises(ValueError, match="output_mode: parser_outputs"):
		run_remap_job({"model_number": "1", "chain_id": "A", "order_id": "1"}, out_dir, "TEST", 1)


@pytest.mark.parametrize(
	("method", "atom_selection"),
	[
		("nmr", "full_chain"),
		("nmr", "backbone_plus_neighbors"),
		("xray", "backbone_plus_hydrogens"),
		("xray", "backbone_plus_neighbors"),
	],
)
def test_remap_rejects_public_layouts_outside_method_policy(
	tmp_path: Path,
	method: str,
	atom_selection: str,
) -> None:
	metadata = parser_metadata_file(tmp_path / "outputs", "TEST", 1, "A")
	metadata.parent.mkdir(parents=True)
	metadata.write_text(
		json.dumps(
			{"method": method, "atom_selection": atom_selection, "output_mode": "parser_outputs"}
		),
		encoding="utf-8",
	)

	with pytest.raises(ValueError, match="Remapping is not defined"):
		run_remap_job(
			{"model_number": "1", "chain_id": "A", "order_id": "1"},
			tmp_path / "outputs",
			"TEST",
			1,
		)


def test_remap_output_transaction_restores_existing_files_on_stage_failure(
	tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
	"""A failed T write must not leave a mixed-generation X/I/T set."""
	output_dir = tmp_path / "outputs" / "TEST"
	remapped_dir = output_dir / "remapped"
	remapped_dir.mkdir(parents=True)
	final_paths = (
		remapped_dir / "X_TEST_model1_chainA_ddgp_hc_order_5.dat",
		remapped_dir / "I_TEST_model1_chainA_ddgp_hc_order_5.dat",
		remapped_dir / "T_TEST_model1_chainA_ddgp_hc_order_5.dat",
	)
	for path in final_paths:
		path.write_text(f"old-{path.name}", encoding="utf-8")

	def fail_t_write(*args: object, **kwargs: object) -> None:
		raise OSError("simulated T write failure")

	monkeypatch.setattr(remapping_pipeline, "save_cliques_from_matrix_T", fail_t_write)
	with pytest.raises(OSError, match="simulated T write failure"):
		remapping_pipeline._save_remapped_outputs(
			pd.DataFrame([[1, 2, 1, 1, 1.0, 1.0, "N", "CA", "ALA", "ALA"]]),
			pd.DataFrame([[1, "N", 1, "ALA", 0.0, 0.0, 0.0]]),
			np.zeros((1, 7)),
			output_dir,
			"TEST",
			1,
			"A",
			5,
			"backbone_plus_hydrogens",
		)

	assert [path.read_text(encoding="utf-8") for path in final_paths] == [
		f"old-{path.name}" for path in final_paths
	]
	assert not list(remapped_dir.glob(".remap-staging-*"))


def test_parser_output_transaction_restores_existing_generation_on_failure(
	tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
	"""A failed parser stage cannot mix new metadata with old parser outputs."""
	output_dir = tmp_path / "outputs"
	tmp_dir = output_dir / "TEST" / "tmp"
	final_dir = output_dir / "TEST"
	final_dir.mkdir(parents=True)
	tmp_dir.mkdir()
	output_files = {
		"filtered_atoms": final_dir / "X_filtered.dat",
		"metadata": final_dir / "metadata.json",
		"x": final_dir / "X.dat",
		"a": final_dir / "A.dat",
		"i": final_dir / "I.dat",
		"a_tmp": tmp_dir / "A.tmp",
		"dc_1": tmp_dir / "dc1",
		"dc_2": tmp_dir / "dc2",
		"dc_3": tmp_dir / "dc3",
		"dc_4": tmp_dir / "dc4",
		"dc_5": tmp_dir / "dc5",
		"tsv_structure": tmp_dir / "filtered_atoms.dat",
	}
	for path in output_files.values():
		if path.parent == final_dir:
			path.write_text("old-generation", encoding="utf-8")

	def prepare(*args: object, **kwargs: object) -> None:
		output_files["filtered_atoms"].write_text("new-generation", encoding="utf-8")

	def fail_outputs(*args: object, **kwargs: object) -> None:
		raise OSError("simulated parser failure")

	monkeypatch.setattr(artifacts_pipeline, "_prepare_filtered_structure_file", prepare)
	monkeypatch.setattr(artifacts_pipeline, "build_covalent_and_planar_topology", lambda path: {})
	monkeypatch.setattr(artifacts_pipeline, "_write_parser_outputs", fail_outputs)
	context = {
		"output_files": output_files,
		"output_mode": "parser_outputs",
		"method": "nmr",
		"atom_selection": "backbone",
		"protein_chain": None,
		"tsv_structure_file": output_files["tsv_structure"],
		"tmp_dir": tmp_dir,
	}

	with pytest.raises(OSError, match="simulated parser failure"):
		artifacts_pipeline.write_parser_artifacts(context)  # type: ignore[arg-type]

	assert all(
		path.read_text(encoding="utf-8") == "old-generation"
		for path in (
			output_files["filtered_atoms"],
			output_files["metadata"],
			output_files["x"],
			output_files["a"],
			output_files["i"],
		)
	)
	assert not tmp_dir.exists()


def test_parser_rollback_reports_cleanup_failure_without_masking_original(
	tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
	output_dir = tmp_path / "outputs" / "TEST"
	tmp_dir = output_dir / "tmp"
	tmp_dir.mkdir(parents=True)
	output_files = {
		"filtered_atoms": output_dir / "X_filtered.dat",
		"metadata": output_dir / "metadata.json",
		"x": output_dir / "X.dat",
		"a": output_dir / "A.dat",
		"i": output_dir / "I.dat",
		"tsv_structure": tmp_dir / "filtered_atoms.dat",
	}
	for path in output_files.values():
		if path.parent == output_dir:
			path.write_text("old", encoding="utf-8")

	monkeypatch.setattr(
		artifacts_pipeline,
		"_prepare_filtered_structure_file",
		lambda *args: (_ for _ in ()).throw(OSError("original parser failure")),
	)
	monkeypatch.setattr(
		Path,
		"write_bytes",
		lambda self, data: (_ for _ in ()).throw(OSError("restore failed")),
	)
	context = {
		"output_files": output_files,
		"output_mode": "parser_outputs",
		"method": "nmr",
		"atom_selection": "backbone",
		"protein_chain": None,
		"tsv_structure_file": output_files["tsv_structure"],
		"tmp_dir": tmp_dir,
	}

	with pytest.raises(OSError, match="original parser failure") as failure:
		artifacts_pipeline.write_parser_artifacts(context)  # type: ignore[arg-type]

	notes = getattr(failure.value, "__notes__", ())
	assert any("Rollback cleanup failed" in note for note in notes) or "Rollback cleanup failed" in capsys.readouterr().err


def test_reexecuting_parser_preserves_previous_generation_on_failure(
	tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
	"""A failed rerun leaves every previously published parser file unchanged."""
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	_write_complete_nmr_pdb(pdb_dir / "TEST.pdb")
	out_dir = tmp_path / "outputs"
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

	parser(params, pdb_dir, tmp_path / "seeds", out_dir, "TEST", True)
	published = {
		path: path.read_bytes()
		for path in (out_dir / "TEST").glob("*")
		if path.is_file()
	}

	def fail_outputs(*args: object, **kwargs: object) -> None:
		raise OSError("simulated rerun failure")

	monkeypatch.setattr(artifacts_pipeline, "_write_parser_outputs", fail_outputs)
	with pytest.raises(OSError, match="simulated rerun failure"):
		parser(params, pdb_dir, tmp_path / "seeds", out_dir, "TEST", True)

	assert {path: path.read_bytes() for path in published} == published


@pytest.mark.parametrize(("output_mode", "warned"), [("parser_outputs", True), ("filtered_atoms", False)])
def test_xray_constraint_conversion_warnings_only_when_constraints_are_written(
	tmp_path: Path, capsys: pytest.CaptureFixture[str], output_mode: str, warned: bool
) -> None:
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	_write_complete_xray_pdb(pdb_dir / "XRAY.pdb")
	params = {
		"method": "auto",
		"output_mode": output_mode,
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone_plus_hydrogens",
		"distance_constraints": "interval_experimental",
		"vdw_constraints": "no",
		"torsion_angle_width": "60.0",
		"percentage_backbone_torsion_angles": "100.0",
		"percentage_hh_distance_constraints": "50",
	}

	parser(params, pdb_dir, tmp_path / "seeds", tmp_path / "outputs", "XRAY", True)
	err = capsys.readouterr().err
	# The selection conversion changes the filtered atoms, so it is always reported.
	assert "unavailable for X-ray" in err
	assert ("is NMR-specific; using 'interval_centered'" in err) is warned
	assert ("ignored for X-ray CA--CA constraints" in err) is warned

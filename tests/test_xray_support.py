from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pdb_parser.io import detect_experimental_method, xray_ca_distance_constraints
from pdb_parser.io.validate import ensure_model_chain_ready


def _write_atoms(path: Path) -> None:
	pd.DataFrame(
		[
			(1, "CA", 1, "ALA", 0.0, 0.0, 0.0),
			(2, "CA", 2, "ALA", 3.8, 0.0, 0.0),
			(3, "CA", 3, "ALA", 6.5, 0.0, 0.0),
			(4, "CA", 4, "ALA", 10.0, 0.0, 0.0),
		],
		columns=["atom_id", "atom_name", "resid", "resname", "x", "y", "z"],
	).to_csv(path, sep="\t", index=False)


def test_xray_detection_and_single_model_fallback(tmp_path: Path, capsys) -> None:
	pdb_path = tmp_path / "minimal_xray.pdb"
	pdb_path.write_text(
		"EXPDTA    X-RAY DIFFRACTION\n"
		"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N  \n"
		"ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C  \n"
		"ATOM      3  C   ALA A   1       1.000   1.000   0.000  1.00 20.00           C  \n"
		"TER\nEND\n",
		encoding="utf-8",
	)

	assert detect_experimental_method(pdb_path) == "xray"
	atoms, effective_model, method = ensure_model_chain_ready(
		pdb_path,
		2,
		"A",
		method="xray",
	)

	assert atoms.n_atoms > 0
	assert effective_model == 1
	assert method == "xray"
	assert "replaced by model 1" in capsys.readouterr().err


def test_xray_ca_distances_skip_adjacent_pairs_and_use_precise_mode(tmp_path: Path) -> None:
	atoms = tmp_path / "atoms.tsv"
	output = tmp_path / "xray.tsv"
	_write_atoms(atoms)

	xray_ca_distance_constraints(atoms, output, "precise", None, 20.0)
	result = pd.read_csv(output, sep="\t")

	assert set(zip(result.resid_i, result.resid_j)) == {(3, 1), (4, 1), (4, 2)}
	assert np.allclose(result.d_l, result.d_u)


def test_xray_ca_centered_intervals_use_long_epsilon_and_vdw_lower_bound(tmp_path: Path) -> None:
	atoms = tmp_path / "atoms.tsv"
	output = tmp_path / "xray.tsv"
	_write_atoms(atoms)

	xray_ca_distance_constraints(
		atoms,
		output,
		"interval_centered",
		2.0,
		20.0,
		np.random.default_rng(123),
	)
	result = pd.read_csv(output, sep="\t")

	assert (result.d_l < result.d_u).all()
	assert (result.d_u <= 20.0).all()
	assert (result.d_l < np.array([6.5, 10.0, 6.2])).all()


def test_vectorized_xray_ca_sampling_consumes_the_random_stream_like_the_scalar_loop() -> None:
	from pdb_parser.io.distance_intervals import (
		get_xray_ca_centered_interval,
		sample_centered_intervals,
	)
	from pdb_parser.io.vdw import VDW_CARBON_RADIUS, lower_bounds_below

	count = 200_000
	dist = np.random.default_rng(5).uniform(3.5, 4.95, count)
	scalar_rng = np.random.default_rng(31)
	expected, rejected_pairs = [], 0
	for value in dist.tolist():
		before = scalar_rng.bit_generator.state
		expected.append(get_xray_ca_centered_interval(value, 2.0, 5.0, scalar_rng))
		probe = np.random.default_rng()
		probe.bit_generator.state = before
		probe.normal()
		rejected_pairs += probe.bit_generator.state != scalar_rng.bit_generator.state

	vector_rng = np.random.default_rng(31)
	lower = lower_bounds_below(dist, np.full(count, 2.0 * VDW_CARBON_RADIUS))
	d_l, d_u, accepted = sample_centered_intervals(dist, np.full(count, 2.0), lower, 5.0, vector_rng)

	assert rejected_pairs > 0
	assert accepted.all() and all(interval is not None for interval in expected)
	assert np.array_equal(d_l, np.array([interval[0] for interval in expected]))
	assert np.array_equal(d_u, np.array([interval[1] for interval in expected]))
	assert vector_rng.bit_generator.state == scalar_rng.bit_generator.state


def test_xray_ca_constraints_skip_same_and_adjacent_residues_in_table_order(tmp_path: Path) -> None:
	atoms = tmp_path / "atoms.tsv"
	_write_atoms(atoms)
	output = tmp_path / "out.dat"
	table = xray_ca_distance_constraints(atoms, output, "precise", None, 20.0)

	# Residues 1-4: only the non-adjacent pairs (3,1), (4,1), (4,2), later atom first.
	assert list(zip(table["atom_id_i"], table["atom_id_j"])) == [(3, 1), (4, 1), (4, 2)]
	assert output.read_text(encoding="utf-8").splitlines()[1].split("\t")[:6] == [
		"3", "1", "3", "1", "6.500000000000", "6.500000000000",
	]

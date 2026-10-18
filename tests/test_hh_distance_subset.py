"""percentage_hh_distance_constraints: random subset of the NMR long-range restraints."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from test_pipeline_integration import _write_complete_nmr_pdb

from pdb_parser.io.distance_constraints import nmr_distance_constraints
from pdb_parser.pipeline.parser_pipeline import parser
from pdb_parser.pipeline.validation import (
	_preflight_parser_params,
	required_seed_fields_for_execution,
)
from pdb_parser.utils.inputs import read_params

N_ATOMS = 6
N_PAIRS = N_ATOMS * (N_ATOMS - 1) // 2


def _write_atoms(path: Path, atom_names: tuple[str, ...]) -> Path:
	rows = ["atom_id\tatom_name\tresid\tresname\tx\ty\tz"]
	atom_id = 1
	for resid in range(1, N_ATOMS + 1):
		for atom_name in atom_names:
			rows.append(f"{atom_id}\t{atom_name}\t{resid}\tALA\t{0.5 * atom_id:.3f}\t0.0\t0.0")
			atom_id += 1
	path.write_text("\n".join(rows) + "\n", encoding="utf-8")
	return path


def _run(
	tmp_path: Path,
	mode: str,
	percentage: float,
	*,
	subset_seed: int | None = 7,
	atom_selection: str = "backbone_plus_hydrogens",
	atom_names: tuple[str, ...] = ("H",),
	name: str = "out",
) -> list[str]:
	input_file = _write_atoms(tmp_path / "filtered_atoms.dat", atom_names)
	output_file = tmp_path / f"{name}.dat"
	nmr_distance_constraints(
		input_file, output_file, mode, atom_selection,
		1.0, 2.0, 20.0, 2.5, 3.5, 5.0, 0.4,
		np.random.default_rng(11),
		hh_distance_percentage=percentage,
		rng_hh_distance_subset=None if subset_seed is None else np.random.default_rng(subset_seed),
	)
	return output_file.read_text(encoding="utf-8").splitlines()[1:]


@pytest.mark.parametrize("mode", ["precise", "interval_centered", "interval_experimental"])
def test_subset_keeps_the_requested_share_of_the_full_restraint_list(tmp_path: Path, mode: str) -> None:
	full = _run(tmp_path, mode, 100.0, name="full")
	subset = _run(tmp_path, mode, 40.0, name="subset")

	assert len(full) == N_PAIRS
	assert len(subset) == round(N_PAIRS * 0.4)
	# Same rows (including the centred intervals) as the 100% run, in the same order.
	assert subset == [row for row in full if row in subset]


@pytest.mark.parametrize("mode", ["precise", "interval_centered", "interval_experimental"])
def test_subset_is_reproducible_for_a_seed_and_varies_between_seeds(tmp_path: Path, mode: str) -> None:
	first = _run(tmp_path, mode, 50.0, subset_seed=3, name="a")
	again = _run(tmp_path, mode, 50.0, subset_seed=3, name="b")
	other = _run(tmp_path, mode, 50.0, subset_seed=4, name="c")

	assert first == again
	assert first != other


def test_zero_and_full_percentages_need_no_random_generator(tmp_path: Path) -> None:
	assert _run(tmp_path, "precise", 0.0, subset_seed=None, name="none") == []
	assert len(_run(tmp_path, "precise", 100.0, subset_seed=None, name="all")) == N_PAIRS


def test_partial_percentage_requires_a_random_generator(tmp_path: Path) -> None:
	with pytest.raises(ValueError, match="random generator"):
		_run(tmp_path, "precise", 50.0, subset_seed=None)


def test_backbone_proxy_restraints_are_subsampled_too(tmp_path: Path) -> None:
	full = _run(tmp_path, "precise", 100.0, atom_selection="backbone", atom_names=("N", "CA", "C"), name="full")
	subset = _run(tmp_path, "precise", 50.0, atom_selection="backbone", atom_names=("N", "CA", "C"), name="half")

	assert len(subset) == round(len(full) * 0.5)
	assert set(subset) <= set(full)


def _params(**overrides: str) -> dict[str, str]:
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
		"vdw_constraints": "no",
		"torsion_angle_width": "0.0",
		"percentage_backbone_torsion_angles": "100.0",
	}
	params.update(overrides)
	return params


def test_percentage_is_optional_and_defaults_to_all_restraints() -> None:
	assert _preflight_parser_params(_params())["percentage_hh_distance_constraints"] == 100.0


def test_precise_mode_keeps_the_requested_percentage() -> None:
	parsed = _preflight_parser_params(
		_params(distance_constraints="precise", percentage_hh_distance_constraints="30")
	)
	assert parsed["percentage_hh_distance_constraints"] == 30.0


@pytest.mark.parametrize("value", ["-1", "100.5", "nan", "abc"])
def test_percentage_outside_zero_to_hundred_is_rejected(value: str) -> None:
	with pytest.raises(ValueError, match="percentage_hh_distance_constraints"):
		_preflight_parser_params(_params(percentage_hh_distance_constraints=value))


@pytest.mark.parametrize(
	("method", "percentage", "needs_seed"),
	[("nmr", 50.0, True), ("nmr", 100.0, False), ("nmr", 0.0, False), ("xray", 50.0, False)],
)
def test_subset_seed_is_required_only_for_a_random_nmr_subset(
	method: str, percentage: float, needs_seed: bool
) -> None:
	fields = required_seed_fields_for_execution(
		"precise", 0.0, 100.0, 3, method=method, hh_distance_percentage=percentage
	)
	assert ("seed_hh_distance_subset" in fields) is needs_seed


def test_parser_subsamples_nmr_restraints_reproducibly(tmp_path: Path) -> None:
	pdb_dir = tmp_path / "pdb"
	pdb_dir.mkdir()
	_write_complete_nmr_pdb(pdb_dir / "TEST.pdb")
	seed_dir = tmp_path / "seeds"

	def noe_rows(out_name: str, percentage: str) -> list[str]:
		out_dir = tmp_path / out_name
		parser(_params(percentage_hh_distance_constraints=percentage), pdb_dir, seed_dir, out_dir, "TEST", False)
		lines = (out_dir / "TEST" / "tmp" / "distance_constraints_4.dat").read_text(encoding="utf-8")
		return lines.splitlines()[1:]

	full = noe_rows("full", "100")
	half = noe_rows("half", "50")
	half_again = noe_rows("half_again", "50")
	edges = [row for row in full if _is_discretization_edge(row)]

	# H_i--HA_i for the 3 residues and HA_{i-1}--H_i for residues 2 and 3.
	assert len(edges) == 5
	assert len(half) == len(edges) + round((len(full) - len(edges)) * 0.5)
	assert set(edges) <= set(half) <= set(full)
	assert half == half_again
	assert "seed_hh_distance_subset" in read_params(seed_dir / "TEST" / "seed_TEST.dat")


def _is_discretization_edge(row: str) -> bool:
	fields = row.split()
	(resid_i, name_i), (resid_j, name_j) = sorted(((int(fields[2]), fields[6]), (int(fields[3]), fields[7])))
	if resid_i == resid_j:
		return {name_i, name_j} == {"H", "HA"}
	return resid_j == resid_i + 1 and (name_i, name_j) == ("HA", "H")


def _write_named_atoms(path: Path, residues: list[tuple[str, tuple[str, ...]]]) -> Path:
	rows = ["atom_id\tatom_name\tresid\tresname\tx\ty\tz"]
	atom_id = 1
	for resid, (resname, atom_names) in enumerate(residues, start=1):
		for atom_name in atom_names:
			rows.append(f"{atom_id}\t{atom_name}\t{resid}\t{resname}\t{0.5 * atom_id:.3f}\t0.0\t0.0")
			atom_id += 1
	path.write_text("\n".join(rows) + "\n", encoding="utf-8")
	return path


def _name_pairs(rows: list[str]) -> set[tuple[str, str]]:
	pairs = set()
	for row in rows:
		fields = row.split()
		pairs.add(tuple(sorted((f"{fields[2]}{fields[6]}", f"{fields[3]}{fields[7]}"))))
	return pairs


@pytest.mark.parametrize("percentage", [0.0, 30.0])
def test_discretization_edges_are_always_kept(tmp_path: Path, percentage: float) -> None:
	# N-terminal H1/H2/H3, a glycine (HA2/HA3) and a proline (HD2/HD3 as HN).
	input_file = _write_named_atoms(tmp_path / "filtered_atoms.dat", [
		("ALA", ("H1", "H2", "H3", "HA")),
		("GLY", ("H", "HA2", "HA3")),
		("PRO", ("HD2", "HD3", "HA")),
		("SER", ("H", "HA", "HB2")),
	])

	def run(share: float) -> list[str]:
		output_file = tmp_path / f"out_{share:g}.dat"
		nmr_distance_constraints(
			input_file, output_file, "precise", "backbone_plus_hydrogens",
			1.0, 2.0, 100.0, 2.5, 3.5, 5.0, 0.4,
			hh_distance_percentage=share, rng_hh_distance_subset=np.random.default_rng(5),
		)
		return output_file.read_text(encoding="utf-8").splitlines()[1:]

	full = run(100.0)
	rows = run(percentage)

	# HN and HA are resolved as in the remap stage: H1 > H > HD2 > HD3 and HA > HA2 > HA3.
	edges = {
		("1H1", "1HA"), ("1HA", "2H"), ("2H", "2HA2"), ("2HA2", "3HD2"),
		("3HA", "3HD2"), ("3HA", "4H"), ("4H", "4HA"),
	}
	assert edges <= _name_pairs(rows)
	assert set(rows) <= set(full)
	assert len(rows) == len(edges) + round((len(full) - len(edges)) * percentage / 100.0)

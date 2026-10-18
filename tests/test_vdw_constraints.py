from __future__ import annotations

from pathlib import Path

import pytest

from pdb_parser.io.vdw_constraints import (
	build_vdw_torsion_pairs,
	is_backbone_dihedral_path,
	vdw_distances,
)


def _topology() -> dict[str, object]:
	graph = {
		(1, "C"): {(1, "N"), (2, "N")},
		(1, "N"): {(1, "C")},
		(2, "N"): {(1, "C"), (2, "CA")},
		(2, "CA"): {(2, "N"), (2, "C")},
		(2, "C"): {(2, "CA")},
	}
	return {
		"atom_index": {(1, "C"): 0, (1, "N"): 1, (2, "N"): 2, (2, "CA"): 3, (2, "C"): 4},
		"logical_graph": graph,
		"bonded_pairs": {(1, 0), (2, 0), (3, 2), (4, 3)},
		"angle_pairs": {(2, 1), (3, 0), (4, 2)},
		"excluded_pairs": {(1, 0), (2, 0), (3, 2), (4, 3), (2, 1), (3, 0), (4, 2)},
	}


def test_backbone_dihedral_paths_are_recognized_in_both_directions() -> None:
	path = ((1, "C"), (2, "N"), (2, "CA"), (2, "C"))
	assert is_backbone_dihedral_path(path)
	assert is_backbone_dihedral_path(tuple(reversed(path)))
	assert not is_backbone_dihedral_path(((1, "N"), (1, "CA"), (2, "N"), (2, "CA")))


def test_vdw_torsion_pairs_separate_backbone_dihedrals() -> None:
	relaxed, backbone = build_vdw_torsion_pairs(_topology())  # type: ignore[arg-type]
	assert (4, 0) in backbone
	assert (4, 0) not in relaxed


def test_vdw_distances_rejects_empty_input_without_a_residue_sequence(tmp_path: Path) -> None:
	input_file = tmp_path / "filtered_atoms.dat"
	input_file.write_text("atom_id\tatom_name\tresid\tresname\tx\ty\tz\n", encoding="utf-8")
	output_file = tmp_path / "vdw.dat"

	with pytest.raises(ValueError, match="valid residue sequence"):
		vdw_distances(input_file, _topology(), output_file)  # type: ignore[arg-type]


def test_vectorized_lower_bounds_match_the_scalar_rule() -> None:
	import numpy as np

	from pdb_parser.io.vdw import lower_bound_below, lower_bounds_below

	rng = np.random.default_rng(3)
	base = np.concatenate([rng.uniform(0.5, 3.5, 5000), [2.4, 2.4, 2.4, 2.4, 3.0]])
	# Below, equal to and above the base, every percentage step, and the 0.9*distance fallback.
	distance = np.concatenate([rng.uniform(0.0, 4.0, 5000), [2.5, 2.4, 0.24, 0.1, 0.3]])

	expected = np.array([lower_bound_below(d, b) for d, b in zip(distance.tolist(), base.tolist())])
	assert np.array_equal(lower_bounds_below(distance, base), expected)

"""Direct-import regression coverage for the specialized pipeline modules."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pdb_parser.io.files import validate_positional_schema
from pdb_parser.pipeline import parser as package_parser
from pdb_parser.pipeline.artifacts import write_parser_artifacts
from pdb_parser.pipeline.context import prepare_parser_execution_context
from pdb_parser.pipeline.parser_pipeline import parser
from pdb_parser.pipeline.remapping import write_remapped_artifacts
from pdb_parser.pipeline.validation import validate_public_parser_params
from pdb_parser.remapping.clique_patterns import (
	build_clique_pairs as specialized_build_clique_pairs,
)
from pdb_parser.remapping.cliques import build_backbone_torsion_matrix
from pdb_parser.remapping.dataframes import (
	remap_distance_dataframe as specialized_remap_distance_dataframe,
)
from pdb_parser.remapping.layout import detect_remap_atom_layout
from pdb_parser.remapping.ordering import get_backbone_residue_numeric_order
from pdb_parser.remapping.vertex_ordering import (
	build_backbone_torsion_matrix as legacy_build_backbone_torsion_matrix,
)


def test_specialized_modules_are_directly_importable_without_the_facades() -> None:
	assert prepare_parser_execution_context.__module__ == "pdb_parser.pipeline.context"
	assert write_parser_artifacts.__module__ == "pdb_parser.pipeline.artifacts"
	assert write_remapped_artifacts.__module__ == "pdb_parser.pipeline.remapping"
	assert validate_public_parser_params.__module__ == "pdb_parser.pipeline.validation"
	assert detect_remap_atom_layout.__module__ == "pdb_parser.remapping.layout"
	assert get_backbone_residue_numeric_order.__module__ == "pdb_parser.remapping.ordering"
	assert build_backbone_torsion_matrix.__module__ == "pdb_parser.remapping.torsions"


def test_legacy_facades_preserve_the_specialized_public_entrypoints() -> None:
	assert package_parser is parser
	assert legacy_build_backbone_torsion_matrix is build_backbone_torsion_matrix


def test_legacy_clique_facade_matches_specialized_behavior() -> None:
	from pdb_parser.remapping.cliques import build_clique_pairs, remap_distance_dataframe

	available_atoms = [(4, "C"), (5, "N"), (5, "CA"), (5, "HA2"), (5, "H")]
	assert build_clique_pairs(5, 5, available_atoms) == specialized_build_clique_pairs(
		5,
		5,
		available_atoms,
	)
	distances = pd.DataFrame([[2, 1, 2, 1, 1.0, 2.0, "CA", "N", "GLY", "ALA"]])
	assert remap_distance_dataframe(distances, [1, 2]).equals(
		specialized_remap_distance_dataframe(distances, [1, 2])
	)


@pytest.mark.parametrize("schema, columns", [("X", 6), ("I", 9), ("A", 7)])
def test_positional_schema_validation_reports_file_and_expected_format(
	schema: str, columns: int
) -> None:
	with pytest.raises(ValueError, match=rf"Invalid {schema} schema.*expected"):
		validate_positional_schema(pd.DataFrame([[0] * columns]), schema, "broken.dat")


def test_positional_schema_validation_identifies_invalid_column() -> None:
	with pytest.raises(ValueError, match="column 5 'x'"):
		validate_positional_schema(
			pd.DataFrame([[1, "N", 1, "ALA", "bad", 0.0, 0.0]]),
			"X",
			"broken-x.dat",
		)


@pytest.mark.parametrize(
	("schema", "row", "message"),
	[
		("X", [1.5, "CA", 1, "ALA", 0.0, 0.0, 0.0], "positive integers"),
		("I", [1, 2, 1, 2, 3.0, 2.0, "CA", "CA", "ALA", "ALA"], "lower distance"),
		("A", [1, "ALA", 0.0, float("inf"), 0.0, 1.0, 0.0, 1.0], "finite"),
	],
)
def test_positional_schema_validation_enforces_semantic_invariants(
	schema: str, row: list[object], message: str
) -> None:
	with pytest.raises(ValueError, match=message):
		validate_positional_schema(pd.DataFrame([row]), schema, "invalid.dat")


@pytest.mark.parametrize(
	"row",
	[
		[1, 2, 1, 1, -0.1, 1.0, "N", "CA", "ALA", "ALA"],
		[1, 1, 1, 1, 0.0, 1.0, "N", "N", "ALA", "ALA"],
	]
)
def test_positional_schema_validation_rejects_invalid_physical_distance_pairs(
	row: list[object],
) -> None:
	with pytest.raises(ValueError, match="(non-negative|distinct atoms)"):
		validate_positional_schema(pd.DataFrame([row]), "I", "invalid-I.dat")


def test_positional_schema_validation_rejects_inconsistent_residue_metadata() -> None:
	df = pd.DataFrame(
		[
			[1, "CA", 1, "ALA", 0.0, 0.0, 0.0],
			[2, "N", 1, "GLY", 1.0, 0.0, 0.0],
		]
	)

	with pytest.raises(ValueError, match="consistent residue/atom metadata"):
		validate_positional_schema(df, "X", "inconsistent.dat")


def test_positional_schema_validation_requires_text_columns() -> None:
	with pytest.raises(ValueError, match="must contain text"):
		validate_positional_schema(
			pd.DataFrame([[1, 7, 1, "ALA", 0.0, 0.0, 0.0]]), "X", "typed.dat"
		)


def test_legacy_facades_match_specialized_results_for_representative_helpers() -> None:
	from pdb_parser.io import distance_constraints as legacy_io
	from pdb_parser.io.distance_intervals import get_precise_interval
	from pdb_parser.remapping import vertex_ordering as legacy_vertex
	from pdb_parser.remapping.ordering import get_internal_residue_ordering_template
	from pdb_parser.remapping.torsions import get_interval_dl_du

	assert legacy_io._get_precise_interval(2.5, 5.0) == get_precise_interval(2.5, 5.0)
	assert legacy_vertex.get_internal_residue_ordering_template(3, "GLY") == (
		get_internal_residue_ordering_template(3, "GLY")
	)

	intervals = {(1, 2): (4.0, 6.0)}
	assert legacy_vertex.get_interval_dl_du(1, 2, intervals) == get_interval_dl_du(
		1,
		2,
		intervals,
	)


def test_legacy_distance_facade_preserves_interval_behavior() -> None:
	from pdb_parser.io import distance_constraints as legacy_io
	from pdb_parser.io.distance_intervals import (
		get_centered_interval,
		get_experimental_interval,
		get_xray_ca_centered_interval,
	)

	assert legacy_io._get_experimental_interval(4.0, 1.5, 3.0, 5.0, 1.0) == (
		get_experimental_interval(4.0, 1.5, 3.0, 5.0, 1.0)
	)
	assert legacy_io._get_centered_interval(
		4.0,
		1,
		2,
		0.5,
		1.0,
		1.0,
		10.0,
		np.random.default_rng(12),
	) == get_centered_interval(
		4.0,
		1,
		2,
		0.5,
		1.0,
		1.0,
		10.0,
		np.random.default_rng(12),
	)
	assert legacy_io._get_xray_ca_centered_interval(
		4.0,
		0.5,
		5.0,
		np.random.default_rng(12),
	) == get_xray_ca_centered_interval(
		4.0,
		0.5,
		5.0,
		np.random.default_rng(12),
	)


def test_io_facade_matches_specialized_angular_conversion_behavior() -> None:
	import pandas as pd
	from numpy.testing import assert_allclose

	from pdb_parser.io import distance_constraints as legacy_io
	from pdb_parser.io.angular_constraints import (
		build_backbone_residue_map,
		torsion_interval_distance_bounds,
	)

	atoms = pd.DataFrame(
		[
			(1, 1, "ALA", "N", 0.0, 0.0, 0.0),
			(2, 1, "ALA", "CA", 1.0, 0.0, 0.0),
			(3, 1, "ALA", "C", 2.0, 0.0, 0.0),
		],
		columns=["atom_id", "resid", "resname", "atom_name", "x", "y", "z"],
	)

	assert legacy_io.wrap_angle_deg(540.0) == -180.0
	legacy_map = legacy_io.build_backbone_residue_map(atoms)
	specialized_map = build_backbone_residue_map(atoms)
	assert legacy_map.keys() == specialized_map.keys()
	for resid in legacy_map:
		assert legacy_map[resid]["resname"] == specialized_map[resid]["resname"]
		for atom_name in legacy_map[resid]["atoms"]:
			assert legacy_map[resid]["atoms"][atom_name]["atom_id"] == (
				specialized_map[resid]["atoms"][atom_name]["atom_id"]
			)
			assert_allclose(
				legacy_map[resid]["atoms"][atom_name]["coord"],
				specialized_map[resid]["atoms"][atom_name]["coord"],
			)
	legacy_bounds = legacy_io.torsion_interval_distance_bounds(
		atoms.iloc[0][["x", "y", "z"]].to_numpy(),
		atoms.iloc[1][["x", "y", "z"]].to_numpy(),
		atoms.iloc[2][["x", "y", "z"]].to_numpy(),
		atoms.iloc[1][["x", "y", "z"]].to_numpy(),
		0.0,
		10.0,
	)
	specialized_bounds = torsion_interval_distance_bounds(
		atoms.iloc[0][["x", "y", "z"]].to_numpy(),
		atoms.iloc[1][["x", "y", "z"]].to_numpy(),
		atoms.iloc[2][["x", "y", "z"]].to_numpy(),
		atoms.iloc[1][["x", "y", "z"]].to_numpy(),
		0.0,
		10.0,
	)
	assert_allclose(legacy_bounds, specialized_bounds)


def test_pipeline_facade_remap_wrapper_preserves_specialized_result(monkeypatch) -> None:
	from pdb_parser.pipeline import parser_pipeline

	expected = object()

	def fake_write(*args, **kwargs):
		return expected

	monkeypatch.setattr(parser_pipeline, "write_remapped_artifacts", fake_write)
	result = parser_pipeline.remap_instance({}, "outputs", "TEST", [])

	assert result is expected


def test_pipeline_parser_facade_preserves_stage_order(monkeypatch) -> None:
	from pdb_parser.pipeline import parser_pipeline

	context = {"tmp_dir": object()}
	called: list[object] = []

	monkeypatch.setattr(
		parser_pipeline,
		"prepare_parser_execution_context",
		lambda *args: context,
	)
	monkeypatch.setattr(
		parser_pipeline,
		"write_parser_artifacts",
		lambda value: called.append(("write", value)),
	)
	monkeypatch.setattr(
		parser_pipeline,
		"cleanup_parser_artifacts",
		lambda value, remove: called.append(("cleanup", value, remove)),
	)

	parser_pipeline.parser({}, "pdb", "seed", "out", "TEST", remove_tmp_dir=True)

	assert called == [
		("write", context),
		("cleanup", context, True),
	]


def test_legacy_vertex_facade_preserves_dataframe_transformation() -> None:
	import pandas as pd

	from pdb_parser.remapping import vertex_ordering as legacy_vertex
	from pdb_parser.remapping.dataframes import remap_structure_dataframe

	structure = pd.DataFrame(
		[
			[10, 1, 1, "N", 0.0, 0.0, 0.0],
			[20, 1, 1, "CA", 1.0, 0.0, 0.0],
		],
	)

	assert legacy_vertex.remap_structure_dataframe(structure, [20, 10]).equals(
		remap_structure_dataframe(structure, [20, 10])
	)


def test_legacy_vertex_facade_preserves_all_split_public_helpers() -> None:
	from pdb_parser import remapping
	from pdb_parser.remapping import vertex_ordering

	for name in (
		"print_skip",
		"get_missing_backbone_atoms",
		"get_missing_backbone_plus_hydrogen_atoms",
		"ensure_backbone_plus_hydrogens_dataframe",
		"get_numeric_atom_order_from_named_order",
		"remap_distance_dataframe",
		"remap_structure_dataframe",
	):
		assert name in vertex_ordering.__all__
		assert hasattr(vertex_ordering, name)

	assert remapping.remap_distance_dataframe is vertex_ordering.remap_distance_dataframe
	assert remapping.remap_structure_dataframe is vertex_ordering.remap_structure_dataframe

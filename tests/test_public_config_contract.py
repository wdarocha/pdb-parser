from __future__ import annotations

from pathlib import Path

import pytest

from pdb_parser.public_config import (
	INTERNAL_ATOM_SELECTIONS,
	PUBLIC_INPUT_ATOM_SELECTIONS,
	PUBLIC_PARSER_PARAM_KEYS,
	PUBLIC_REMAP_PARAM_KEYS,
	SUPPORTED_ORDER_IDS,
	VALID_DISTANCE_MODES,
	VALID_METHODS,
	VALID_OUTPUT_MODES,
	VALID_VDW_OPTIONS,
	parse_public_order_choice,
	parse_public_order_id,
	validate_public_remap_params,
)


def test_public_atom_selection_exposes_all_parser_atom_selections() -> None:
	assert PUBLIC_INPUT_ATOM_SELECTIONS <= INTERNAL_ATOM_SELECTIONS
	assert PUBLIC_INPUT_ATOM_SELECTIONS == {
		"full_chain",
		"backbone",
		"backbone_plus_hydrogens",
		"backbone_plus_neighbors",
	}


def test_all_public_atom_selections_are_accepted_by_configuration() -> None:
	from pdb_parser.pipeline import validate_public_parser_params

	for atom_selection in sorted(PUBLIC_INPUT_ATOM_SELECTIONS):
		params = {
			"model_number": "1", "chain_id": "A", "atom_selection": atom_selection,
			"distance_constraints": "precise", "vdw_constraints": "no",
			"torsion_angle_width": "0", "percentage_backbone_torsion_angles": "100",
			"max_distance": "20",
		}
		validate_public_parser_params(params)


def test_public_method_and_output_modes_match_contract() -> None:
	assert VALID_METHODS == {"auto", "nmr", "xray"}
	assert VALID_OUTPUT_MODES == {"filtered_atoms", "parser_outputs"}


def test_public_distance_modes_match_current_release_contract() -> None:
	assert VALID_DISTANCE_MODES == {
		"precise",
		"interval_centered",
		"interval_experimental",
	}


def test_public_vdw_options_match_current_release_contract() -> None:
	assert VALID_VDW_OPTIONS == {"yes", "no"}


def test_public_parser_keys_include_current_release_surface() -> None:
	assert {
		"method",
		"output_mode",
		"model_number",
		"chain_id",
		"atom_selection",
		"distance_constraints",
		"vdw_constraints",
		"tightening_upper_bounds",
		"torsion_angle_width",
		"percentage_backbone_torsion_angles",
		"percentage_hh_distance_constraints",
	} <= PUBLIC_PARSER_PARAM_KEYS


def test_public_remap_keys_match_current_release_contract() -> None:
	assert PUBLIC_REMAP_PARAM_KEYS == {"model_number", "chain_id", "order_id", "order_type"}


def test_supported_order_ids_match_current_release_contract() -> None:
	assert SUPPORTED_ORDER_IDS == set(range(1, 11))


def test_parse_public_order_id_accepts_supported_values() -> None:
	assert parse_public_order_id("9") == 9


def test_parse_public_order_id_rejects_unsupported_values() -> None:
	with pytest.raises(ValueError, match="Invalid order_id"):
		parse_public_order_id("11")


def test_parse_public_order_choice_rejects_non_numeric_order() -> None:
	with pytest.raises(ValueError, match="order_id"):
		parse_public_order_choice("random")


def test_validate_public_remap_params_accepts_current_public_surface(tmp_path: Path) -> None:
	params = {
		"model_number": "1",
		"chain_id": "A",
		"order_id": "9",
	}

	validate_public_remap_params(params, tmp_path / "instance_remap.cfg")


def test_validate_public_remap_params_rejects_unknown_keys(tmp_path: Path) -> None:
	params = {
		"model_number": "1",
		"chain_id": "A",
		"order_id": "9",
		"secret_mode": "yes",
	}

	with pytest.raises(ValueError, match="Unsupported remapping parameter\\(s\\): secret_mode"):
		validate_public_remap_params(params, tmp_path / "instance_remap.cfg")

from __future__ import annotations

from pathlib import Path

import pytest

from pdb_parser.instance_remap import read_instance_remap_params
from pdb_parser.pipeline.parser_pipeline import _preflight_parser_params
from pdb_parser.public_config import (
	parse_public_order_choice,
	parse_public_order_id,
)
from pdb_parser.utils.fs import ensure_existing_dir
from pdb_parser.utils.inputs import (
	normalize_nonempty_config_value,
	parse_int_config_value,
	read_params,
	read_pdb_ids,
	reject_unknown_keys,
	require_keys,
)
from pdb_parser.utils.seeds import load_or_create_seed_bundle


def write_text_file(path: Path, content: str) -> Path:
	path.write_text(content, encoding="utf-8")
	return path


def test_read_pdb_ids_deduplicates_normalizes_and_ignores_comments(tmp_path: Path) -> None:
	pdb_ids_file = write_text_file(
		tmp_path / "pdb_ids.txt",
		"""
		# comment
		1abc
		2def   # inline comment
		1ABC
		
		3GhI
		""",
	)

	assert read_pdb_ids(pdb_ids_file) == ["1ABC", "2DEF", "3GHI"]


def test_read_pdb_ids_rejects_invalid_identifier(tmp_path: Path) -> None:
	pdb_ids_file = write_text_file(tmp_path / "pdb_ids.txt", "ABCDE\n")

	with pytest.raises(ValueError, match="Invalid PDB id"):
		read_pdb_ids(pdb_ids_file)


def test_read_params_parses_key_value_file_with_comments(tmp_path: Path) -> None:
	params_file = write_text_file(
		tmp_path / "params.cfg",
		"""
		# comment
		Model_Number: 1
		chain_id: A # inline comment
		vdw_constraints: yes
		""",
	)

	assert read_params(params_file) == {
		"model_number": "1",
		"chain_id": "A",
		"vdw_constraints": "yes",
	}


def test_read_params_rejects_duplicate_keys(tmp_path: Path) -> None:
	params_file = write_text_file(
		tmp_path / "params.cfg",
		"model_number: 1\nmodel_number: 2\n",
	)

	with pytest.raises(ValueError, match="Duplicate parameter 'model_number'"):
		read_params(params_file)


def test_reject_unknown_keys_reports_public_context() -> None:
	with pytest.raises(ValueError, match="Unsupported parser parameter\\(s\\): hidden"):
		reject_unknown_keys({"hidden": "1"}, {"model_number"}, "parser")


def test_require_keys_reports_missing_public_keys(tmp_path: Path) -> None:
	with pytest.raises(ValueError, match="Missing required remapping parameter\\(s\\)"):
		require_keys({"model_number": "1"}, ("model_number", "chain_id"), tmp_path / "cfg", "remapping")


def test_parse_int_config_value_accepts_valid_integer() -> None:
	assert parse_int_config_value("7", "model_number", "remapping") == 7


def test_parse_int_config_value_rejects_invalid_integer() -> None:
	with pytest.raises(ValueError, match="model_number must be an integer"):
		parse_int_config_value("abc", "model_number", "remapping")


def test_normalize_nonempty_config_value_strips_value() -> None:
	assert normalize_nonempty_config_value(" A ", "chain_id", "remapping") == "A"


def test_normalize_nonempty_config_value_rejects_blank_string() -> None:
	with pytest.raises(ValueError, match="chain_id cannot be empty"):
		normalize_nonempty_config_value("   ", "chain_id", "remapping")


def test_ensure_existing_dir_accepts_existing_directory(tmp_path: Path) -> None:
	assert ensure_existing_dir(tmp_path, "Output directory") == tmp_path


def test_ensure_existing_dir_rejects_missing_directory(tmp_path: Path) -> None:
	with pytest.raises(ValueError, match="Output directory does not exist"):
		ensure_existing_dir(tmp_path / "missing", "Output directory")


def test_parser_preflight_accepts_public_atom_selection() -> None:
	params = {
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone_plus_hydrogens",
		"distance_constraints": "precise",
		"vdw_constraints": "yes",
		"torsion_angle_width": "0",
		"percentage_backbone_torsion_angles": "100",
		"max_distance": "5",
	}

	parsed = _preflight_parser_params(params)

	assert parsed["atom_selection"] == "backbone_plus_hydrogens"


def test_parser_preflight_forces_vdw_constraints_to_no_for_precise_mode() -> None:
	params = {
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone_plus_hydrogens",
		"distance_constraints": "precise",
		"vdw_constraints": "yes",
		"tightening_upper_bounds": "yes",
		"torsion_angle_width": "0",
		"percentage_backbone_torsion_angles": "100",
		"max_distance": "5",
	}

	parsed = _preflight_parser_params(params)

	assert parsed["distance_constraints"] == "precise"
	assert parsed["vdw_constraints"] == "no"
	assert parsed["tightening_upper_bounds"] == "no"


def _interval_centered_params(**overrides: str) -> dict[str, str]:
	params = {
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone_plus_hydrogens",
		"distance_constraints": "interval_centered",
		"epsilon_long": "2.0",
		"max_distance": "5.0",
		"vdw_constraints": "no",
		"torsion_angle_width": "0",
		"percentage_backbone_torsion_angles": "100",
	}
	params.update(overrides)
	return params


def test_parser_preflight_tightening_upper_bounds_defaults_to_no() -> None:
	assert _preflight_parser_params(_interval_centered_params())["tightening_upper_bounds"] == "no"


def test_parser_preflight_tightening_upper_bounds_is_independent_of_vdw() -> None:
	parsed = _preflight_parser_params(_interval_centered_params(tightening_upper_bounds="yes"))

	assert parsed["vdw_constraints"] == "no"
	assert parsed["tightening_upper_bounds"] == "yes"


def test_parser_preflight_rejects_invalid_tightening_upper_bounds() -> None:
	with pytest.raises(ValueError, match="tightening_upper_bounds"):
		_preflight_parser_params(_interval_centered_params(tightening_upper_bounds="maybe"))


def test_parser_preflight_accepts_backbone_atom_selection() -> None:
	params = {
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone",
		"distance_constraints": "precise",
		"vdw_constraints": "yes",
		"torsion_angle_width": "0",
		"percentage_backbone_torsion_angles": "100",
		"max_distance": "5",
	}

	parsed = _preflight_parser_params(params)

	assert parsed["atom_selection"] == "backbone"


def test_parser_preflight_accepts_full_chain_atom_selection() -> None:
	params = {
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "full_chain",
		"distance_constraints": "precise",
		"vdw_constraints": "yes",
		"torsion_angle_width": "0",
		"percentage_backbone_torsion_angles": "100",
		"max_distance": "5",
	}

	_preflight_parser_params(params)


def test_parser_preflight_rejects_missing_nmr_only_experimental_parameters() -> None:
	params = {
		"method": "nmr",
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone",
		"distance_constraints": "interval_experimental",
		"vdw_constraints": "no",
		"torsion_angle_width": "0",
		"percentage_backbone_torsion_angles": "100",
	}

	with pytest.raises(ValueError, match="noe_strong"):
		_preflight_parser_params(params)


def test_parser_preflight_rejects_nonfinite_and_nonpositive_interval_values() -> None:
	params = {
		"method": "nmr",
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone",
		"distance_constraints": "interval_centered",
		"epsilon_short": "0",
		"epsilon_long": "2",
		"max_distance": "5",
		"vdwr_hh": "1.8",
		"vdw_constraints": "no",
		"torsion_angle_width": "0",
		"percentage_backbone_torsion_angles": "100",
	}

	with pytest.raises(ValueError, match="epsilon_short.*greater than zero"):
		_preflight_parser_params(params)

	params["epsilon_short"] = "1"
	params["max_distance"] = "nan"
	with pytest.raises(ValueError, match="max_distance.*finite"):
		_preflight_parser_params(params)


def test_parser_preflight_rejects_invalid_nmr_noe_bound_order() -> None:
	params = {
		"method": "nmr",
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone",
		"distance_constraints": "interval_experimental",
		"noe_strong": "3.5",
		"noe_medium": "2.5",
		"noe_weak": "5.0",
		"vdwr_hh": "1.8",
		"vdw_constraints": "no",
		"torsion_angle_width": "0",
		"percentage_backbone_torsion_angles": "100",
	}

	with pytest.raises(ValueError, match="vdwr_hh < noe_strong"):
		_preflight_parser_params(params)


def _auto_experimental_params(**overrides: str) -> dict[str, str]:
	params = {
		"method": "auto",
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone_plus_hydrogens",
		"distance_constraints": "interval_experimental",
		"epsilon_short": "1.0",
		"epsilon_long": "2.0",
		"max_distance": "5.0",
		"noe_strong": "2.5",
		"noe_medium": "3.5",
		"noe_weak": "5.0",
		"vdwr_hh": "1.8",
		"vdw_constraints": "no",
		"torsion_angle_width": "60",
		"percentage_backbone_torsion_angles": "100",
	}
	params.update(overrides)
	return params


@pytest.mark.parametrize("overrides", [{"noe_strong": "4.0"}, {"vdwr_hh": "2.5"}, {"noe_medium": "6.0"}])
def test_parser_preflight_checks_noe_order_for_auto_when_all_bounds_are_given(overrides: dict[str, str]) -> None:
	with pytest.raises(ValueError, match="vdwr_hh < noe_strong"):
		_preflight_parser_params(_auto_experimental_params(**overrides))


@pytest.mark.parametrize(
	("key", "value", "message"),
	[
		("epsilon_long", "-1", "epsilon_long.*greater than zero"),
		("epsilon_short", "0", "epsilon_short.*greater than zero"),
		("max_distance", "inf", "max_distance.*finite"),
	],
)
def test_parser_preflight_rejects_invalid_distances_unused_by_the_mode(key: str, value: str, message: str) -> None:
	# interval_experimental does not use epsilon_* or max_distance, but the file must still be valid.
	with pytest.raises(ValueError, match=message):
		_preflight_parser_params(_auto_experimental_params(**{key: value}))


def test_parser_preflight_rejects_negative_torsion_angle_width() -> None:
	with pytest.raises(ValueError, match="torsion_angle_width.*zero or greater"):
		_preflight_parser_params(_auto_experimental_params(torsion_angle_width="-10"))


def test_parser_preflight_accepts_the_shipped_auto_experimental_values() -> None:
	assert _preflight_parser_params(_auto_experimental_params())["noe_weak"] == 5.0


def test_parser_preflight_accepts_xray_experimental_fallback_without_nmr_parameters() -> None:
	params = {
		"method": "xray",
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone",
		"distance_constraints": "interval_experimental",
		"vdw_constraints": "no",
		"torsion_angle_width": "0",
		"percentage_backbone_torsion_angles": "100",
	}

	parsed = _preflight_parser_params(params)

	assert parsed["method"] == "xray"


def test_parser_preflight_rejects_unknown_public_parameter_key() -> None:
	params = {
		"model_number": "1",
		"chain_id": "A",
		"atom_selection": "backbone_plus_hydrogens",
		"distance_constraints": "precise",
		"vdw_constraints": "yes",
		"torsion_angle_width": "0",
		"percentage_backbone_torsion_angles": "100",
		"max_distance": "5",
		"hidden_feature": "1",
	}

	with pytest.raises(ValueError, match="Unsupported parser parameter\\(s\\): hidden_feature"):
		_preflight_parser_params(params)


def test_read_instance_remap_params_rejects_unknown_public_parameter_key(tmp_path: Path) -> None:
	params_file = write_text_file(
		tmp_path / "instance_remap.cfg",
		"""
		model_number: 1
		chain_id: A
		order_id: 9
		secret_mode: yes
		""",
	)

	with pytest.raises(ValueError, match="Unsupported remapping parameter\\(s\\): secret_mode"):
		read_instance_remap_params(params_file)


def test_parse_public_order_id_accepts_supported_values() -> None:
	assert parse_public_order_id("1") == 1
	assert parse_public_order_id("10") == 10


def test_parse_public_order_id_rejects_unsupported_value() -> None:
	with pytest.raises(ValueError, match="Invalid order_id"):
		parse_public_order_id("11")


def test_parse_public_order_choice_rejects_non_numeric_order() -> None:
	with pytest.raises(ValueError, match="order_id"):
		parse_public_order_choice("random")


def test_read_instance_remap_params_accepts_public_config(tmp_path: Path) -> None:
	params_file = write_text_file(
		tmp_path / "instance_remap.cfg",
		"""
		model_number: 1
		chain_id: A
		order_id: 9
		""",
	)

	params = read_instance_remap_params(params_file)

	assert params == {
		"model_number": "1",
		"chain_id": "A",
		"order_id": "9",
	}


def test_read_instance_remap_params_rejects_non_numeric_order(tmp_path: Path) -> None:
	params_file = write_text_file(
		tmp_path / "instance_remap.cfg",
		"""
		model_number: 1
		chain_id: A
		order_id: random
		""",
	)

	with pytest.raises(ValueError, match="order_id"):
		read_instance_remap_params(params_file)


def test_load_or_create_seed_bundle_creates_and_reuses_required_seed_file(tmp_path: Path) -> None:
	seed_dir = tmp_path / "seeds"

	first_bundle = load_or_create_seed_bundle(
		seed_dir=seed_dir,
		pdb_id="1ABC",
		model_number=1,
		chain_id="A",
		required_seed_fields={"seed_talos_angle_centers"},
	)

	assert first_bundle["seed_file"].exists()
	assert isinstance(first_bundle["seed_talos_angle_centers"], int)
	assert first_bundle["seed_interval_centered_distances"] is None
	assert first_bundle["seed_phi_psi_mask"] is None

	second_bundle = load_or_create_seed_bundle(
		seed_dir=seed_dir,
		pdb_id="1ABC",
		model_number=1,
		chain_id="A",
		required_seed_fields={"seed_talos_angle_centers"},
	)

	assert second_bundle["seed_talos_angle_centers"] == first_bundle["seed_talos_angle_centers"]


def test_load_or_create_seed_bundle_removes_seed_file_when_no_fields_are_required(tmp_path: Path) -> None:
	seed_dir = tmp_path / "seeds"

	initial_bundle = load_or_create_seed_bundle(
		seed_dir=seed_dir,
		pdb_id="1ABC",
		model_number=1,
		chain_id="A",
		required_seed_fields={"seed_phi_psi_mask"},
	)

	assert initial_bundle["seed_file"].exists()

	final_bundle = load_or_create_seed_bundle(
		seed_dir=seed_dir,
		pdb_id="1ABC",
		model_number=1,
		chain_id="A",
		required_seed_fields=set(),
	)

	assert not final_bundle["seed_file"].exists()


@pytest.mark.parametrize("value", ["0", "-1", str(2**63)])
def test_load_or_create_seed_bundle_rejects_invalid_seed_values(tmp_path: Path, value: str) -> None:
	seed_file = tmp_path / "seeds" / "1ABC" / "seed_1ABC.dat"
	seed_file.parent.mkdir(parents=True)
	seed_file.write_text(
		f"pdb_id: 1ABC\nmodel_number: 1\nchain_id: A\nseed_talos_angle_centers: {value}\n",
		encoding="utf-8",
	)
	with pytest.raises(ValueError, match="invalid value"):
		load_or_create_seed_bundle(tmp_path / "seeds", "1ABC", 1, "A", {"seed_talos_angle_centers"})

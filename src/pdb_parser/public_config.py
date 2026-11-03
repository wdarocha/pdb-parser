"""Public configuration surface for the current release.

This module centralizes the values intentionally exposed through
public CLIs and configuration files. Internal or future options may still
exist elsewhere in the codebase, but they are not part of the supported
release contract unless listed here.
"""

from __future__ import annotations

from pathlib import Path

from pdb_parser.utils.inputs import (
	normalize_nonempty_config_value,
	parse_int_config_value,
	parse_supported_int_config_value,
	reject_unknown_keys,
	require_keys,
)

VALID_DISTANCE_MODES = {
	"precise",
	"interval_centered",
	"interval_experimental",
}

VALID_METHODS = {"auto", "nmr", "xray"}
VALID_OUTPUT_MODES = {"filtered_atoms", "parser_outputs"}

INTERNAL_ATOM_SELECTIONS = {
	"full_chain",
	"backbone",
	"backbone_plus_hydrogens",
	"backbone_plus_neighbors",
}

PUBLIC_INPUT_ATOM_SELECTIONS = {
	"full_chain",
	"backbone",
	"backbone_plus_hydrogens",
	"backbone_plus_neighbors",
}

VALID_VDW_OPTIONS = {"yes", "no"}

PUBLIC_PARSER_PARAM_KEYS = {
	"method",
	"output_mode",
	"model_number",
	"chain_id",
	"atom_selection",
	"distance_constraints",
	"epsilon_short",
	"epsilon_long",
	"max_distance",
	"noe_strong",
	"noe_medium",
	"noe_weak",
	"vdwr_hh",
	"vdw_constraints",
	"tightening_upper_bounds",
	"torsion_angle_width",
	"percentage_backbone_torsion_angles",
	"percentage_hh_distance_constraints",
}

SUPPORTED_ORDER_IDS = set(range(1, 11))
DDGP_ORDER_TYPE = "ddgp"
REDMDGP_ORDER_TYPE = "redmdgp"
DMDGP_ORDER_TYPE = "dmdgp"
VALID_ORDER_TYPES = {DDGP_ORDER_TYPE, REDMDGP_ORDER_TYPE, DMDGP_ORDER_TYPE}
# order_type: dmdgp has a single fixed order (N CA C), so order_id is optional; this value
# stands in for it internally and never appears in an output name.
DMDGP_ORDER_ID = 1
PUBLIC_REMAP_PARAM_KEYS = {"model_number", "chain_id", "order_id", "order_type"}


def validate_public_remap_params(params: dict[str, str], path: Path) -> None:
	"""Validate the current public remapping configuration surface."""
	reject_unknown_keys(params, PUBLIC_REMAP_PARAM_KEYS, "remapping")
	order_types = parse_public_order_types(params)
	required = ("model_number", "chain_id") if order_types == [DMDGP_ORDER_TYPE] else ("model_number", "chain_id", "order_id")
	require_keys(params, required, path, "remapping")
	parse_int_config_value(params["model_number"], "model_number", "remapping")
	normalize_nonempty_config_value(params["chain_id"], "chain_id", "remapping")
	parse_public_remap_order_choices(params)


def parse_public_remap_order_choices(params: dict[str, str]) -> list[int]:
	"""Return the `order_id` values of a remap configuration.

	With `order_type: dmdgp` alone, `order_id` is optional (the DMDGP order is fixed); a
	given value is still validated, then ignored.
	"""
	if parse_public_order_types(params) == [DMDGP_ORDER_TYPE]:
		if str(params.get("order_id", "")).strip():
			parse_public_order_choices(params["order_id"])
		return [DMDGP_ORDER_ID]
	return parse_public_order_choices(params["order_id"])


def parse_public_order_type(params: dict[str, str]) -> str:
	"""Return the optional `order_type` (`ddgp` by default) for the remapping stage.

	`redmdgp` rewrites the same `order_id` (1 to 10) as an iDMDGP re-order; `dmdgp` requests
	the fixed DMDGP order N CA C of `backbone` outputs.
	"""
	value = str(params.get("order_type", "")).strip().lower() or DDGP_ORDER_TYPE
	if value not in VALID_ORDER_TYPES:
		raise ValueError(
			f"Invalid order_type {value!r}. Supported values in the remapping params file are: "
			f"{', '.join(sorted(VALID_ORDER_TYPES))}."
		)
	return value


def parse_public_order_types(params: dict[str, str]) -> list[str]:
	"""Return the `order_type` values, given one or comma-separated (`ddgp` by default).

	Repeated values are kept once, in their first position.
	"""
	raw_value = str(params.get("order_type", ""))
	if "," not in raw_value:
		return [parse_public_order_type(params)]
	order_types: list[str] = []
	for item in raw_value.split(","):
		order_type = parse_public_order_type({"order_type": item})
		if not item.strip():
			raise ValueError("order_type cannot contain an empty item in the remapping params file.")
		if order_type not in order_types:
			order_types.append(order_type)
	if DMDGP_ORDER_TYPE in order_types:
		raise ValueError(
			"order_type: dmdgp cannot be combined with ddgp or redmdgp: it applies to backbone "
			"outputs, the others to backbone_plus_hydrogens outputs."
		)
	return order_types


def parse_public_order_choices(raw_value: str) -> list[int]:
	"""Parse one `order_id`, or a list of values and ranges such as `1-10` or `1,3,9`.

	Repeated values are kept once, in their first position. A single value is validated
	exactly as by :func:`parse_public_order_choice`.
	"""
	order_ids: list[int] = []
	for item in raw_value.split(","):
		bounds = item.split("-")
		if len(bounds) == 2 and bounds[0].strip() and bounds[1].strip():
			start, end = (parse_public_order_choice(bound) for bound in bounds)
			if start > end:
				raise ValueError(
					f"Invalid order_id range {item.strip()!r} in the remapping params file: "
					"the start must not exceed the end."
				)
			values = range(start, end + 1)
		else:
			if not item.strip() and "," in raw_value:
				raise ValueError("order_id cannot contain an empty item in the remapping params file.")
			values = [parse_public_order_choice(item)]
		order_ids.extend(value for value in values if value not in order_ids)
	return order_ids


def parse_public_order_choice(raw_value: str) -> int:
	"""Parse and validate a public `order_id` value for the current release."""
	return parse_supported_int_config_value(
		raw_value.strip(),
		"order_id",
		"remapping",
		SUPPORTED_ORDER_IDS,
	)


def parse_public_order_id(raw_value: str) -> int:
	"""Parse and validate a numeric public `order_id` value."""
	return parse_public_order_choice(raw_value)

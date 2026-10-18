"""Parser parameter validation and policy helpers."""

from __future__ import annotations

import math

from pdb_parser.public_config import (
    PUBLIC_INPUT_ATOM_SELECTIONS,
    PUBLIC_PARSER_PARAM_KEYS,
    VALID_DISTANCE_MODES,
    VALID_METHODS,
    VALID_OUTPUT_MODES,
    VALID_VDW_OPTIONS,
)

from .contracts import ParserParams


def resolve_effective_angular_width(
    distance_constraints: str,
    torsion_angle_width: float,
) -> float:
    """Return the torsion width effective under the selected distance mode."""
    return 0.0 if distance_constraints == "precise" else max(torsion_angle_width, 360.0 if torsion_angle_width < 0.0 else torsion_angle_width)


def resolve_effective_backbone_torsion_percentage(
    distance_constraints: str, percentage_backbone_torsion_angles: float
) -> float:
    """Return the phi/psi selection percentage effective under parser policy."""
    if distance_constraints == "precise":
        return 100.0
    return min(max(percentage_backbone_torsion_angles, 0.0), 100.0)


def required_seed_fields_for_execution(
    distance_constraints: str,
    angular_width: float,
    backbone_torsion_percentage: float,
    n_residues: int,
    method: str = "nmr",
    hh_distance_percentage: float = 100.0,
) -> set[str]:
    """Return exactly the seed fields that can affect one parser execution."""
    fields: set[str] = set()
    if distance_constraints == "interval_centered":
        fields.add("seed_interval_centered_distances")
    if angular_width > 0.0:
        fields.add("seed_talos_angle_centers")
    detected = int(round(n_residues * backbone_torsion_percentage / 100.0))
    if 0 < n_residues - detected < n_residues:
        fields.add("seed_phi_psi_mask")
    if method == "nmr" and 0.0 < hh_distance_percentage < 100.0:
        fields.add("seed_hh_distance_subset")
    return fields


def reject_unknown_parser_param_keys(params: dict[str, object], allowed_keys: set[str]) -> None:
    """Reject configuration keys outside the public parser contract."""
    unknown_keys = sorted(set(params) - allowed_keys)
    if unknown_keys:
        raise ValueError(
            f"Unsupported parser parameter(s): {', '.join(unknown_keys)}. "
            f"Supported keys in this release are: {', '.join(sorted(allowed_keys))}."
        )

NOE_KEYS = ("noe_strong", "noe_medium", "noe_weak", "vdwr_hh")
PHYSICAL_DISTANCE_KEYS = ("epsilon_short", "epsilon_long", "max_distance", *NOE_KEYS)


def _require_param(params: ParserParams, key: str) -> str:
    """Return a stripped parameter value or raise a clear configuration error."""
    raw_value = params.get(key)
    if raw_value is None:
        raise ValueError(f"Missing required parser parameter: {key}")

    value = str(raw_value).strip()
    if not value:
        raise ValueError(f"Parser parameter '{key}' cannot be empty.")

    return value


def _parse_int_param(params: ParserParams, key: str) -> int:
    """Parse an integer parser parameter with a user-facing error message."""
    value = _require_param(params, key)
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"Parser parameter '{key}' must be an integer, got {value!r}.") from exc


def _parse_float_param(params: ParserParams, key: str) -> float:
    """Parse a float parser parameter with a user-facing error message."""
    value = _require_param(params, key)
    try:
        parsed_value = float(value)
    except ValueError as exc:
        raise ValueError(f"Parser parameter '{key}' must be a number, got {value!r}.") from exc
    if not math.isfinite(parsed_value):
        raise ValueError(f"Parser parameter '{key}' must be a finite number, got {value!r}.")
    return parsed_value


def _require_positive_param(parsed: ParserParams, key: str) -> None:
    """Require a parsed physical distance parameter to be strictly positive."""
    value = float(parsed[key])
    if value <= 0.0:
        raise ValueError(f"Parser parameter '{key}' must be greater than zero.")


def _validate_noe_parameters(parsed: ParserParams) -> None:
    """Validate the ordered physical bounds used by NMR experimental intervals."""
    for key in ("noe_strong", "noe_medium", "noe_weak", "vdwr_hh"):
        _require_positive_param(parsed, key)
    if not (
        float(parsed["vdwr_hh"]) < float(parsed["noe_strong"])
        <= float(parsed["noe_medium"])
        <= float(parsed["noe_weak"])
    ):
        raise ValueError(
        "NMR experimental bounds must satisfy "
        "vdwr_hh < noe_strong <= noe_medium <= noe_weak."
    )


def _parse_choice_param(
    params: ParserParams,
    key: str,
    valid_values: set[str],
) -> str:
    """Parse a normalized choice parameter."""
    value = _require_param(params, key).lower()
    if value not in valid_values:
        raise ValueError(
            f"Unsupported value for parser parameter '{key}': {value!r}. "
            f"Expected one of {sorted(valid_values)}."
    )
    return value


def _parse_optional_choice_param(
    params: ParserParams,
    key: str,
    valid_values: set[str],
    default: str,
) -> str:
    """Parse an optional choice while retaining old configuration files."""
    if key not in params or not str(params[key]).strip():
        return default
    return _parse_choice_param(params, key, valid_values)




def _parse_float_or_default(params: ParserParams, key: str, default: float) -> float:
    """Read an optional numeric field used by X-ray mode fallback."""
    if key not in params or not str(params[key]).strip():
        return default
    return _parse_float_param(params, key)


def _reject_unknown_param_keys(
    params: ParserParams,
    allowed_keys: set[str],
    context_label: str,
) -> None:
    """Reject unknown keys so internal features stay unavailable via public config."""
    if context_label == "parser":
        reject_unknown_parser_param_keys(params, allowed_keys)
        return
    unknown_keys = sorted(set(params) - allowed_keys)
    if unknown_keys:
        raise ValueError(f"Unsupported {context_label} parameter(s): {', '.join(unknown_keys)}. " f"Supported keys in this release are: {', '.join(sorted(allowed_keys))}.")


def _preflight_parser_params(params: ParserParams) -> ParserParams:
    """Validate parser parameters before creating output artifacts."""
    _reject_unknown_param_keys(params, PUBLIC_PARSER_PARAM_KEYS, "parser")

    parsed: ParserParams = {
        "method": _parse_optional_choice_param(params, "method", VALID_METHODS, "auto"),
        "output_mode": _parse_optional_choice_param(
            params, "output_mode", VALID_OUTPUT_MODES, "parser_outputs"
        ),
        "model_number": _parse_int_param(params, "model_number"),
        "chain_id": _require_param(params, "chain_id"),
        "atom_selection": _parse_choice_param(
            params,
            "atom_selection",
            PUBLIC_INPUT_ATOM_SELECTIONS,
        ),
        "distance_constraints": _parse_choice_param(
            params, "distance_constraints", VALID_DISTANCE_MODES
        ),
        "vdw_constraints": _parse_choice_param(params, "vdw_constraints", VALID_VDW_OPTIONS),
        "tightening_upper_bounds": _parse_optional_choice_param(
            params, "tightening_upper_bounds", VALID_VDW_OPTIONS, "no"
        ),
        "torsion_angle_width": _parse_float_param(params, "torsion_angle_width"),
        "percentage_backbone_torsion_angles": _parse_float_param(
            params, "percentage_backbone_torsion_angles"
        ),
        "percentage_hh_distance_constraints": _parse_float_or_default(
            params, "percentage_hh_distance_constraints", 100.0
        ),
    }

    # Every physical distance present in the file must be valid, even when the
    # selected mode does not use it.
    for key in PHYSICAL_DISTANCE_KEYS:
        if key in params and str(params[key]).strip():
            if _parse_float_param(params, key) <= 0.0:
                raise ValueError(f"Parser parameter '{key}' must be greater than zero.")
    if float(parsed["torsion_angle_width"]) < 0.0:
        raise ValueError("Parser parameter 'torsion_angle_width' must be zero or greater.")

    distance_constraints = parsed["distance_constraints"]
    if distance_constraints == "precise":
        # Precise distances are exact: no VDW bounds to add, no upper bound to tighten.
        parsed["vdw_constraints"] = "no"
        parsed["tightening_upper_bounds"] = "no"

    if distance_constraints in {"precise", "interval_centered"}:
        parsed["max_distance"] = _parse_float_param(params, "max_distance")
        _require_positive_param(parsed, "max_distance")

    if distance_constraints == "interval_centered":
        if "epsilon_short" in params and str(params["epsilon_short"]).strip():
            parsed["epsilon_short"] = _parse_float_param(params, "epsilon_short")
        parsed["epsilon_long"] = _parse_float_param(params, "epsilon_long")
        _require_positive_param(parsed, "epsilon_long")
        if "vdwr_hh" in params and str(params["vdwr_hh"]).strip():
            parsed["vdwr_hh"] = _parse_float_param(params, "vdwr_hh")
            _require_positive_param(parsed, "vdwr_hh")
        if "epsilon_short" in parsed:
            _require_positive_param(parsed, "epsilon_short")

    elif distance_constraints == "interval_experimental":
        for key in ("noe_strong", "noe_medium", "noe_weak", "vdwr_hh"):
            if key in params and str(params[key]).strip():
                parsed[key] = _parse_float_param(params, key)

    # An explicit NMR configuration can be fully validated before opening a
    # structure.  ``auto`` remains intentionally deferred to EXPDTA detection.
    if parsed["method"] == "nmr" and distance_constraints == "interval_centered":
        parsed["epsilon_short"] = _parse_float_param(params, "epsilon_short")
        parsed["vdwr_hh"] = _parse_float_param(params, "vdwr_hh")
        _require_positive_param(parsed, "epsilon_short")
        _require_positive_param(parsed, "vdwr_hh")
    if parsed["method"] == "nmr" and distance_constraints == "interval_experimental":
        for key in NOE_KEYS:
            parsed[key] = _parse_float_param(params, key)
        _validate_noe_parameters(parsed)
    elif (
        parsed["method"] == "auto"
        and distance_constraints == "interval_experimental"
        and all(key in parsed for key in NOE_KEYS)
    ):
        # ``auto`` may still receive X-ray entries, which do not need these keys;
        # when they are all given, their order is checked before any structure.
        _validate_noe_parameters(parsed)

    percentage = float(parsed["percentage_backbone_torsion_angles"])
    if not 0.0 <= percentage <= 100.0:
        raise ValueError(
            "Parser parameter 'percentage_backbone_torsion_angles' must be between 0 and 100."
        )
    if not 0.0 <= float(parsed["percentage_hh_distance_constraints"]) <= 100.0:
        raise ValueError(
            "Parser parameter 'percentage_hh_distance_constraints' must be between 0 and 100."
        )

    return parsed


def validate_public_parser_params(params: ParserParams) -> None:
    """Validate the current public parser configuration surface.

    This helper exists so public CLIs can fail fast on unsupported parser
    configuration before any per-structure work starts.
    """
    _preflight_parser_params(params)

"""Input readers for pdb-parser: PDB id list and params file."""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path

from pdb_parser.messages import print_warning

_PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")


def read_pdb_ids(path: Path) -> list[str]:
	"""Read PDB ids from a text file, validate, and deduplicate preserving order.
	
	Rules:
	- Lines starting with '#' are comments.
	- Inline comments after '#' are removed.
	- Blank lines are ignored.
	- IDs are uppercased and must match ^[A-Za-z0-9]{4}$.
	"""
	if not path.exists():
		raise FileNotFoundError(f"PDB id list not found: {path}")

	seen: set[str] = set()
	result: list[str] = []

	with path.open("r", encoding="utf-8", errors="replace") as fh:
		for lineno, raw in enumerate(fh, start=1):
			line = raw.strip()
			if not line or line.startswith("#"):
				continue
			if "#" in line:
				line = line.split("#", 1)[0].strip()
			if not line:
				continue

			pid = line.upper()
			if not _PDB_ID_RE.match(pid):
				raise ValueError(
					f"Invalid PDB id at {path}:{lineno} -> '{line}'. "
					"Expected exactly 4 alphanumeric characters."
				)
			if pid not in seen:
				seen.add(pid)
				result.append(pid)

	if not result:
		raise ValueError(f"No valid PDB ids found in {path}.")
	return result


def read_params(path: Path) -> dict[str, str]:
	"""Read a key: value parameter file.

	Rules
	- Lines starting with '#' are comments.
	- Inline comments after '#' are removed.
	- Keys are lowercased and stripped; values are stripped (kept as strings).
	- Empty lines are ignored.
	- Only the ':' separator is accepted.
	"""
	if not path.exists():
		raise FileNotFoundError(f"Params file not found: {path}")

	params: dict[str, str] = {}
	with path.open("r", encoding="utf-8", errors="replace") as fh:
		for lineno, raw in enumerate(fh, start=1):
			line = raw.strip()
			if not line or line.startswith("#"):
				continue
			if "#" in line:
				line = line.split("#", 1)[0].strip()
			if not line:
				continue
			if ":" not in line:
				raise ValueError(
					f"Invalid params line at {path}:{lineno} -> '{raw.rstrip()}'. "
					"Expected 'key: value'."
				)
			key, _, value = line.partition(":")
			k = key.strip().lower()
			v = value.strip()
			if not k:
				raise ValueError(f"Empty key at {path}:{lineno}.")
			if k in params:
				raise ValueError(
					f"Duplicate parameter '{k}' at {path}:{lineno}. "
					"Each configuration key may be specified only once."
				)
			params[k] = v

	if not params:
		print_warning(f"Params file {path} parsed empty.")
	return params


def reject_unknown_keys(
	params: dict[str, str],
	allowed_keys: Iterable[str],
	context_label: str,
) -> None:
	"""Reject unsupported public config keys with a clear user-facing error."""
	allowed_key_set = set(allowed_keys)
	unknown_keys = sorted(set(params) - allowed_key_set)
	if unknown_keys:
		raise ValueError(
			f"Unsupported {context_label} parameter(s): {', '.join(unknown_keys)}. "
			f"Supported keys in this release are: {', '.join(sorted(allowed_key_set))}."
		)


def require_keys(
	params: dict[str, str],
	required_keys: Iterable[str],
	path: Path,
	context_label: str,
) -> None:
	"""Require the given config keys to exist in one public params file."""
	missing_keys = [key for key in required_keys if key not in params]
	if missing_keys:
		raise ValueError(
			f"Missing required {context_label} parameter(s) in {path}: {', '.join(missing_keys)}"
		)


def parse_int_config_value(
	raw_value: str,
	field_name: str,
	context_label: str,
) -> int:
	"""Parse one public config value as integer with a user-facing error."""
	try:
		return int(raw_value)
	except ValueError as exc:
		raise ValueError(f"{field_name} must be an integer in the {context_label} params file.") from exc


def parse_supported_int_config_value(
	raw_value: str,
	field_name: str,
	context_label: str,
	supported_values: Iterable[int],
) -> int:
	"""Parse one integer config value and require it to belong to the supported set."""
	value = parse_int_config_value(raw_value, field_name, context_label)
	supported_value_set = set(supported_values)
	if value not in supported_value_set:
		supported_text = ", ".join(str(item) for item in sorted(supported_value_set))
		raise ValueError(
			f"Invalid {field_name}. Supported values in the {context_label} params file are: {supported_text}."
		)
	return value


def normalize_nonempty_config_value(
	raw_value: str,
	field_name: str,
	context_label: str,
) -> str:
	"""Normalize one public config value and reject empty strings."""
	value = raw_value.strip()
	if not value:
		raise ValueError(f"{field_name} cannot be empty in the {context_label} params file.")
	return value

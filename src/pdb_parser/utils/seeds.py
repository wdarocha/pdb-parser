"""Seed-file helpers for reproducible per-PDB random streams."""

from __future__ import annotations

import os
import secrets
import tempfile
from pathlib import Path

import numpy as np

from pdb_parser.messages import print_ok

from .fs import ensure_dir
from .inputs import read_params

SeedRecord = dict[str, str]
SeedBundle = dict[str, Path | int | None]


SEED_FIELD_ORDER = (
	"seed_interval_centered_distances",
	"seed_talos_angle_centers",
	"seed_phi_psi_mask",
	"seed_hh_distance_subset",
)
def _generate_seed_value() -> int:
	"""Return a positive 63-bit seed suitable for NumPy Generator."""
	return secrets.randbelow((1 << 63) - 1) + 1


def _seed_file_path(seed_dir: str | Path, pdb_id: str) -> Path:
	"""Return the canonical per-PDB seed file path inside the seed directory."""
	return Path(seed_dir) / pdb_id / f"seed_{pdb_id}.dat"


def _validate_seed_file_metadata(
	record: SeedRecord,
	seed_file: Path,
	pdb_id: str,
	model_number: int,
	chain_id: str,
) -> None:
	"""Validate that one existing seed file belongs to the requested run context."""
	required_keys = ("pdb_id", "model_number", "chain_id")
	missing_keys = [key for key in required_keys if key not in record]
	if missing_keys:
		raise ValueError(
			f"Seed file {seed_file} is missing required field(s): {', '.join(missing_keys)}."
		)

	stored_pdb_id = record["pdb_id"].strip().upper()
	stored_chain_id = record["chain_id"].strip()

	try:
		stored_model_number = int(record["model_number"])
	except ValueError as exc:
		raise ValueError(
			f"Seed file {seed_file} has invalid model_number: {record['model_number']!r}."
		) from exc

	if stored_pdb_id != pdb_id.upper():
		raise ValueError(
			f"Seed file {seed_file} belongs to pdb_id={stored_pdb_id}, "
			f"but the current run requested {pdb_id.upper()}."
		)

	if stored_model_number != model_number or stored_chain_id != chain_id:
		raise ValueError(
			f"Seed file {seed_file} belongs to model={stored_model_number}, chain={stored_chain_id!r}, "
			f"but the current run requested model={model_number}, chain={chain_id!r}."
		)


def _write_seed_file(seed_file: Path, record: SeedRecord) -> None:
	"""Write the normalized seed record to disk in the public key-value format."""
	ensure_dir(seed_file.parent)
	fd, temporary_name = tempfile.mkstemp(prefix=f".{seed_file.name}.", dir=seed_file.parent)
	try:
		with os.fdopen(fd, "w", encoding="utf-8") as fh:
			fh.write(f"pdb_id: {record['pdb_id']}\n")
			fh.write(f"model_number: {record['model_number']}\n")
			fh.write(f"chain_id: {record['chain_id']}\n")
			for key in SEED_FIELD_ORDER:
				if key in record:
					fh.write(f"{key}: {record[key]}\n")
			fh.flush()
			os.fsync(fh.fileno())
		os.replace(temporary_name, seed_file)
	finally:
		Path(temporary_name).unlink(missing_ok=True)


def _validate_seed_value(raw_value: str, key: str, seed_file: Path) -> int:
	"""Validate the positive integer range accepted by NumPy's seed handling."""
	try:
		value = int(raw_value)
	except ValueError as exc:
		raise ValueError(
			f"Seed file {seed_file} has invalid integer value for {key}: {raw_value!r}."
		) from exc
	if not 1 <= value <= np.iinfo(np.int64).max:
		raise ValueError(
			f"Seed file {seed_file} has invalid value for {key}: expected an integer in 1..{np.iinfo(np.int64).max}."
		)
	return value


def load_or_create_seed_bundle(
	seed_dir: str | Path,
	pdb_id: str,
	model_number: int,
	chain_id: str,
	required_seed_fields: set[str] | list[str] | tuple[str, ...],
) -> SeedBundle:
	"""Load or create only the per-PDB random seeds required by this execution."""
	seed_file = _seed_file_path(seed_dir, pdb_id)
	seed_file_existed = seed_file.exists()
	required_fields = {
		field
		for field in required_seed_fields
		if field in SEED_FIELD_ORDER
	}

	invalid_fields = sorted(set(required_seed_fields) - set(SEED_FIELD_ORDER))
	if invalid_fields:
		raise ValueError(
			f"Unsupported seed field(s) requested for {pdb_id}: {', '.join(invalid_fields)}."
		)

	if seed_file_existed:
		record = read_params(seed_file)
		_validate_seed_file_metadata(record, seed_file, pdb_id, model_number, chain_id)
	else:
		record: SeedRecord = {
			"pdb_id": pdb_id.upper(),
			"model_number": str(model_number),
			"chain_id": chain_id,
		}

	missing_fields: list[str] = []
	existing_fields: set[str] = set()
	for key in SEED_FIELD_ORDER:
		raw_value = record.get(key, "").strip()
		if not raw_value:
			record.pop(key, None)
			continue

		_validate_seed_value(raw_value, key, seed_file)

		existing_fields.add(key)

	for key in sorted(required_fields - existing_fields):
		record[key] = str(_generate_seed_value())
		missing_fields.append(key)

	removed_fields = sorted(existing_fields - required_fields)
	for key in removed_fields:
		record.pop(key, None)

	if required_fields:
		if not seed_file_existed or missing_fields or removed_fields:
			_write_seed_file(seed_file, record)
			if seed_file_existed:
				changes: list[str] = []
				if missing_fields:
					changes.append(f"added: {', '.join(sorted(missing_fields))}")
				if removed_fields:
					changes.append(f"removed: {', '.join(removed_fields)}")
				change_summary = "; ".join(changes)
				print_ok(
					f"Seed file updated at: {seed_file}"
					f"{f' ({change_summary})' if change_summary else ''}"
				)
			else:
				print_ok(f"Seed file created at: {seed_file}")
		else:
			print_ok(f"Seed file reused from: {seed_file}")
	elif seed_file_existed:
		seed_file.unlink()
		print_ok(f"Seed file removed at: {seed_file} (no seeds needed for this run)")
	else:
		print_ok(f"No seed file needed for this run: {pdb_id}")

	return {
		"seed_file": seed_file,
		"seed_interval_centered_distances": (
			int(record["seed_interval_centered_distances"])
			if "seed_interval_centered_distances" in record
			else None
		),
		"seed_talos_angle_centers": (
			int(record["seed_talos_angle_centers"])
			if "seed_talos_angle_centers" in record
			else None
		),
		"seed_phi_psi_mask": (
			int(record["seed_phi_psi_mask"])
			if "seed_phi_psi_mask" in record
			else None
		),
		"seed_hh_distance_subset": (
			int(record["seed_hh_distance_subset"])
			if "seed_hh_distance_subset" in record
			else None
		),
	}

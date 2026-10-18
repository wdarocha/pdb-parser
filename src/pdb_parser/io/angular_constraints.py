"""Convert backbone torsion intervals into distance constraints."""

from __future__ import annotations

from pathlib import Path
from typing import TypedDict

import numpy as np
import pandas as pd

from pdb_parser.geometry.distance_geometry import torsion_angle_2_endpoint_distance
from pdb_parser.io.atom_tables import load_filtered_atoms_table, validate_filtered_atoms_table
from pdb_parser.io.files import distance_table
from pdb_parser.io.sequence import residue_sequence_from_ids
from pdb_parser.messages import print_ok


class BackboneAtomRecord(TypedDict):
	atom_id: int
	coord: np.ndarray


class BackboneResidueRecord(TypedDict):
	resname: str
	atoms: dict[str, BackboneAtomRecord]


BackboneResidueMap = dict[int, BackboneResidueRecord]
DistanceRecord = tuple[int, int, int, int, float, float, str, str, str, str]


def wrap_angle_deg(angle: float) -> float:
	"""Wrap an angle in degrees to the range [-180, 180)."""
	return (angle + 180.0) % 360.0 - 180.0


def full_range_angle() -> tuple[float, float]:
	"""Return a full-range torsion interval representing [-180, 180]."""
	return 0.0, 180.0


def circular_distance_deg(angle_a: float, angle_b: float) -> float:
	"""Return the shortest circular distance between two angles in degrees."""
	return abs(wrap_angle_deg(angle_a - angle_b))


def angle_belongs_to_interval(target: float, center: float, radius: float, tol: float = 1e-12) -> bool:
	"""Check whether a target angle belongs to a circular interval."""
	return circular_distance_deg(target, center) <= radius + tol


def torsion_interval_distance_bounds(
	x1: np.ndarray,
	x2: np.ndarray,
	x3: np.ndarray,
	x4: np.ndarray,
	center_deg: float,
	radius_deg: float,
) -> tuple[float, float]:
	"""Return endpoint-distance bounds induced by a circular torsion interval."""
	center_deg = wrap_angle_deg(float(center_deg))
	radius_deg = float(radius_deg)
	if radius_deg < 0.0:
		raise ValueError(f"radius_deg must be nonnegative, got {radius_deg}.")

	if radius_deg >= 180.0:
		candidate_angles_deg = [0.0, 180.0]
	else:
		candidate_angles_deg = [
			wrap_angle_deg(center_deg - radius_deg),
			wrap_angle_deg(center_deg + radius_deg),
		]
		for critical_angle in (0.0, 180.0, -180.0):
			if angle_belongs_to_interval(critical_angle, center_deg, radius_deg):
				candidate_angles_deg.append(critical_angle)

	distances = [
		torsion_angle_2_endpoint_distance(x1, x2, x3, x4, np.deg2rad(angle_deg))
		for angle_deg in candidate_angles_deg
	]
	return min(distances), max(distances)


def build_backbone_residue_map(df_atoms: pd.DataFrame) -> BackboneResidueMap:
	"""Build a residue-indexed dictionary with N, CA, and C coordinates."""
	validate_filtered_atoms_table(df_atoms, "<in-memory filtered atom table>", require_nonempty=True)
	residue_map: BackboneResidueMap = {}
	for row in df_atoms.itertuples(index=False):
		atom_name = str(row.atom_name).strip().upper()
		if atom_name not in {"N", "CA", "C"}:
			continue
		resid = int(row.resid)
		if resid not in residue_map:
			residue_map[resid] = {"resname": str(row.resname).strip().upper(), "atoms": {}}
		residue_map[resid]["atoms"][atom_name] = {
			"atom_id": int(row.atom_id),
			"coord": np.array([float(row.x), float(row.y), float(row.z)], dtype=float),
		}
	return residue_map


def _adjacent_residue(residue_map: BackboneResidueMap, resid: int, step: int) -> int | None:
	"""Return the actual neighboring residue, independent of the starting number."""
	residue_ids = residue_sequence_from_ids(residue_map, context="filtered atom table")
	try:
		index = residue_ids.index(int(resid)) + step
	except ValueError:
		return None
	return residue_ids[index] if 0 <= index < len(residue_ids) else None


def get_backbone_atom(
	residue_map: BackboneResidueMap, resid: int, atom_name: str
) -> BackboneAtomRecord | None:
	"""Return the stored atom record for a residue and backbone atom name."""
	residue_data = residue_map.get(int(resid))
	if residue_data is None:
		return None
	return residue_data["atoms"].get(atom_name)


def omega_interval_to_distance_record(
	residue_map: BackboneResidueMap, resid: int, resname: str, omega_center: float, omega_radius: float
) -> DistanceRecord | None:
	"""Convert omega_i to the CA_(i-1)--CA_i distance interval."""
	previous = _adjacent_residue(residue_map, resid, -1)
	if previous is None:
		return None
	ca_prev = get_backbone_atom(residue_map, previous, "CA")
	c_prev = get_backbone_atom(residue_map, previous, "C")
	n_i = get_backbone_atom(residue_map, resid, "N")
	ca_i = get_backbone_atom(residue_map, resid, "CA")
	if ca_prev is None or c_prev is None or n_i is None or ca_i is None:
		return None
	d_l, d_u = torsion_interval_distance_bounds(ca_prev["coord"], c_prev["coord"], n_i["coord"], ca_i["coord"], omega_center, omega_radius)
	return (int(ca_i["atom_id"]), int(ca_prev["atom_id"]), resid, previous, d_l, d_u, "CA", "CA", str(resname).upper(), residue_map[previous]["resname"].upper())


def phi_interval_to_distance_record(
	residue_map: BackboneResidueMap, resid: int, resname: str, phi_center: float, phi_radius: float
) -> DistanceRecord | None:
	"""Convert phi_i to the C_(i-1)--C_i distance interval."""
	previous = _adjacent_residue(residue_map, resid, -1)
	if previous is None:
		return None
	c_prev = get_backbone_atom(residue_map, previous, "C")
	n_i = get_backbone_atom(residue_map, resid, "N")
	ca_i = get_backbone_atom(residue_map, resid, "CA")
	c_i = get_backbone_atom(residue_map, resid, "C")
	if c_prev is None or n_i is None or ca_i is None or c_i is None:
		return None
	d_l, d_u = torsion_interval_distance_bounds(c_prev["coord"], n_i["coord"], ca_i["coord"], c_i["coord"], phi_center, phi_radius)
	return (int(c_i["atom_id"]), int(c_prev["atom_id"]), resid, previous, d_l, d_u, "C", "C", str(resname).upper(), residue_map[previous]["resname"].upper())


def psi_interval_to_distance_record(
	residue_map: BackboneResidueMap, resid: int, resname: str, psi_center: float, psi_radius: float
) -> DistanceRecord | None:
	"""Convert psi_i to the N_i--N_(i+1) distance interval."""
	n_i = get_backbone_atom(residue_map, resid, "N")
	ca_i = get_backbone_atom(residue_map, resid, "CA")
	c_i = get_backbone_atom(residue_map, resid, "C")
	next_residue = _adjacent_residue(residue_map, resid, 1)
	if next_residue is None:
		return None
	n_next = get_backbone_atom(residue_map, next_residue, "N")
	if n_i is None or ca_i is None or c_i is None or n_next is None:
		return None
	d_l, d_u = torsion_interval_distance_bounds(n_i["coord"], ca_i["coord"], c_i["coord"], n_next["coord"], psi_center, psi_radius)
	return (int(n_next["atom_id"]), int(n_i["atom_id"]), next_residue, resid, d_l, d_u, "N", "N", residue_map[next_residue]["resname"].upper(), str(resname).upper())


def backbone_angular_interval_to_distance_interval(
	tsv_structure_file: str | Path,
	tsv_angular_file: str | Path,
	output_distance_file: str | Path,
) -> pd.DataFrame:
	"""Convert omega, phi and psi intervals into backbone distance intervals (also returned)."""
	df_atoms = load_filtered_atoms_table(tsv_structure_file)
	validate_filtered_atoms_table(df_atoms, tsv_structure_file, require_nonempty=True)
	angular_path = Path(tsv_angular_file)
	output_path = Path(output_distance_file)
	if not angular_path.exists():
		raise ValueError(f"Angular-interval file not found: {angular_path}")
	df_angles = pd.read_csv(angular_path, sep="\t", dtype={"resname": str}, float_precision="round_trip")
	df_angles.columns = df_angles.columns.str.strip()
	residue_map = build_backbone_residue_map(df_atoms)
	records: list[DistanceRecord] = []

	for row in df_angles.itertuples(index=False):
		resid = int(row.resid)
		resname = str(row.resname).strip().upper()
		for record in (
			omega_interval_to_distance_record(residue_map, resid, resname, float(row.omega_center), float(row.omega_radius)),
			phi_interval_to_distance_record(residue_map, resid, resname, float(row.phi_center), float(row.phi_radius)),
			psi_interval_to_distance_record(residue_map, resid, resname, float(row.psi_center), float(row.psi_radius)),
		):
			if record is not None:
				records.append(record)

	with output_path.open("w", encoding="utf-8") as output_file:
		output_file.write("atom_id_i\tatom_id_j\tresid_i\tresid_j\td_l\td_u\tatom_name_i\tatom_name_j\tresname_i\tresname_j\n")
		for record in records:
			output_file.write(
				f"{record[0]}\t{record[1]}\t{record[2]}\t{record[3]}\t{record[4]:.12f}\t{record[5]:.12f}\t"
				f"{record[6]}\t{record[7]}\t{record[8]}\t{record[9]}\n"
			)

	print_ok(f"Backbone angular intervals converted to distance intervals saved to: {output_path}")
	return distance_table(
		(int(r[0]), int(r[1]), int(r[2]), int(r[3]), float(r[4]), float(r[5]), r[6], r[7], r[8], r[9]) for r in records
	)

"""Compatibility facade for distance-constraint helpers.

New code should import each helper from its dedicated module. This module keeps
the historical import path stable for downstream callers.
"""

from pdb_parser.io.angular_constraints import (
	BackboneAtomRecord,
	BackboneResidueMap,
	BackboneResidueRecord,
	DistanceRecord,
	angle_belongs_to_interval,
	backbone_angular_interval_to_distance_interval,
	build_backbone_residue_map,
	circular_distance_deg,
	full_range_angle,
	get_backbone_atom,
	omega_interval_to_distance_record,
	phi_interval_to_distance_record,
	psi_interval_to_distance_record,
	torsion_interval_distance_bounds,
	wrap_angle_deg,
)
from pdb_parser.io.angular_sampling import talos_n_like
from pdb_parser.io.atom_tables import load_filtered_atoms_table
from pdb_parser.io.constraint_files import (
	MergedDistanceRecord,
	merge_distance_constraint_files,
	tightening_upper_bounds,
)
from pdb_parser.io.distance_intervals import (
	get_centered_interval as _get_centered_interval,
)
from pdb_parser.io.distance_intervals import (
	get_experimental_interval as _get_experimental_interval,
)
from pdb_parser.io.distance_intervals import (
	get_precise_interval as _get_precise_interval,
)
from pdb_parser.io.distance_intervals import (
	get_xray_ca_centered_interval as _get_xray_ca_centered_interval,
)
from pdb_parser.io.nmr_constraints import nmr_distance_constraints
from pdb_parser.io.planar_constraints import planar_peptide_distances
from pdb_parser.io.topology import (
	AtomKey,
	AtomPair,
	TopologyData,
	build_covalent_and_planar_topology,
	covalent_and_planar_distances,
)
from pdb_parser.io.vdw_constraints import (
	build_vdw_torsion_pairs,
	is_backbone_dihedral_path,
	vdw_distances,
)
from pdb_parser.io.xray_constraints import xray_ca_distance_constraints

__all__ = [
	"AtomKey",
	"AtomPair",
	"TopologyData",
	"BackboneAtomRecord",
	"BackboneResidueRecord",
	"BackboneResidueMap",
	"DistanceRecord",
	"MergedDistanceRecord",
	"load_filtered_atoms_table",
	"build_covalent_and_planar_topology",
	"covalent_and_planar_distances",
	"build_vdw_torsion_pairs",
	"is_backbone_dihedral_path",
	"vdw_distances",
	"planar_peptide_distances",
	"nmr_distance_constraints",
	"xray_ca_distance_constraints",
	"talos_n_like",
	"wrap_angle_deg",
	"full_range_angle",
	"circular_distance_deg",
	"angle_belongs_to_interval",
	"torsion_interval_distance_bounds",
	"build_backbone_residue_map",
	"get_backbone_atom",
	"omega_interval_to_distance_record",
	"phi_interval_to_distance_record",
	"psi_interval_to_distance_record",
	"backbone_angular_interval_to_distance_interval",
	"merge_distance_constraint_files",
	"tightening_upper_bounds",
	"_get_centered_interval",
	"_get_experimental_interval",
	"_get_precise_interval",
	"_get_xray_ca_centered_interval",
]

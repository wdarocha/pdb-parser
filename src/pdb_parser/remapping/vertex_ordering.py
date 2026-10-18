"""Stable compatibility facade for legacy vertex-ordering imports."""

from pdb_parser.messages import print_skip

from .clique_patterns import (
	build_available_atoms,
	build_clique_pairs,
	build_ddgp_pattern_entries,
	get_atoms_of_residue,
	get_clique_pattern_ids,
	map_clique_type,
	resolve_atom_name,
)
from .dataframes import remap_distance_dataframe, remap_structure_dataframe
from .first_residue import build_first_residue_pattern
from .layout import (
    detect_remap_atom_layout,
    ensure_backbone_plus_hydrogens_dataframe,
    get_missing_backbone_atoms,
    get_missing_backbone_plus_hydrogen_atoms,
    validate_backbone_plus_hydrogens_residue,
)
from .ordering import (
    choose_first_existing_atom,
    first_residue_order,
    get_backbone_residue_numeric_order,
    get_first_residue_ha_atom,
    get_ha_group,
    get_internal_residue_numeric_order,
    get_internal_residue_ordering_template,
    get_nonterminal_h_group,
    get_numeric_atom_order_from_named_order,
    resolve_atom_name_choice,
)
from .torsions import (
	build_atom_clique_index_matrix,
	build_atom_lookup,
	build_backbone_torsion_matrix,
	build_interval_lookup,
	get_interval_dl_du,
	update_torsion_matrix,
)

__all__ = [
	"build_atom_clique_index_matrix",
	"build_atom_lookup",
	"build_backbone_torsion_matrix",
	"build_available_atoms",
	"build_clique_pairs",
	"build_ddgp_pattern_entries",
	"build_first_residue_pattern",
	"build_interval_lookup",
    "choose_first_existing_atom",
    "detect_remap_atom_layout",
    "first_residue_order",
	"get_backbone_residue_numeric_order",
	"get_atoms_of_residue",
	"get_clique_pattern_ids",
    "get_first_residue_ha_atom",
	"get_ha_group",
	"get_interval_dl_du",
    "get_internal_residue_numeric_order",
    "get_internal_residue_ordering_template",
    "get_nonterminal_h_group",
    "get_missing_backbone_atoms",
    "get_missing_backbone_plus_hydrogen_atoms",
    "get_numeric_atom_order_from_named_order",
	"map_clique_type",
	"resolve_atom_name",
	"resolve_atom_name_choice",
    "remap_distance_dataframe",
    "remap_structure_dataframe",
    "ensure_backbone_plus_hydrogens_dataframe",
	"validate_backbone_plus_hydrogens_residue",
	"update_torsion_matrix",
    "print_skip",
]

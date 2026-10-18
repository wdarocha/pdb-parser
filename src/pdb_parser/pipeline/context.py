"""Parser configuration and execution-context construction."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pdb_parser.io import ensure_model_chain_ready, resolve_effective_atom_selection
from pdb_parser.messages import print_warning
from pdb_parser.output_layout import (
    filtered_atoms_file,
    normalize_layout_components,
    parser_angular_file,
    parser_distance_file,
    parser_metadata_file,
    parser_structure_file,
    pdb_output_dir,
)
from pdb_parser.utils import ensure_dir, load_or_create_seed_bundle

from .contracts import (
    ParserContext,
    ParserExecutionContext,
    ParserParams,
    PathMap,
    SeedBundle,
)
from .validation import (
    _parse_float_or_default,
    _parse_float_param,
    _preflight_parser_params,
    _require_positive_param,
    _validate_noe_parameters,
    required_seed_fields_for_execution,
    resolve_effective_angular_width,
    resolve_effective_backbone_torsion_percentage,
)


def _prepare_parser_output_dirs(out_dir: str | Path, pdb_id: str) -> tuple[Path, Path]:
    """Create and return the per-PDB output and temporary directories."""
    out_dir_i = pdb_output_dir(out_dir, pdb_id)
    ensure_dir(out_dir_i)

    tmp_dir = out_dir_i / "tmp"
    ensure_dir(tmp_dir)

    return out_dir_i, tmp_dir


def _prepare_parser_context(
    params: ParserParams,
    pdb_data_dir: str | Path,
    pdb_id: str,
) -> ParserContext:
    """Resolve the validated parser context for one PDB/model/chain combination."""
    pdb_id, _, _, _ = normalize_layout_components(pdb_id, 1, "A")
    pdb_file_name = Path(pdb_data_dir) / f"{pdb_id}.pdb"
    parsed_params = _preflight_parser_params(params)
    requested_model = int(parsed_params["model_number"])
    chosen_chain = str(parsed_params["chain_id"])
    atom_selection = str(parsed_params["atom_selection"])
    distance_constraints = str(parsed_params["distance_constraints"])
    protein_chain, chosen_model, method = ensure_model_chain_ready(
        pdb_file_name,
        requested_model,
        chosen_chain,
        method=str(parsed_params["method"]),
    )
    if method == "nmr" and distance_constraints == "interval_centered":
        # These are NMR-specific even though X-ray also has a centred mode.
        parsed_params["epsilon_short"] = _parse_float_param(params, "epsilon_short")
        parsed_params["vdwr_hh"] = _parse_float_param(params, "vdwr_hh")
        _require_positive_param(parsed_params, "epsilon_short")
        _require_positive_param(parsed_params, "vdwr_hh")
    if method == "nmr" and distance_constraints == "interval_experimental":
        for key in ("noe_strong", "noe_medium", "noe_weak", "vdwr_hh"):
            parsed_params[key] = _parse_float_param(params, key)
        _validate_noe_parameters(parsed_params)
    if method == "xray" and atom_selection == "backbone_plus_hydrogens":
        atom_selection = resolve_effective_atom_selection(method, atom_selection)
        parsed_params["atom_selection"] = atom_selection
    # filtered_atoms writes no constraints: only the atom-selection conversion concerns it.
    writes_constraints = str(parsed_params.get("output_mode", "parser_outputs")) == "parser_outputs"
    if method == "xray" and distance_constraints == "interval_experimental":
        if writes_constraints:
            print_warning(
                "distance_constraints='interval_experimental' is NMR-specific; "
                "using 'interval_centered' for X-ray CA--CA constraints."
            )
        distance_constraints = "interval_centered"
        parsed_params["distance_constraints"] = distance_constraints
        parsed_params["epsilon_long"] = _parse_float_or_default(
            params, "epsilon_long", 2.0
        )
        parsed_params["max_distance"] = _parse_float_or_default(
            params, "max_distance", 5.0
        )
        _require_positive_param(parsed_params, "epsilon_long")
        _require_positive_param(parsed_params, "max_distance")
    angular_width = resolve_effective_angular_width(
        distance_constraints,
        float(parsed_params["torsion_angle_width"]),
    )
    backbone_torsion_percentage = resolve_effective_backbone_torsion_percentage(
        distance_constraints,
        float(parsed_params["percentage_backbone_torsion_angles"]),
    )
    hh_distance_percentage = float(parsed_params["percentage_hh_distance_constraints"])
    if method == "xray" and hh_distance_percentage < 100.0 and writes_constraints:
        print_warning(
            "percentage_hh_distance_constraints is NMR-specific; "
            "ignored for X-ray CA--CA constraints."
        )
    required_seed_fields = required_seed_fields_for_execution(
        distance_constraints,
        angular_width,
        backbone_torsion_percentage,
        len(protein_chain.residues),
        method=method,
        hh_distance_percentage=hh_distance_percentage,
    )

    return {
        "pdb_file_name": pdb_file_name,
        "parsed_params": parsed_params,
        "chosen_model": chosen_model,
        "chosen_chain": chosen_chain,
        "atom_selection": atom_selection,
        "distance_constraints": distance_constraints,
        "method": method,
        "output_mode": str(parsed_params["output_mode"]),
        "protein_chain": protein_chain,
        "angular_width": angular_width,
        "backbone_torsion_percentage": backbone_torsion_percentage,
        "required_seed_fields": required_seed_fields,
    }


def _load_parser_seed_generators(
    seed_dir: str | Path,
    pdb_id: str,
    chosen_model: int,
    chosen_chain: str,
    required_seed_fields: set[str],
) -> tuple[
    SeedBundle,
    np.random.Generator | None,
    np.random.Generator | None,
    np.random.Generator | None,
    np.random.Generator | None,
]:
    """Load or create the per-PDB seed bundle and the RNGs derived from it."""
    seed_bundle = load_or_create_seed_bundle(
        seed_dir,
        pdb_id,
        chosen_model,
        chosen_chain,
        required_seed_fields,
    )
    return (seed_bundle, *_build_seed_generators(seed_bundle))


def _build_seed_generators(seed_bundle: SeedBundle) -> tuple[
    np.random.Generator | None,
    np.random.Generator | None,
    np.random.Generator | None,
    np.random.Generator | None,
]:
    """Build the RNGs used by the parser from the resolved seed bundle."""
    rng_interval_centered = (
        np.random.default_rng(seed_bundle["seed_interval_centered_distances"])
        if seed_bundle["seed_interval_centered_distances"] is not None
        else None
    )
    rng_talos_angles = (
        np.random.default_rng(seed_bundle["seed_talos_angle_centers"])
        if seed_bundle["seed_talos_angle_centers"] is not None
        else None
    )
    rng_phi_psi_mask = (
        np.random.default_rng(seed_bundle["seed_phi_psi_mask"])
        if seed_bundle["seed_phi_psi_mask"] is not None
        else None
    )
    rng_hh_distance_subset = (
        np.random.default_rng(seed_bundle["seed_hh_distance_subset"])
        if seed_bundle.get("seed_hh_distance_subset") is not None
        else None
    )

    return rng_interval_centered, rng_talos_angles, rng_phi_psi_mask, rng_hh_distance_subset


def _build_parser_output_files(
    out_dir_i: Path,
    tmp_dir: Path,
    pdb_id: str,
    chosen_model: int,
    chosen_chain: str,
) -> PathMap:
    """Return the canonical parser output paths for one structure."""
    return {
        "tsv_structure": tmp_dir / "filtered_atoms.dat",
        "filtered_atoms": filtered_atoms_file(
            out_dir_i.parent, pdb_id, chosen_model, chosen_chain
        ),
        "metadata": parser_metadata_file(
            out_dir_i.parent, pdb_id, chosen_model, chosen_chain
        ),
        "x": parser_structure_file(
            out_dir_i.parent, pdb_id, chosen_model, chosen_chain
        ),
        "a_tmp": tmp_dir / "angular_constraints_5.dat",
        "a": parser_angular_file(out_dir_i.parent, pdb_id, chosen_model, chosen_chain),
        "i": parser_distance_file(out_dir_i.parent, pdb_id, chosen_model, chosen_chain),
        "dc_1": tmp_dir / "distance_constraints_1.dat",
        "dc_2": tmp_dir / "distance_constraints_2.dat",
        "dc_3": tmp_dir / "distance_constraints_3.dat",
        "dc_4": tmp_dir / "distance_constraints_4.dat",
        "dc_5": tmp_dir / "distance_constraints_5.dat",
    }


def _prepare_parser_execution_context(
    params: ParserParams,
    pdb_data_dir: str | Path,
    seed_dir: str | Path,
    out_dir: str | Path,
    pdb_id: str,
) -> ParserExecutionContext:
    """Resolve the full execution context for one parser job."""
    pdb_id, _, _, _ = normalize_layout_components(pdb_id, 1, "A")
    parser_context = _prepare_parser_context(params, pdb_data_dir, pdb_id)
    chosen_model = parser_context["chosen_model"]
    chosen_chain = parser_context["chosen_chain"]
    output_mode = parser_context["output_mode"]
    required_seed_fields = parser_context["required_seed_fields"]
    out_dir_i, tmp_dir = _prepare_parser_output_dirs(out_dir, pdb_id)
    output_files = _build_parser_output_files(
        out_dir_i,
        tmp_dir,
        pdb_id,
        chosen_model,
        chosen_chain,
    )
    tsv_structure_file = (
        output_files["filtered_atoms"]
        if output_mode == "filtered_atoms"
        else output_files["tsv_structure"]
    )
    if output_mode == "filtered_atoms":
        rng_interval_centered = None
        rng_talos_angles = None
        rng_phi_psi_mask = None
        rng_hh_distance_subset = None
    else:
        _, rng_interval_centered, rng_talos_angles, rng_phi_psi_mask, rng_hh_distance_subset = (
            _load_parser_seed_generators(
                seed_dir,
                pdb_id,
                chosen_model,
                chosen_chain,
                required_seed_fields,
            )
        )

    return {
        **parser_context,
        "tmp_dir": tmp_dir,
        "output_files": output_files,
        "tsv_structure_file": tsv_structure_file,
        "rng_interval_centered": rng_interval_centered,
        "rng_talos_angles": rng_talos_angles,
        "rng_phi_psi_mask": rng_phi_psi_mask,
        "rng_hh_distance_subset": rng_hh_distance_subset,
    }


def prepare_parser_context(
    params: ParserParams,
    pdb_data_dir: str | Path,
    pdb_id: str,
) -> ParserContext:
    """Public compatibility entry point for parser context construction."""
    return _prepare_parser_context(params, pdb_data_dir, pdb_id)


def prepare_parser_execution_context(
    params: ParserParams,
    pdb_data_dir: str | Path,
    seed_dir: str | Path,
    out_dir: str | Path,
    pdb_id: str,
) -> ParserExecutionContext:
    """Public compatibility entry point for parser execution context construction."""
    return _prepare_parser_execution_context(
        params, pdb_data_dir, seed_dir, out_dir, pdb_id
    )

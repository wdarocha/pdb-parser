"""Remapping-stage context, transformations, and artifact writing."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np

from pdb_parser.io import (
	read_validated_positional_file,
	save_cliques_from_matrix_T,
	save_coordinates_from_df_structure,
	save_distances_from_df_structure,
	validate_positional_cross_references,
	validate_positional_schema,
)
from pdb_parser.io.sequence import validate_residue_sequence
from pdb_parser.messages import print_warning
from pdb_parser.output_layout import (
    parser_angular_file,
    parser_distance_file,
    parser_metadata_file,
    parser_structure_file,
    pdb_output_dir,
    remapped_clique_file,
    remapped_distance_file,
    remapped_output_dir,
    remapped_structure_file,
)
from pdb_parser.public_config import (
    DDGP_ORDER_TYPE,
    DMDGP_ORDER_ID,
    DMDGP_ORDER_TYPE,
    REDMDGP_ORDER_TYPE,
    parse_public_order_choice,
    parse_public_order_type,
)
from pdb_parser.remapping import (
    build_atom_clique_index_matrix,
    build_available_atoms,
    build_backbone_torsion_matrix,
    build_ddgp_pattern_entries,
    build_first_residue_pattern,
    detect_remap_atom_layout,
    first_residue_order,
    get_backbone_residue_numeric_order,
    get_internal_residue_numeric_order,
    remap_distance_dataframe,
    remap_structure_dataframe,
    validate_backbone_plus_hydrogens_residue,
)
from pdb_parser.remapping.redmdgp import (
    build_distance_lookup,
    build_redmdgp_atom_sequence,
    build_redmdgp_clique_matrix,
    build_redmdgp_distance_dataframe,
    validate_redmdgp_sequence,
)
from pdb_parser.utils import ensure_dir

from .contracts import (
    AngularFrame,
    AtomOrder,
    CliqueEntry,
    DistanceFrame,
    OrderVector,
    ParserParams,
    RemapContext,
    RemapExecutionContext,
    StructureFrame,
)


def _validate_remap_metadata(
    out_dir: str | Path,
    pdb_id: str,
    model_number: int,
    chain_id: str,
) -> None:
    """Reject remapping layouts not defined for the parser method/selection pair."""
    metadata_path = parser_metadata_file(out_dir, pdb_id, model_number, chain_id)
    if not metadata_path.exists():
        return
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    allowed = {
        ("nmr", "backbone"),
        ("nmr", "backbone_plus_hydrogens"),
        ("xray", "backbone"),
    }
    pair = (metadata.get("method"), metadata.get("atom_selection"))
    if metadata.get("output_mode") != "parser_outputs":
        raise ValueError("Remapping requires parser metadata with output_mode: parser_outputs.")
    if pair not in allowed:
        raise ValueError(
            "Remapping is not defined for parser method/atom_selection "
            f"{pair!r}. Allowed combinations are: {sorted(allowed)}."
        )

def _build_remap_input_files(
    out_dir: str | Path,
    pdb_id: str,
    chosen_model: int,
    chosen_chain: str,
) -> tuple[Path, Path, Path, Path]:
    """Return the remapping directory and its three input files."""
    return (
        pdb_output_dir(out_dir, pdb_id),
        parser_structure_file(out_dir, pdb_id, chosen_model, chosen_chain),
        parser_angular_file(out_dir, pdb_id, chosen_model, chosen_chain),
        parser_distance_file(out_dir, pdb_id, chosen_model, chosen_chain),
    )


def _prepare_remap_context(
    params: ParserParams,
    out_dir: str | Path,
    pdb_id: str,
) -> RemapContext:
    """Resolve the basic context used by the remapping stage."""
    chosen_model = int(params["model_number"])
    chosen_chain = str(params["chain_id"])
    out_dir_pdb_id, xfile, afile, ifile = _build_remap_input_files(
        out_dir,
        pdb_id,
        chosen_model,
        chosen_chain,
    )

    return {
        "chosen_model": chosen_model,
        "chosen_chain": chosen_chain,
        "out_dir_pdb_id": out_dir_pdb_id,
        "xfile": xfile,
        "afile": afile,
        "ifile": ifile,
    }


def _build_remap_output_files(
    out_dir_pdb_id: Path,
    pdb_id: str,
    chosen_model: int,
    chosen_chain: str,
    chosen_order_id: int | str,
    atom_layout: str,
    order_type: str = DDGP_ORDER_TYPE,
) -> tuple[Path, Path, Path, Path]:
    """Create and return the remapped output directory and files for one structure."""
    out_dir_pdb_id_remapped = remapped_output_dir(out_dir_pdb_id.parent, pdb_id)
    ensure_dir(out_dir_pdb_id_remapped)

    return (
        out_dir_pdb_id_remapped,
        remapped_structure_file(
            out_dir_pdb_id.parent,
            pdb_id,
            chosen_model,
            chosen_chain,
            chosen_order_id,
            atom_layout,
            order_type,
        ),
        remapped_distance_file(
            out_dir_pdb_id.parent,
            pdb_id,
            chosen_model,
            chosen_chain,
            chosen_order_id,
            atom_layout,
            order_type,
        ),
        remapped_clique_file(
            out_dir_pdb_id.parent,
            pdb_id,
            chosen_model,
            chosen_chain,
            chosen_order_id,
            atom_layout,
            order_type,
        ),
    )


def _load_remap_input_dataframes(
    xfile: Path,
    afile: Path,
    ifile: Path,
) -> tuple[StructureFrame, DistanceFrame, AngularFrame]:
    """Read the structure, angular and distance files for the remapping stage."""
    df_X = read_validated_positional_file(xfile, "X")
    # X is the authority for the residue/atom universe.  Validate it before
    # opening the dependent tables so a gapped PDB cannot reach any writer.
    _validate_remap_x_table(df_X)
    df_I = read_validated_positional_file(ifile, "I")
    df_A = read_validated_positional_file(afile, "A")
    _validate_remap_tables(df_X, df_I, df_A)
    return df_X, df_I, df_A


def _validate_remap_x_table(df_X: StructureFrame) -> list[int]:
	"""Validate X_* continuity and return its ordered residue IDs."""
	validate_positional_schema(df_X, "X", "<remapping X data>")
	ordered_residues = list(dict.fromkeys(df_X.iloc[:, 2].astype(int).tolist()))
	try:
		return validate_residue_sequence(ordered_residues, context="remapping X file")
	except ValueError as exc:
		raise ValueError(
			"Remap cannot proceed because residue gaps were found in X_*: "
			f"{exc}"
		) from exc


def _validate_remap_tables(df_X: StructureFrame, df_I: DistanceFrame, df_A: AngularFrame) -> None:
	"""Validate cross-file residue and atom references before remapping."""
	_validate_remap_x_table(df_X)
	validate_positional_cross_references(df_X, df_A, df_I, "<remapping X/A/I data>")


def _resolve_remap_output_order_id(params: ParserParams) -> int | str:
    """Resolve the public order choice used in remapped output filenames."""
    if parse_public_order_type(params) == DMDGP_ORDER_TYPE and not str(params.get("order_id", "")).strip():
        return DMDGP_ORDER_ID
    return parse_public_order_choice(str(params["order_id"]))


def _prepare_remap_execution_context(
    params: ParserParams,
    out_dir: str | Path,
    pdb_id: str,
) -> RemapExecutionContext:
    """Resolve the full execution context for one remap job."""
    remap_context = _prepare_remap_context(params, out_dir, pdb_id)
    chosen_model = remap_context["chosen_model"]
    chosen_chain = remap_context["chosen_chain"]
    _validate_remap_metadata(out_dir, pdb_id, chosen_model, chosen_chain)
    xfile = remap_context["xfile"]
    afile = remap_context["afile"]
    ifile = remap_context["ifile"]
    df_X, df_I, df_A = _load_remap_input_dataframes(xfile, afile, ifile)
    res_ids, nres = _extract_residue_ids_and_count(df_X)
    atom_layout = detect_remap_atom_layout(df_X, res_ids)

    return {
        **remap_context,
        "df_X": df_X,
        "df_I": df_I,
        "df_A": df_A,
        "res_ids": res_ids,
        "nres": nres,
        "chosen_order_id": _resolve_remap_output_order_id(params),
        "order_type": parse_public_order_type(params),
        "atom_layout": atom_layout,
    }


def prepare_remap_execution_context(
    params: ParserParams,
    out_dir: str | Path,
    pdb_id: str,
) -> RemapExecutionContext:
    """Public wrapper for resolving one remap execution context."""
    return _prepare_remap_execution_context(params, out_dir, pdb_id)


def _build_residue_atom_order(
    df_X: StructureFrame,
    res_ids: list[int],
    pdb_id: str,
    chosen_model: int,
    chosen_chain: str,
    ddgp_order_vec: OrderVector,
    atom_layout: str,
    residues_validated: bool = False,
) -> tuple[AtomOrder, int]:
    """Build the residue-by-residue atom order and return it with the skip flag.

    ``residues_validated`` skips the per-residue atom check already passed by an
    earlier order of the same structure (see ``shared_cache``).
    """
    validate_positional_schema(df_X, "X", "<in-memory remapping X data>")
    new_order: AtomOrder = []
    skip_flag = 0

    validate_residue_sequence(res_ids, context="remapping X file")
    nres = len(res_ids)

    for k in range(nres):
        resnum_i = res_ids[k]
        df_X_i = df_X[df_X[2] == resnum_i]

        if atom_layout == "backbone_plus_hydrogens" and not residues_validated:
            skip_flag = validate_backbone_plus_hydrogens_residue(
                df_X_i,
                resnum_i,
                pdb_id,
                chosen_model,
                chosen_chain,
                res_ids[0],
            )

        if skip_flag == 1:
            break

        if atom_layout == "backbone":
            new_order.extend(get_backbone_residue_numeric_order(df_X_i))
        elif k == 0:
            new_order.extend(first_residue_order(df_X_i, res_ids[k]))
        else:
            new_order.extend(
                get_internal_residue_numeric_order(df_X_i, ddgp_order_vec[k])
            )

    return new_order, skip_flag


def _build_remapped_clique_matrix(
    df_Xremapped: StructureFrame,
    df_A: AngularFrame,
    df_Iremapped: DistanceFrame,
    res_ids: list[int],
    nres: int,
    ddgp_order_vec: OrderVector,
    atom_layout: str,
) -> np.ndarray:
    """Build the clique matrix T for the remapped structure."""
    validate_positional_schema(df_Xremapped, "X", "<in-memory remapped X data>")
    validate_positional_schema(df_A, "A", "<in-memory remapping A data>")
    validate_positional_schema(df_Iremapped, "I", "<in-memory remapped I data>")
    if atom_layout == "backbone":
        return build_backbone_torsion_matrix(df_Xremapped, df_A)

    atom_cliques: list[CliqueEntry] = []
    available_atoms = build_available_atoms(df_Xremapped)

    for k in range(nres):
        if k == 0:
            atom_cliques.extend(
                build_first_residue_pattern(available_atoms, res_ids[k])
            )
        else:
            atom_cliques.extend(
                build_ddgp_pattern_entries(
                    ddgp_order_vec[k], res_ids[k], available_atoms
                )
            )

    return build_atom_clique_index_matrix(
        atom_cliques,
        df_Xremapped,
        df_A,
        df_Iremapped,
    )


def _extract_residue_ids_and_count(df_X: StructureFrame) -> tuple[list[int], int]:
    """Return the ordered residue ids and the current residue-count convention."""
    validate_positional_schema(df_X, "X", "<in-memory remapping X data>")
    res_ids = validate_residue_sequence(
        list(dict.fromkeys(df_X.iloc[:, 2].astype(int).tolist())),
        context="remapping X file",
    )
    nres = len(res_ids)
    return res_ids, nres


def _remap_structure_and_distance_dataframes(
    df_X: StructureFrame,
    df_I: DistanceFrame,
    new_order: AtomOrder,
) -> tuple[StructureFrame, DistanceFrame]:
    """Remap the structure and distance dataframes according to the new atom order."""
    validate_positional_schema(df_X, "X", "<in-memory remapping X data>")
    validate_positional_schema(df_I, "I", "<in-memory remapping I data>")
    atom_ids_new_order = [atom_id for atom_id, _ in new_order]
    return (
        remap_structure_dataframe(df_X, atom_ids_new_order),
        remap_distance_dataframe(df_I, atom_ids_new_order),
    )


def _save_remapped_outputs(
    df_Iremapped: DistanceFrame,
    df_Xremapped: StructureFrame,
    T: np.ndarray,
    out_dir_pdb_id: Path,
    pdb_id: str,
    chosen_model: int,
    chosen_chain: str,
    chosen_order_id: int | str,
    atom_layout: str,
    order_type: str = DDGP_ORDER_TYPE,
) -> None:
	"""Publish remapped X, I and T files as one recoverable transaction."""
	output_dir, xfile, ifile, tfile = _build_remap_output_files(
		out_dir_pdb_id,
        pdb_id,
        chosen_model,
        chosen_chain,
        chosen_order_id,
		atom_layout,
		order_type,
	)
	final_files = (xfile, ifile, tfile)
	staging_dir = Path(tempfile.mkdtemp(prefix=".remap-staging-", dir=output_dir))
	staged_files = tuple(staging_dir / path.name for path in final_files)
	backups = {
		path: path.read_bytes() for path in final_files if path.exists() and path.is_file()
	}
	try:
		save_distances_from_df_structure(df_Iremapped, staged_files[1])
		save_coordinates_from_df_structure(df_Xremapped, staged_files[0])
		save_cliques_from_matrix_T(T, staged_files[2])
		for staged_file, final_file in zip(staged_files, final_files, strict=True):
			os.replace(staged_file, final_file)
	except Exception:
		original = __import__("sys").exc_info()[1]
		cleanup_errors: list[Exception] = []
		for final_file in final_files:
			try:
				if final_file in backups:
					final_file.write_bytes(backups[final_file])
				elif final_file.exists() and final_file.is_file():
					final_file.unlink()
			except Exception as exc:
				cleanup_errors.append(exc)
		if cleanup_errors and original is not None:
			diagnostic = "Rollback cleanup failed: " + "; ".join(
				str(error) for error in cleanup_errors
			)
			if hasattr(original, "add_note"):
				original.add_note(diagnostic)
			else:
				print_warning(diagnostic)
		raise
	finally:
		shutil.rmtree(staging_dir, ignore_errors=True)


def _write_remapped_outputs_from_vector(
    context: RemapExecutionContext,
    pdb_id: str,
    ddgp_order_vec: OrderVector,
) -> int:
    """Build remapped dataframes and write the final outputs for one structure."""
    # Backbone outputs keep the fixed DMDGP order whatever the order_type.
    if (
        context.get("order_type", DDGP_ORDER_TYPE) == REDMDGP_ORDER_TYPE
        and context["atom_layout"] == "backbone_plus_hydrogens"
    ):
        return _write_redmdgp_outputs(context, pdb_id)
    cache = context.get("shared_cache")
    new_order, skip_flag = _build_residue_atom_order(
        context["df_X"],
        context["res_ids"],
        pdb_id,
        context["chosen_model"],
        context["chosen_chain"],
        ddgp_order_vec,
        context["atom_layout"],
        residues_validated=bool(cache and cache.get("residues_validated")),
    )
    if skip_flag == 0 and cache is not None:
        cache["residues_validated"] = True
    if skip_flag == 1:
        return skip_flag

    df_Xremapped, df_Iremapped = _remap_structure_and_distance_dataframes(
        context["df_X"],
        context["df_I"],
        new_order,
    )
    T = _build_remapped_clique_matrix(
        df_Xremapped,
        context["df_A"],
        df_Iremapped,
        context["res_ids"],
        context["nres"],
        ddgp_order_vec,
        context["atom_layout"],
    )
    _save_remapped_outputs(
        df_Iremapped,
        df_Xremapped,
        T,
        context["out_dir_pdb_id"],
        pdb_id,
        context["chosen_model"],
        context["chosen_chain"],
        context["chosen_order_id"],
        context["atom_layout"],
    )

    return skip_flag



def _write_redmdgp_outputs(context: RemapExecutionContext, pdb_id: str) -> int:
    """Write the iDMDGP re-order X/I/T files (``order_type: redmdgp``) for one
    ``backbone_plus_hydrogens`` structure."""
    df_X, df_I, df_A, res_ids = context["df_X"], context["df_I"], context["df_A"], context["res_ids"]
    # Orders of one structure share the per-residue check and the parser distance lookup.
    cache = context.get("shared_cache")
    if not (cache and cache.get("residues_validated")):
        for residue_id in res_ids:
            if validate_backbone_plus_hydrogens_residue(
                df_X[df_X[2] == residue_id],
                residue_id,
                pdb_id,
                context["chosen_model"],
                context["chosen_chain"],
                res_ids[0],
            ):
                return 1
        if cache is not None:
            cache["residues_validated"] = True

    sequence = build_redmdgp_atom_sequence(df_X, res_ids, int(context["chosen_order_id"]))
    distances = cache.get("parser_distance_lookup") if cache is not None else None
    if distances is None:
        distances = build_distance_lookup(df_I)
        if cache is not None:
            cache["parser_distance_lookup"] = distances
    validate_redmdgp_sequence(sequence, df_X, distances)
    _save_remapped_outputs(
        build_redmdgp_distance_dataframe(sequence, df_X, df_I, distances),
        remap_structure_dataframe(df_X, sequence),
        build_redmdgp_clique_matrix(sequence, df_X, df_A, distances),
        context["out_dir_pdb_id"],
        pdb_id,
        context["chosen_model"],
        context["chosen_chain"],
        context["chosen_order_id"],
        context["atom_layout"],
        REDMDGP_ORDER_TYPE,
    )
    return 0


def prepare_remap_context(
    params: ParserParams,
    out_dir: str | Path,
    pdb_id: str,
) -> RemapExecutionContext:
    """Public compatibility entry point for remap context construction."""
    return _prepare_remap_execution_context(params, out_dir, pdb_id)


def write_remapped_artifacts(
    params: ParserParams,
    out_dir: str | Path,
    pdb_id: str,
    ddgp_order_vec: OrderVector,
    context: RemapExecutionContext | None = None,
) -> int:
    """Create remapped X/I/T artifacts, returning the historical skip flag."""
    execution_context = context or _prepare_remap_execution_context(
        params, out_dir, pdb_id
    )
    return _write_remapped_outputs_from_context(
        execution_context,
        pdb_id,
        ddgp_order_vec,
    )


def _write_remapped_outputs_from_context(
    context: RemapExecutionContext,
    pdb_id: str,
    ddgp_order_vec: OrderVector,
) -> int:
    """Write outputs using one typed execution context."""
    return _write_remapped_outputs_from_vector(context, pdb_id, ddgp_order_vec)

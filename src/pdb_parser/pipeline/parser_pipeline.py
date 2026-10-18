"""Compatibility facade for the parser and remapping pipeline.

Implementation lives in the specialized pipeline modules.
"""

from __future__ import annotations

from pathlib import Path

from .artifacts import (  # noqa: F401
    _cleanup_tmp_dir,
    _prepare_filtered_structure_file,
    _write_parser_metadata,
    _write_parser_outputs,
    cleanup_parser_artifacts,
    write_parser_artifacts,
)
from .context import (  # noqa: F401
    _prepare_parser_context,
    _prepare_parser_execution_context,
    prepare_parser_context,
    prepare_parser_execution_context,
)
from .contracts import OrderVector, ParserParams, RemapExecutionContext
from .remapping import (  # noqa: F401
    _prepare_remap_execution_context,
    _write_remapped_outputs_from_vector,
    prepare_remap_context,
    write_remapped_artifacts,
)
from .validation import _preflight_parser_params, validate_public_parser_params  # noqa: F401


def parser(
    params: ParserParams,
    pdb_data_dir: str | Path,
    seed_dir: str | Path,
    out_dir: str | Path,
    pdb_id: str,
    remove_tmp_dir: bool = False,
) -> None:
	"""Run the full parser pipeline for one PDB/model/chain combination."""
	context = prepare_parser_execution_context(params, pdb_data_dir, seed_dir, out_dir, pdb_id)
	write_parser_artifacts(context)
	cleanup_parser_artifacts(context, remove_tmp_dir)


def remap_instance(
    params: ParserParams,
    out_dir: str | Path,
    pdb_id: str,
    ddgp_order_vec: OrderVector,
    remap_context: RemapExecutionContext | None = None,
) -> int:
    """Remap parser artifacts and return the historical skip flag."""
    return write_remapped_artifacts(params, out_dir, pdb_id, ddgp_order_vec, remap_context)


def prepare_remap_execution_context(
    params: ParserParams,
    out_dir: str | Path,
    pdb_id: str,
) -> RemapExecutionContext:
    """Public compatibility entry point for remap execution contexts."""
    return prepare_remap_context(params, out_dir, pdb_id)


__all__ = [
    "parser",
    "prepare_remap_execution_context",
    "remap_instance",
    "validate_public_parser_params",
]

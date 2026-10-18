"""Parser artifact generation and cleanup."""

from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup

from pdb_parser.io import (
    TopologyData,
    backbone_angular_interval_to_distance_interval,
    build_covalent_and_planar_topology,
    convert_tsv_angular_to_pdb_format,
    convert_tsv_structure_to_pdb_format,
    covalent_and_planar_distances,
    merge_distance_constraint_tables,
    nmr_distance_constraints,
    planar_peptide_distances,
    save_filtered_atoms,
    talos_n_like,
    vdw_distances,
    xray_ca_distance_constraints,
)
from pdb_parser.messages import print_warning
from pdb_parser.public_config import VALID_DISTANCE_MODES

from .contracts import ParserExecutionContext, ParserParams


def _prepare_filtered_structure_file(
    protein_chain: AtomGroup,
    atom_selection: str,
    tsv_structure_file: Path,
    method: str,
) -> None:
    """Write the filtered structure TSV used by downstream parser stages."""
    save_filtered_atoms(
        tsv_structure_file,
        protein_chain,
        atom_selection,
        method=method,
    )


def _cleanup_tmp_dir(tmp_dir: Path, remove_tmp_dir: bool) -> None:
    """Remove the temporary parser directory when cleanup is enabled."""
    if remove_tmp_dir and tmp_dir.exists():
        shutil.rmtree(tmp_dir)


def _write_parser_metadata(
    output_file: Path,
    method: str,
    atom_selection: str,
    output_mode: str,
) -> None:
    """Persist the effective parser choices consumed by the remapping stage."""
    output_file.write_text(
        json.dumps(
            {
                "method": method,
                "atom_selection": atom_selection,
                "output_mode": output_mode,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _write_nmr_distance_constraints(
    tsv_structure_file: Path,
    output_distance_file: Path,
    distance_constraints: str,
    atom_selection: str,
    parsed_params: ParserParams,
    rng_interval_centered: np.random.Generator | None,
    rng_hh_distance_subset: np.random.Generator | None = None,
):
    """Write the NMR-derived distance-constraint file for the selected public mode; return its table."""
    subset = {
        "hh_distance_percentage": float(parsed_params.get("percentage_hh_distance_constraints", 100.0)),
        "rng_hh_distance_subset": rng_hh_distance_subset,
    }
    if distance_constraints == "interval_centered":
        return nmr_distance_constraints(
            tsv_structure_file,
            output_distance_file,
            distance_constraints,
            atom_selection,
            float(parsed_params["epsilon_short"]),
            float(parsed_params["epsilon_long"]),
            float(parsed_params["max_distance"]),
            None,
            None,
            None,
            float(parsed_params["vdwr_hh"]),
            rng_interval_centered,
            **subset,
        )

    if distance_constraints == "precise":
        return nmr_distance_constraints(
            tsv_structure_file,
            output_distance_file,
            distance_constraints,
            atom_selection,
            None,
            None,
            float(parsed_params["max_distance"]),
            None,
            None,
            None,
            0.0,
            None,
            **subset,
        )

    if distance_constraints == "interval_experimental":
        return nmr_distance_constraints(
            tsv_structure_file,
            output_distance_file,
            distance_constraints,
            atom_selection,
            None,
            None,
            None,
            float(parsed_params["noe_strong"]),
            float(parsed_params["noe_medium"]),
            float(parsed_params["noe_weak"]),
            float(parsed_params["vdwr_hh"]),
            None,
            **subset,
        )

    raise ValueError(
        f"Unsupported distance_constraints mode: {distance_constraints!r}. "
        f"Expected one of {sorted(VALID_DISTANCE_MODES)}."
    )


def _write_long_range_distance_constraints(
    tsv_structure_file: Path,
    output_distance_file: Path,
    method: str,
    distance_constraints: str,
    atom_selection: str,
    parsed_params: ParserParams,
    rng_interval_centered: np.random.Generator | None,
    rng_hh_distance_subset: np.random.Generator | None = None,
):
    """Dispatch long-range restraints according to the experimental method; return their table."""
    if method == "nmr":
        return _write_nmr_distance_constraints(
            tsv_structure_file,
            output_distance_file,
            distance_constraints,
            atom_selection,
            parsed_params,
            rng_interval_centered,
            rng_hh_distance_subset,
        )

    return xray_ca_distance_constraints(
        tsv_structure_file,
        output_distance_file,
        distance_constraints,
        (
            float(parsed_params["epsilon_long"])
            if distance_constraints == "interval_centered"
            else None
        ),
        float(parsed_params["max_distance"]),
        rng_interval_centered,
    )


def _write_parser_outputs(
    context: ParserExecutionContext,
    topology: TopologyData,
) -> None:
    """Run output-generation stages after the filtered structure is available."""
    convert_tsv_structure_to_pdb_format(
        context["tsv_structure_file"], context["output_files"]["x"]
    )

    # Each stage writes its tmp/ file and returns its full-precision table for the merge.
    tables = [covalent_and_planar_distances(
        context["tsv_structure_file"], topology, context["output_files"]["dc_1"]
    )]

    if context["parsed_params"]["vdw_constraints"] == "yes":
        tables.append(vdw_distances(
            context["tsv_structure_file"], topology, context["output_files"]["dc_2"]
        ))

    tables.append(planar_peptide_distances(
        context["tsv_structure_file"],
        context["output_files"]["dc_3"],
        context["atom_selection"],
    ))

    tables.append(_write_long_range_distance_constraints(
        context["tsv_structure_file"],
        context["output_files"]["dc_4"],
        context["method"],
        context["distance_constraints"],
        context["atom_selection"],
        context["parsed_params"],
        context["rng_interval_centered"],
        context.get("rng_hh_distance_subset"),
    ))

    talos_n_like(
        context["pdb_file_name"],
        context["chosen_model"],
        context["chosen_chain"],
        context["output_files"]["a_tmp"],
        0.0,
        context["angular_width"],
        context["angular_width"],
        context["backbone_torsion_percentage"],
        context["rng_talos_angles"],
        context["rng_phi_psi_mask"],
    )
    convert_tsv_angular_to_pdb_format(
        context["output_files"]["a_tmp"], context["output_files"]["a"]
    )

    tables.append(backbone_angular_interval_to_distance_interval(
        context["tsv_structure_file"],
        context["output_files"]["a_tmp"],
        context["output_files"]["dc_5"],
    ))

    merge_distance_constraint_tables(
        tables,
        context["output_files"]["i"],
        tighten_upper_bounds=context["parsed_params"]["tightening_upper_bounds"] == "yes",
    )


def _remove_staged_artifacts(
    paths: tuple[Path, ...], staging_dir: Path
) -> tuple[tuple[Path, Path], ...]:
    """Remove paths transactionally, retaining a restorable copy while doing so."""
    staged: list[tuple[Path, Path]] = []
    try:
        # Keep preparation and deletion in the same transaction.  In
        # particular, do not let an exception while copying/renaming a later
        # artifact escape before the already moved artifacts are restorable.
        for index, path in enumerate(paths):
            if not path.exists():
                continue
            backup = staging_dir / f"backup-{index}"
            if path.is_dir():
                shutil.copytree(path, backup)
            else:
                shutil.copy2(path, backup)
            path.rename(staging_dir / f"staged-{index}")
            staged.append((path, backup))

        for index, (_path, _backup) in enumerate(staged):
            staged_path = staging_dir / f"staged-{paths.index(_path)}"
            if staged_path.is_dir():
                shutil.rmtree(staged_path)
            elif staged_path.exists():
                staged_path.unlink()
    except Exception:
        _restore_staged_artifacts(tuple(staged), staging_dir)
        raise
    return tuple(staged)


def _restore_staged_artifacts(
    staged: tuple[tuple[Path, Path], ...], staging_dir: Path
) -> None:
    """Restore artifacts from cleanup backups without relying on deletion."""
    for index, (path, backup) in enumerate(staged):
        if path.exists():
            displaced = staging_dir / f"displaced-{index}"
            path.rename(displaced)
        if backup.is_dir():
            shutil.copytree(backup, path)
        elif backup.exists():
            shutil.copy2(backup, path)


def write_parser_artifacts(context: ParserExecutionContext) -> None:
    """Generate all artifacts described by an execution context transactionally.

    The individual writers intentionally retain their historical paths and
    formats.  This boundary is responsible for making the complete parser
    stage recoverable when a later writer fails: existing published files are
    restored and temporary/intermediate files are removed.
    """
    output_files = context["output_files"]
    published_paths = tuple(
        output_files[key]
        for key in ("filtered_atoms", "metadata", "x", "a", "i")
    )
    remapped_path = output_files["metadata"].parent / "remapped"
    rollback_dir = Path(tempfile.mkdtemp(
        prefix=".parser-rollback-", dir=str(output_files["metadata"].parent)
    ))
    backups = {
        path: path.read_bytes()
        for path in published_paths
        if path.exists() and path.is_file()
    }
    staged: tuple[tuple[Path, Path], ...] = ()

    try:
        _prepare_filtered_structure_file(
            context["protein_chain"],
            context["atom_selection"],
            context["tsv_structure_file"],
            context["method"],
        )
        _write_parser_metadata(
            output_files["metadata"],
            context["method"],
            context["atom_selection"],
            context["output_mode"],
        )
        if context["output_mode"] == "filtered_atoms":
            # Remove prior parser/remap products only after the new filtered
            # artifact and its metadata have been written successfully.
            stale_paths = tuple(output_files[key] for key in ("x", "a", "i"))
            stale_paths += (remapped_path,)
            staged = _remove_staged_artifacts(stale_paths, rollback_dir)
            return
        topology = build_covalent_and_planar_topology(context["tsv_structure_file"])
        _write_parser_outputs(context, topology)
        # A parser-output generation supersedes the persistent filtered-only
        # artifact.  Do this last so a failed generation can restore it.
        staged = _remove_staged_artifacts(
            (output_files["filtered_atoms"],), rollback_dir
        )
    except Exception:
        original = __import__("sys").exc_info()[1]
        cleanup_errors: list[Exception] = []
        try:
            _restore_staged_artifacts(staged, rollback_dir)
        except Exception as exc:
            cleanup_errors.append(exc)
        for path in published_paths:
            try:
                if path in backups:
                    path.write_bytes(backups[path])
                elif path.exists() and path.is_file():
                    path.unlink()
            except Exception as exc:
                cleanup_errors.append(exc)
        try:
            _cleanup_tmp_dir(context["tmp_dir"], remove_tmp_dir=True)
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
        shutil.rmtree(rollback_dir, ignore_errors=True)

def cleanup_parser_artifacts(
    context: ParserExecutionContext,
    remove_tmp_dir: bool,
) -> None:
    """Apply the parser temporary-artifact cleanup policy."""
    _cleanup_tmp_dir(context["tmp_dir"], remove_tmp_dir)

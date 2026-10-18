"""Shared type contracts for parser and remapping pipeline stages."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, TypeAlias, TypedDict

import numpy as np
import pandas as pd
from numpy.typing import NDArray

if TYPE_CHECKING:
	from MDAnalysis.core.groups import AtomGroup

ParserConfigValue: TypeAlias = str | int | float
ParserParams: TypeAlias = dict[str, ParserConfigValue]
PathMap: TypeAlias = dict[str, Path]
AtomOrderEntry = tuple[int, str]
AtomOrder = list[AtomOrderEntry]
CliqueEntry = tuple[list[AtomOrderEntry], int]
SeedBundle = dict[str, int | None]
OrderVector = NDArray[np.int_]
StructureFrame = pd.DataFrame
DistanceFrame = pd.DataFrame
AngularFrame = pd.DataFrame


class ParserContext(TypedDict):
	pdb_file_name: Path
	parsed_params: ParserParams
	chosen_model: int
	chosen_chain: str
	atom_selection: str
	distance_constraints: str
	method: str
	output_mode: str
	protein_chain: AtomGroup
	angular_width: float
	backbone_torsion_percentage: float
	required_seed_fields: set[str]


class ParserExecutionContext(ParserContext):
	tmp_dir: Path
	output_files: PathMap
	tsv_structure_file: Path
	rng_interval_centered: np.random.Generator | None
	rng_talos_angles: np.random.Generator | None
	rng_phi_psi_mask: np.random.Generator | None
	rng_hh_distance_subset: np.random.Generator | None


class RemapContext(TypedDict):
	chosen_model: int
	chosen_chain: str
	out_dir_pdb_id: Path
	xfile: Path
	afile: Path
	ifile: Path


class RemapExecutionContext(RemapContext):
	df_X: StructureFrame
	df_I: DistanceFrame
	df_A: AngularFrame
	res_ids: list[int]
	nres: int
	chosen_order_id: int | str
	order_type: str
	atom_layout: str

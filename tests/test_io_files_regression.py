from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pdb_parser.io.files import (
	save_cliques_from_matrix_T,
	save_coordinates_from_df_structure,
	save_distances_from_df_structure,
)


def test_save_distances_from_df_structure_preserves_fixed_width_output(tmp_path: Path) -> None:
	df = pd.DataFrame(
		[
			[1, 2, 10, 11, 1.25, 2.5, "N", "CA", "ALA", "GLY"],
		]
	)

	output_file = tmp_path / "I_test.dat"
	save_distances_from_df_structure(df, output_file)

	assert output_file.read_text(encoding="utf-8") == (
		"    1     2     10     11   1.250000000000   2.500000000000    N   CA ALA GLY\n"
	)


def test_save_coordinates_from_df_structure_writes_only_xyz_columns(tmp_path: Path) -> None:
	df = pd.DataFrame(
		[
			[1, "N", 10, "ALA", 1.2349, 2.3451, 3.4567],
			[2, "CA", 10, "ALA", -1.0, 0.0, 4.0],
		]
	)

	output_file = tmp_path / "X_test.dat"
	save_coordinates_from_df_structure(df, output_file)

	assert output_file.read_text(encoding="utf-8") == (
		"    1.235     2.345     3.457\n"
		"   -1.000     0.000     4.000\n"
	)


def test_save_cliques_from_matrix_T_writes_expected_numeric_format(tmp_path: Path) -> None:
	T = np.array(
		[
			[1, 2, 3, 4, 1, 6.25, 7.5],
			[8, 9, 10, 11, -1, 180.0, 14.125],
		],
		dtype=float,
	)

	output_file = tmp_path / "T_test.dat"
	save_cliques_from_matrix_T(T, str(output_file))

	assert output_file.read_text(encoding="utf-8") == (
		"    1     2     3     4  1   6.250000000000   7.500000000000\n"
		"    8     9    10    11 -1 180.000000000000  14.125000000000\n"
	)


@pytest.mark.parametrize(
	("matrix", "message"),
	[
		(np.zeros((2, 6)), "7 columns"),
		(np.zeros(7), "two-dimensional"),
		(np.array([[1, 2, 3, 4, 5, np.nan, 0.0]]), "finite"),
		(np.array([[1.5, 2, 3, 4, 5, 0.0, 0.0]]), "integers"),
	],
)
def test_save_cliques_from_matrix_T_rejects_malformed_matrices(
	tmp_path: Path, matrix: np.ndarray, message: str
) -> None:
	with pytest.raises(ValueError, match=message):
		save_cliques_from_matrix_T(matrix, tmp_path / "T_invalid.dat")

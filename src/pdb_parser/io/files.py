from pathlib import Path
import pandas as pd

# -----------------------------------------------------------------------------------------------------
def read_space_separated_file(
	input_file: str | Path,
) -> pd.DataFrame:
	"""
	Read a space-separated file without header into a DataFrame.
	"""
	input_file = Path(input_file)

	if not input_file.exists():
		raise ValueError(f"Input file not found: {input_file}")

	dataframe = pd.read_csv(
		input_file,
		sep=r"\s+",
		header=None,
		engine="python",
	)

	return dataframe
# -----------------------------------------------------------------------------------------------------
def save_distances_from_df_structure(
	df_I,
	filepath,
) -> None:
	"""
	Save a distance constraint DataFrame to file using the required
	fixed-width formatting.

	Expected column convention:
	0 : atom_id_i
	1 : atom_id_j
	2 : resid_i
	3 : resid_j
	4 : d_l
	5 : d_u
	6 : atom_name_i
	7 : atom_name_j
	8 : resname_i
	9 : resname_j
	"""

	with open(filepath, "w", encoding="utf-8") as f:

		for _, row in df_I.iterrows():

			f.write(
				"%5d %5d %6d %6d %20.16f %20.16f %4.4s %4.4s %s %s\n"
				% (
					int(row[0]),
					int(row[1]),
					int(row[2]),
					int(row[3]),
					float(row[4]),
					float(row[5]),
					str(row[6]),
					str(row[7]),
					str(row[8]),
					str(row[9]),
				)
			)

	print(f"[OK] Distance constraint file saved to: {filepath}")
# -----------------------------------------------------------------------------------------------------
def save_coordinates_from_df_structure(
	df_X,
	filepath,
) -> None:
	"""
	Save the last three columns (x, y, z) of a structure DataFrame
	using the format: %.3f %.3f %.3f
	"""

	with open(filepath, "w", encoding="utf-8") as f:

		for _, row in df_X.iterrows():

			x, y, z = row.iloc[-3:]

			f.write(
				"%.3f %.3f %.3f\n"
				% (
					float(x),
					float(y),
					float(z),
				)
			)

	print(f"[OK] Coordinate file saved to: {filepath}")
# -----------------------------------------------------------------------------------------------------

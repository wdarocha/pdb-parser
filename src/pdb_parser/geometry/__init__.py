"""
Geometry utilities for the pdb_parser package.

This module provides mathematical tools related to Distance Geometry,
including transformations between torsion angles and distances and
other geometric relations used in protein backbone modeling.

The functions defined here are independent of file I/O and are intended
to support algorithms operating on molecular coordinates and distance
constraints.
"""

from .distance_geometry import (
	distances_2_abs_torsion_angle,
	torsion_angle_2_endpoint_distance,
	torsion_angle_parameters,
	torsion_angle_with_points,
)

__all__ = [
	"distances_2_abs_torsion_angle",
	"torsion_angle_2_endpoint_distance",
	"torsion_angle_parameters",
	"torsion_angle_with_points",
]


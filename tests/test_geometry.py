from __future__ import annotations

import math

import numpy as np
import pytest

from pdb_parser.geometry.distance_geometry import (
	distances_2_abs_torsion_angle,
	lambda_function,
	rho2_function,
	torsion_angle_2_endpoint_distance,
	torsion_angle_parameters,
	torsion_angle_with_points,
)


def test_lambda_function_rejects_zero_denominator() -> None:
	with pytest.raises(ValueError, match="division by zero"):
		lambda_function(0.0, 1.0, 1.0)


def test_rho2_function_clamps_tiny_negative_roundoff() -> None:
	assert rho2_function(1.0, 1.0 + 1e-13) == 0.0


def test_torsion_angle_round_trip_from_points() -> None:
	x1 = np.array([0.0, 0.0, 0.0])
	x2 = np.array([1.0, 0.0, 0.0])
	x3 = np.array([1.0, 1.0, 0.0])
	x4 = np.array([1.0, 1.0, 1.0])

	tau = torsion_angle_with_points(x1, x2, x3, x4)
	d14 = torsion_angle_2_endpoint_distance(x1, x2, x3, x4, tau)

	expected_d14 = float(np.linalg.norm(x1 - x4))
	assert math.isclose(d14, expected_d14, rel_tol=1e-12, abs_tol=1e-12)

	d12 = float(np.linalg.norm(x1 - x2))
	d13 = float(np.linalg.norm(x1 - x3))
	d23 = float(np.linalg.norm(x2 - x3))
	d24 = float(np.linalg.norm(x2 - x4))
	d34 = float(np.linalg.norm(x3 - x4))

	abs_tau = distances_2_abs_torsion_angle(d12, d13, d14, d23, d24, d34)
	assert math.isclose(abs(abs(tau)), abs_tau, rel_tol=1e-12, abs_tol=1e-12)


def test_torsion_angle_parameters_return_finite_values_for_valid_geometry() -> None:
	p, q = torsion_angle_parameters(
		d12=1.0,
		d13=math.sqrt(2.0),
		d23=1.0,
		d24=math.sqrt(2.0),
		d34=1.0,
	)

	assert math.isfinite(p)
	assert math.isfinite(q)
	assert q >= 0.0

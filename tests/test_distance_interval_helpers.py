from __future__ import annotations

import math

import numpy as np
import pytest

from pdb_parser.io.distance_constraints import (
	_get_centered_interval,
	_get_experimental_interval,
	_get_precise_interval,
	_get_xray_ca_centered_interval,
)


def test_get_precise_interval_respects_cutoff() -> None:
	assert _get_precise_interval(3.0, 5.0) == (3.0, 3.0)
	assert _get_precise_interval(5.0, 5.0) is None


def test_get_experimental_interval_classifies_by_noe_threshold() -> None:
	assert _get_experimental_interval(2.0, 2.5, 3.5, 5.0, 1.8) == (1.8, 2.5)
	assert _get_experimental_interval(3.0, 2.5, 3.5, 5.0, 1.8) == (1.8, 3.5)
	assert _get_experimental_interval(4.0, 2.5, 3.5, 5.0, 1.8) == (1.8, 5.0)
	assert _get_experimental_interval(5.0, 2.5, 3.5, 5.0, 1.8) is None


@pytest.mark.parametrize(
	("distance", "expected_lower"),
	[
		(2.0, 1.8),  # vdwr_hh already below the reference distance
		(1.8, 0.9 * 1.8),
		(1.4376, 0.7 * 1.8),  # 1BBA: 80% (1.44) is not below 1.4376
		(0.9, 0.4 * 1.8),  # 50% (0.9) is not strictly below 0.9
		(0.1, 0.9 * 0.1),  # not even 10% (0.18) is below: 90% of the distance
	],
)
def test_experimental_lower_bound_never_exceeds_the_reference_distance(distance: float, expected_lower: float) -> None:
	"""NMR models have H-H pairs closer than vdwr_hh (1.44 A in 1BBA); the reference
	structure must still satisfy its interval."""
	lower, upper = _get_experimental_interval(distance, 2.5, 3.5, 5.0, 1.8)

	assert lower == pytest.approx(expected_lower)
	assert lower < distance < upper


def test_lower_bound_below_steps_by_tenths_of_the_base_then_falls_back_to_the_distance() -> None:
	from pdb_parser.io.vdw import lower_bound_below

	assert lower_bound_below(3.0, 2.4) == 2.4
	assert lower_bound_below(2.2, 2.4) == pytest.approx(0.9 * 2.4)
	assert lower_bound_below(1.0, 2.4) == pytest.approx(0.4 * 2.4)
	assert lower_bound_below(0.2, 2.4) == pytest.approx(0.9 * 0.2)
	for distance in (0.01, 0.5, 1.0, 2.399, 2.4, 5.0):
		assert lower_bound_below(distance, 2.4) < distance


def test_get_centered_interval_returns_interval_containing_distance() -> None:
	rng = np.random.default_rng(12345)

	interval = _get_centered_interval(
		dist=4.0,
		resid_i=10,
		resid_j=12,
		epsilon_short=0.5,
		epsilon_long=2.0,
		vdw_threshold=1.8,
		max_distance=6.0,
		rng=rng,
	)

	assert interval is not None
	d_l, d_u = interval
	assert d_l < 4.0 < d_u
	assert d_l >= 1.8
	assert d_u <= 6.0
	assert math.isclose(d_u - d_l, 2.0, rel_tol=0.0, abs_tol=1e-9)


def test_get_centered_interval_returns_none_when_cutoff_is_not_above_distance() -> None:
	assert _get_centered_interval(
		dist=5.0,
		resid_i=1,
		resid_j=2,
		epsilon_short=0.5,
		epsilon_long=1.0,
		vdw_threshold=1.8,
		max_distance=5.0,
		rng=np.random.default_rng(1),
	) is None


def test_centered_interval_helpers_reject_degenerate_inputs_without_looping() -> None:
	assert _get_centered_interval(
		dist=4.0,
		resid_i=1,
		resid_j=3,
		epsilon_short=1.0,
		epsilon_long=0.0,
		vdw_threshold=1.8,
		max_distance=6.0,
	) is None
	assert _get_xray_ca_centered_interval(4.0, float("nan"), 6.0) is None


def test_vectorized_centered_sampling_consumes_the_random_stream_like_the_scalar_loop() -> None:
	from pdb_parser.io.distance_intervals import get_centered_interval, sample_centered_intervals
	from pdb_parser.io.vdw import lower_bounds_below

	rng = np.random.default_rng(2)
	count = 200_000
	dist = rng.uniform(1.0, 4.9, count)
	adjacent = rng.random(count) < 0.3
	epsilon = np.where(adjacent, 1.0, 2.0)

	scalar_rng = np.random.default_rng(99)
	expected, rejected_pairs = [], 0
	for value, eps in zip(dist.tolist(), epsilon.tolist()):
		before = scalar_rng.bit_generator.state
		expected.append(get_centered_interval(value, 1, 1 if eps == 1.0 else 5, 1.0, 2.0, 1.8, 5.0, scalar_rng))
		probe = np.random.default_rng()
		probe.bit_generator.state = before
		probe.normal()
		rejected_pairs += probe.bit_generator.state != scalar_rng.bit_generator.state

	vector_rng = np.random.default_rng(99)
	d_l, d_u, accepted = sample_centered_intervals(dist, epsilon, lower_bounds_below(dist, np.full(count, 1.8)), 5.0, vector_rng)

	# Rejections (draws beyond 4 sigma) happen, so the rewind path is exercised.
	assert rejected_pairs > 0
	assert accepted.all() and all(interval is not None for interval in expected)
	assert np.array_equal(d_l, np.array([interval[0] for interval in expected]))
	assert np.array_equal(d_u, np.array([interval[1] for interval in expected]))
	assert vector_rng.bit_generator.state == scalar_rng.bit_generator.state

"""Numerically safe interval builders for long-range distance constraints."""

from __future__ import annotations

import math

import numpy as np

from pdb_parser.io.sequence import are_adjacent_residues
from pdb_parser.io.vdw import VDW_CARBON_RADIUS, lower_bound_below

MAX_INTERVAL_SAMPLING_ATTEMPTS = 10_000


def get_centered_interval(
	dist: float,
	resid_i: int,
	resid_j: int,
	epsilon_short: float,
	epsilon_long: float,
	vdw_threshold: float,
	max_distance: float,
	rng: np.random.Generator | None = None,
	residue_ids: tuple[int, ...] | None = None,
) -> tuple[float, float] | None:
	"""Build a valid centred interval for an NMR-derived reference distance."""
	if residue_ids is None:
		adjacent = abs(resid_i - resid_j) <= 1
	else:
		adjacent = resid_i == resid_j or are_adjacent_residues(
			residue_ids, resid_i, resid_j, context="validated residue sequence"
		) or are_adjacent_residues(
			residue_ids, resid_j, resid_i, context="validated residue sequence"
		)
	epsilon = epsilon_short if adjacent else epsilon_long
	if not all(math.isfinite(value) for value in (dist, epsilon, vdw_threshold, max_distance)):
		return None
	if epsilon <= 0.0 or vdw_threshold <= 0.0 or max_distance <= 0.0 or max_distance <= dist:
		return None

	adjusted_vdw_threshold = lower_bound_below(dist, vdw_threshold)

	rng = rng if rng is not None else np.random.default_rng()
	for _ in range(MAX_INTERVAL_SAMPLING_ATTEMPTS):
		sampled_distance = float(rng.normal(loc=dist, scale=epsilon / 8.0))
		d_l = max(sampled_distance - epsilon / 2.0, adjusted_vdw_threshold)
		d_u = min(sampled_distance + epsilon / 2.0, max_distance)
		if d_l < dist < d_u:
			return d_l, d_u
	return None


def sample_centered_intervals(
	dist: np.ndarray,
	epsilon: np.ndarray,
	lower_vdw: np.ndarray,
	max_distance: float,
	rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
	"""Vectorized sampling of ``get_centered_interval`` for pairs that reach the sampling loop.

	Returns ``(d_l, d_u, accepted)``. The random stream is consumed exactly as by one
	``get_centered_interval`` call per pair, in order: the pairs are drawn in a block, and at
	the first rejected draw the generator is rewound, advanced to that pair, and the pair is
	redrawn one sample at a time, as the scalar loop does. Rejections are rare (a draw more
	than 4 standard deviations away), so almost every pair comes from a block draw.
	"""
	count = len(dist)
	d_l = np.zeros(count)
	d_u = np.zeros(count)
	accepted = np.zeros(count, dtype=bool)
	scale = epsilon / 8.0
	half = epsilon / 2.0
	start = 0
	while start < count:
		state = rng.bit_generator.state
		sampled = rng.normal(loc=dist[start:], scale=scale[start:])
		block_l = np.maximum(sampled - half[start:], lower_vdw[start:])
		block_u = np.minimum(sampled + half[start:], max_distance)
		ok = (block_l < dist[start:]) & (dist[start:] < block_u)
		first_rejected = count - start if ok.all() else int(np.argmin(ok))
		stop = start + first_rejected
		d_l[start:stop], d_u[start:stop] = block_l[:first_rejected], block_u[:first_rejected]
		accepted[start:stop] = True
		if stop == count:
			break
		# Rewind and consume the draws up to the rejected pair's first attempt.
		rng.bit_generator.state = state
		rng.normal(loc=dist[start:stop + 1], scale=scale[start:stop + 1])
		for _ in range(MAX_INTERVAL_SAMPLING_ATTEMPTS - 1):
			value = float(rng.normal(loc=dist[stop], scale=scale[stop]))
			lower = max(value - half[stop], lower_vdw[stop])
			upper = min(value + half[stop], max_distance)
			if lower < dist[stop] < upper:
				d_l[stop], d_u[stop], accepted[stop] = lower, upper, True
				break
		start = stop + 1
	return d_l, d_u, accepted


def get_experimental_interval(
	dist: float,
	noe_strong: float,
	noe_medium: float,
	noe_weak: float,
	vdw_threshold: float,
) -> tuple[float, float] | None:
	"""Return the NMR NOESY interval corresponding to a reference distance.

	The upper bound is the NOE class of ``dist``; the lower bound is ``vdw_threshold``,
	brought below the reference distance by ``lower_bound_below`` when needed.
	"""
	upper = next((bound for bound in (noe_strong, noe_medium, noe_weak) if dist < bound), None)
	if upper is None:
		return None
	return lower_bound_below(dist, vdw_threshold), upper


def get_precise_interval(dist: float, max_distance: float) -> tuple[float, float] | None:
	"""Return the exact interval when the reference distance is within the cutoff."""
	if dist >= max_distance:
		return None
	return dist, dist


def get_xray_ca_centered_interval(
	dist: float,
	epsilon_long: float,
	max_distance: float,
	rng: np.random.Generator | None = None,
) -> tuple[float, float] | None:
	"""Build a centred CA--CA interval with the crystallographic VDW policy."""
	if not all(math.isfinite(value) for value in (dist, epsilon_long, max_distance)):
		return None
	if epsilon_long <= 0.0 or max_distance <= 0.0 or max_distance <= dist:
		return None

	lower_vdw = lower_bound_below(dist, 2.0 * VDW_CARBON_RADIUS)

	rng = rng if rng is not None else np.random.default_rng()
	for _ in range(MAX_INTERVAL_SAMPLING_ATTEMPTS):
		sampled_distance = float(rng.normal(loc=dist, scale=epsilon_long / 8.0))
		d_l = max(sampled_distance - epsilon_long / 2.0, lower_vdw)
		d_u = min(sampled_distance + epsilon_long / 2.0, max_distance)
		if d_l < dist < d_u:
			return d_l, d_u
	return None

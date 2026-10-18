"""Shared van der Waals constants used by constraint builders."""

import numpy as np

VDW_RADII = {
	"H": 1.20,
	"C": 1.70,
	"N": 1.55,
	"O": 1.52,
}
VDW_CARBON_RADIUS = VDW_RADII["C"]


def lower_bound_below(distance: float, base: float) -> float:
	"""Return a lower bound below the reference ``distance``, starting from ``base``.

	The reference structure must satisfy every interval built from it. When ``base`` is
	not below ``distance``, try 90%, 80%, ..., 10% of ``base`` and keep the first value
	below ``distance``; if even 10% is not, use 90% of ``distance``.
	"""
	if base < distance:
		return base
	for percentage in range(90, 0, -10):
		candidate = (percentage / 100.0) * base
		if candidate < distance:
			return candidate
	return 0.9 * distance


def lower_bounds_below(distance: np.ndarray, base: np.ndarray) -> np.ndarray:
	"""Vectorized ``lower_bound_below``: the same rule and the same floating-point operations."""
	result = np.array(base, dtype=float, copy=True)
	pending = ~(result < distance)
	base = np.asarray(base, dtype=float)
	for percentage in range(90, 0, -10):
		candidate = (percentage / 100.0) * base
		accepted = pending & (candidate < distance)
		result[accepted] = candidate[accepted]
		pending &= ~accepted
	result[pending] = 0.9 * distance[pending]
	return result

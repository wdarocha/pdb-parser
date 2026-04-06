import numpy as np

# -----------------------------------------------------------------------------------------------------
def lambda_function(dwv: float, dzv: float, dzw: float) -> float:
	"""
	Compute the lambda parameter associated with a triangle projection.

	This quantity appears in the torsion-to-distance relation for four points
	and corresponds to the projected coordinate along the segment of length dwv.
	"""
	if abs(dwv) <= 1e-15:
		raise ValueError(f"Invalid value dwv={dwv:.16e}: division by zero in lambda_function.")

	return (dzv * dzv + dwv * dwv - dzw * dzw) / (2.0 * dwv)

# -----------------------------------------------------------------------------------------------------
def rho2_function(dzv: float, lambda_value: float) -> float:
	"""
	Compute the squared rho parameter.

	In exact arithmetic, this quantity should be nonnegative. Because of floating-point
	errors, small negative values may appear for nearly degenerate configurations.
	"""
	rho2_value = dzv * dzv - lambda_value * lambda_value

	# Clamp tiny negative values caused by numerical roundoff.
	if rho2_value < 0.0 and abs(rho2_value) <= 1e-12:
		return 0.0

	return rho2_value

# -----------------------------------------------------------------------------------------------------
def torsion_angle_parameters(
	d12: float,
	d13: float,
	d23: float,
	d24: float,
	d34: float,
) -> tuple[float, float]:
	"""
	Compute the parameters p and q in the torsion-distance relation.

	For four points x1, x2, x3, x4, the endpoint distance d14 satisfies

		d14^2 = p - 2 q cos(tau),

	where tau is the torsion angle defined by the quadruple.

	Equivalently,

		cos(tau) = (p - d14^2) / (2 q).
	"""
	lambda_bar = lambda_function(d23, d12, d13)
	rho_bar2 = rho2_function(d12, lambda_bar)

	lambda_value = lambda_function(d23, d24, d34)
	rho2_value = rho2_function(d24, lambda_value)

	if rho_bar2 < 0.0:
		raise ValueError(f"Invalid rho_bar2={rho_bar2:.16e}: geometry is not numerically feasible.")

	if rho2_value < 0.0:
		raise ValueError(f"Invalid rho2={rho2_value:.16e}: geometry is not numerically feasible.")

	delta_lambda = lambda_value - lambda_bar
	p = delta_lambda * delta_lambda + rho2_value + rho_bar2
	q = np.sqrt(rho2_value * rho_bar2)

	return float(p), float(q)

# -----------------------------------------------------------------------------------------------------
def torsion_angle_2_endpoint_distance(
	x1: np.ndarray,
	x2: np.ndarray,
	x3: np.ndarray,
	x4: np.ndarray,
	tau_rad: float,
) -> float:
	"""
	Convert a torsion angle tau into the corresponding endpoint distance d14.

	The input points define the local geometry through the five distances
	d12, d13, d23, d24, and d34, and tau determines the sixth distance d14.
	"""
	d12 = float(np.linalg.norm(x1 - x2))
	d13 = float(np.linalg.norm(x1 - x3))
	d23 = float(np.linalg.norm(x2 - x3))
	d24 = float(np.linalg.norm(x2 - x4))
	d34 = float(np.linalg.norm(x3 - x4))

	p, q = torsion_angle_parameters(d12, d13, d23, d24, d34)

	# In exact arithmetic:
	#	d14^2 = p - 2 q cos(tau)
	d14_sq = p - 2.0 * q * np.cos(tau_rad)

	# Clamp tiny negative values caused by numerical roundoff.
	if d14_sq < 0.0 and abs(d14_sq) <= 1e-12:
		d14_sq = 0.0
	elif d14_sq < 0.0:
		raise ValueError(f"Invalid squared distance d14_sq={d14_sq:.16e}.")

	return float(np.sqrt(d14_sq))

# -----------------------------------------------------------------------------------------------------
def distances_2_abs_torsion_angle(
	d12: float,
	d13: float,
	d14: float,
	d23: float,
	d24: float,
	d34: float,
) -> float:
	"""
	Convert the endpoint distance d14 into the corresponding absolute torsion angle.

	The returned value belongs to [0, pi], since acos gives only the absolute torsion.
	If the signed torsion is needed, an additional orientation rule must be used.
	"""
	p, q = torsion_angle_parameters(d12, d13, d23, d24, d34)

	if abs(q) <= 1e-15:
		raise ValueError("Degenerate configuration: q is zero, so the torsion angle is undefined.")

	cos_tau = (p - d14 * d14) / (2.0 * q)

	# Numerical safety: due to roundoff, cos_tau may be slightly outside [-1, 1].
	cos_tau = float(np.clip(cos_tau, -1.0, 1.0))

	return float(np.arccos(cos_tau))
# -----------------------------------------------------------------------------------------------------
def abs_torsion_angle_with_points(x1, x2, x3, x4):
	d12 = np.linalg.norm(x1 - x2)
	d13 = np.linalg.norm(x1 - x3)
	d14 = np.linalg.norm(x1 - x4)
	d23 = np.linalg.norm(x2 - x3)
	d24 = np.linalg.norm(x2 - x4)
	d34 = np.linalg.norm(x3 - x4)
	
	return distances_2_abs_torsion_angle(d12, d13, d14, d23, d24, d34)
# -----------------------------------------------------------------------------------------------------	
def sign_torsion_angle(x1, x2, x3, x4):
	normal_plane = np.cross(x3 - x2, x1 - x2)
	direction = x4 - x2
	
	return np.sign(np.dot(normal_plane, direction))
# -----------------------------------------------------------------------------------------------------
def torsion_angle_with_points(x1, x2, x3, x4):
	
	return sign_torsion_angle(x1, x2, x3, x4) * abs_torsion_angle_with_points(x1, x2, x3, x4)


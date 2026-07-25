"""Target-ionization rate analysis for MarsASPEN crossing estimators.

For target species j, MarsASPEN evaluates

    q_j = n_j sum_i Wn_i |Vr_i| sigma_ion(E_i)

at each spherical altitude surface. Contributions are retained separately for
H ENA and H+ projectiles and for downward and upward crossings.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np


def ionization_rate_from_components(
    target_density_m3: np.ndarray,
    flux_times_cross_section_s1: np.ndarray,
) -> dict[str, np.ndarray]:
    """Return charge-resolved and total ionization rates.

    ``flux_times_cross_section_s1`` has dimensions
    ``(altitude, charge, direction)``. Charge order is H ENA then H+, and
    direction order is downward then upward. Both directions are added because
    either population can ionize the target.
    """
    density = np.asarray(target_density_m3, dtype=float).squeeze()
    coefficient = np.asarray(flux_times_cross_section_s1, dtype=float)
    expected = (density.size, 2, 2)
    if coefficient.shape != expected and coefficient.shape == expected[::-1]:
        coefficient = coefficient.transpose(2, 1, 0)
    if coefficient.shape != expected:
        raise ValueError(
            f"unexpected coefficient shape {coefficient.shape}; "
            f"expected {expected}"
        )
    if np.any(density < 0.0) or np.any(coefficient < 0.0):
        raise ValueError("density and directional coefficients must be nonnegative")

    directional_rate = coefficient * density[:, None, None]
    by_charge = directional_rate.sum(axis=2)
    return {
        "directional_rate_m3_s1": directional_rate,
        "rate_by_charge_m3_s1": by_charge,
        "total_rate_m3_s1": by_charge.sum(axis=1),
    }


def ionization_rate_from_mat(
    data: Mapping[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Extract and verify an ionization-rate result from a Julia MAT mapping."""
    altitude = np.asarray(data["altitude_surfaces_km"], dtype=float).squeeze()
    density = np.asarray(data["target_density_m3"], dtype=float).squeeze()
    coefficient = np.asarray(
        data["flux_times_ionization_cross_section_s1"], dtype=float
    )
    calculated = ionization_rate_from_components(density, coefficient)
    return {
        "altitude_km": altitude,
        "target_density_m3": density,
        "flux_times_cross_section_s1": coefficient,
        **calculated,
    }

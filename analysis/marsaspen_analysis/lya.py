"""H Ly-alpha volume emission, radiative energy, and limb brightness."""

from __future__ import annotations

from typing import Mapping

import numpy as np

MARS_RADIUS_KM = 3388.25
RAYLEIGH_COLUMN_M2 = 1.0e10


def shell_edges_from_centers(altitude_km: np.ndarray) -> np.ndarray:
    """Return shell edges for strictly increasing altitude centers."""
    centers = np.asarray(altitude_km, dtype=float).squeeze()
    if centers.ndim != 1 or centers.size < 2 or np.any(np.diff(centers) <= 0):
        raise ValueError("altitude_km must contain increasing centers")
    edges = np.empty(centers.size + 1)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def optically_thin_limb_brightness_rayleigh(
    altitude_km: np.ndarray,
    volume_emission_rate_photons_m3_s1: np.ndarray,
    tangent_altitude_km: np.ndarray | None = None,
    mars_radius_km: float = MARS_RADIUS_KM,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate a spherical VER profile along limb lines of sight.

    The volume emission rate is assumed constant inside each spherical shell
    and isotropic. No absorption or H Ly-alpha resonant scattering is applied.
    One Rayleigh corresponds to a column emission of 1e10 photons m^-2 s^-1.
    """
    altitude = np.asarray(altitude_km, dtype=float).squeeze()
    emission = np.asarray(volume_emission_rate_photons_m3_s1, dtype=float)
    if emission.shape[0] != altitude.size:
        raise ValueError("the first emission dimension must be altitude")
    tangent = np.atleast_1d(
        altitude.copy()
        if tangent_altitude_km is None
        else np.asarray(tangent_altitude_km, dtype=float).squeeze()
    )
    edges_m = (mars_radius_km + shell_edges_from_centers(altitude)) * 1000.0
    tangent_radius_m = (mars_radius_km + tangent) * 1000.0
    flat_emission = emission.reshape(altitude.size, -1)
    brightness = np.zeros((tangent.size, flat_emission.shape[1]))

    for it, impact_parameter in enumerate(tangent_radius_m):
        for shell in range(altitude.size):
            inner_radius = edges_m[shell]
            outer_radius = edges_m[shell + 1]
            if outer_radius <= impact_parameter:
                continue
            visible_inner = max(inner_radius, impact_parameter)
            outer_term = np.sqrt(max(outer_radius**2 - impact_parameter**2, 0.0))
            inner_term = np.sqrt(max(visible_inner**2 - impact_parameter**2, 0.0))
            path_length_m = 2.0 * (outer_term - inner_term)
            brightness[it] += flat_emission[shell] * path_length_m

    brightness /= RAYLEIGH_COLUMN_M2
    output_shape = (tangent.size,) + emission.shape[1:]
    return tangent, brightness.reshape(output_shape)


def lya_profiles_from_mat(data: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Extract charge-resolved Ly-alpha profiles from a Julia MAT result."""
    altitude = np.asarray(data["altitude_surfaces_km"], dtype=float).squeeze()
    by_charge = np.asarray(
        data["volume_emission_rate_by_charge_photons_m3_s1"], dtype=float
    )
    expected = (altitude.size, 2)
    if by_charge.shape != expected and by_charge.shape == expected[::-1]:
        by_charge = by_charge.T
    if by_charge.shape != expected:
        raise ValueError(
            f"unexpected charge-resolved VER shape {by_charge.shape}; "
            f"expected {expected}"
        )
    total = by_charge.sum(axis=1)
    photon_energy_j = float(np.asarray(data["photon_energy_j"]).squeeze())
    tangent, limb_by_charge = optically_thin_limb_brightness_rayleigh(
        altitude, by_charge, altitude
    )
    return {
        "altitude_km": altitude,
        "volume_emission_rate_by_charge_photons_m3_s1": by_charge,
        "total_volume_emission_rate_photons_m3_s1": total,
        "photon_energy_j": photon_energy_j,
        "radiative_energy_rate_by_charge_w_m3": by_charge * photon_energy_j,
        "total_radiative_energy_rate_w_m3": total * photon_energy_j,
        "tangent_altitude_km": tangent,
        "limb_brightness_by_charge_rayleigh": limb_by_charge,
        "total_limb_brightness_rayleigh": limb_by_charge.sum(axis=1),
    }

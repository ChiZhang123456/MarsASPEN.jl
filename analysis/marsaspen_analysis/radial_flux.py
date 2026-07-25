"""Radial number-flux analysis for MarsASPEN altitude-surface output.

MarsASPEN currently diagnoses flux through spherical altitude surfaces. For a
particle at position ``r`` with velocity ``v``, the surface-normal direction
is radial:

    Vr = dot(r, v) / norm(r)

and its signed outward number-flux contribution is ``Wn * Vr``. ``Wn`` has
units m^-3 and ``Vr`` has units m s^-1, giving flux in m^-2 s^-1.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np


def local_radial_velocity(
    position_m: np.ndarray,
    velocity_m_s: np.ndarray,
) -> np.ndarray:
    """Return the local radial velocity in m s^-1.

    Positive values point away from Mars and negative values point toward
    Mars. The final dimension of both inputs must contain the three Cartesian
    components.
    """
    position = np.asarray(position_m, dtype=float)
    velocity = np.asarray(velocity_m_s, dtype=float)
    if position.shape != velocity.shape or position.shape[-1] != 3:
        raise ValueError(
            "position_m and velocity_m_s must have matching (..., 3) shapes"
        )
    radius = np.linalg.norm(position, axis=-1)
    if np.any(radius == 0.0):
        raise ValueError("local radial velocity is undefined at r = 0")
    return np.sum(position * velocity, axis=-1) / radius


def particle_radial_flux(
    density_weight_m3: np.ndarray | float,
    radial_velocity_m_s: np.ndarray | float,
) -> np.ndarray:
    """Return signed particle radial-flux contributions in m^-2 s^-1."""
    weight = np.asarray(density_weight_m3, dtype=float)
    velocity = np.asarray(radial_velocity_m_s, dtype=float)
    return weight * velocity


def radial_flux_profiles(
    radial_flux_m2_s: np.ndarray,
    n_altitudes: int | None = None,
) -> dict[str, np.ndarray]:
    """Calculate directional and signed profiles from Julia crossing output.

    Parameters
    ----------
    radial_flux_m2_s
        Array with dimensions ``(altitude, charge, direction)``. Charge order
        is H ENA then H+, and direction order is downward then upward. Arrays
        loaded from MAT v7.3 with reversed dimensions are corrected when
        ``n_altitudes`` is provided.
    n_altitudes
        Expected number of altitude surfaces.

    Returns
    -------
    dict
        ``downward`` and ``upward`` are positive magnitudes.
        ``signed_outward`` equals upward minus downward.
        ``net_downward`` equals downward minus upward.
    """
    flux = np.asarray(radial_flux_m2_s, dtype=float)
    if n_altitudes is not None:
        expected = (int(n_altitudes), 2, 2)
        if flux.shape != expected and flux.shape == expected[::-1]:
            flux = flux.transpose(2, 1, 0)
        if flux.shape != expected:
            raise ValueError(
                f"unexpected radial flux shape {flux.shape}; expected {expected}"
            )
    elif flux.ndim != 3 or flux.shape[1:] != (2, 2):
        raise ValueError(
            "radial_flux_m2_s must have shape (altitude, 2, 2)"
        )
    if np.any(flux < 0.0):
        raise ValueError("directional radial-flux magnitudes cannot be negative")

    downward = flux[:, :, 0]
    upward = flux[:, :, 1]
    return {
        "radial_flux_m2_s": flux,
        "downward": downward,
        "upward": upward,
        "signed_outward": upward - downward,
        "net_downward": downward - upward,
    }


def radial_flux_from_mat(
    data: Mapping[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Extract altitude and radial-flux profiles from a loaded MAT mapping."""
    altitude = np.asarray(data["altitude_surfaces_km"], dtype=float).squeeze()
    if altitude.ndim != 1:
        raise ValueError("altitude_surfaces_km must be one-dimensional")
    profiles = radial_flux_profiles(
        data["radial_flux_m2_s"], n_altitudes=altitude.size
    )
    return {"altitude_km": altitude, **profiles}

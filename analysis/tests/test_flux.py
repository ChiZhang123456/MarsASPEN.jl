"""Tests for local vertical number-flux analysis."""

import numpy as np

from marsaspen_analysis.flux import (
    local_vertical_velocity,
    particle_vertical_flux,
    vertical_flux_from_mat,
)


def test_local_vertical_velocity_uses_radial_direction():
    position = np.array([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0]])
    velocity = np.array([[-4.0, 8.0, 0.0], [5.0, 6.0, 0.0]])
    np.testing.assert_allclose(
        local_vertical_velocity(position, velocity), [-4.0, 6.0]
    )


def test_particle_vertical_flux_units_and_sign():
    flux = particle_vertical_flux(
        np.array([2.0e6, 3.0e6]), np.array([-4.0e5, 1.0e5])
    )
    np.testing.assert_allclose(flux, [-8.0e11, 3.0e11])


def test_vertical_flux_from_mat_direction_and_charge_order():
    raw = np.zeros((3, 2, 2))
    raw[:, 0, 0] = [4.0, 5.0, 6.0]
    raw[:, 0, 1] = [1.0, 2.0, 3.0]
    raw[:, 1, 0] = [8.0, 9.0, 10.0]
    raw[:, 1, 1] = [2.0, 3.0, 4.0]
    result = vertical_flux_from_mat({
        "altitude_surfaces_km": np.array([80.5, 81.5, 82.5]),
        "radial_flux_m2_s": raw,
    })
    np.testing.assert_allclose(result["signed_outward"][:, 0], [-3, -3, -3])
    np.testing.assert_allclose(result["net_downward"][:, 1], [6, 6, 6])

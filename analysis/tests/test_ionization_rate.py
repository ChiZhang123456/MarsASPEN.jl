"""Tests for the target-ionization crossing estimator."""

import numpy as np

from marsaspen_analysis import ionization_rate_from_components


def test_ionization_rate_units_and_charge_sum():
    density = np.array([2.0e10, 3.0e10])
    coefficient = np.zeros((2, 2, 2))
    coefficient[:, 0, :] = [[1.0e-3, 2.0e-3], [2.0e-3, 3.0e-3]]
    coefficient[:, 1, :] = [[4.0e-3, 5.0e-3], [6.0e-3, 7.0e-3]]
    result = ionization_rate_from_components(density, coefficient)

    np.testing.assert_allclose(
        result["rate_by_charge_m3_s1"][0], [6.0e7, 1.8e8]
    )
    np.testing.assert_allclose(
        result["total_rate_m3_s1"], [2.4e8, 5.4e8]
    )

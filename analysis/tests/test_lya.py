"""Tests for Ly-alpha profile and Rayleigh conversion."""

import numpy as np

from marsaspen_analysis import optically_thin_limb_brightness_rayleigh


def test_single_visible_shell_limb_path():
    altitude = np.array([100.5, 101.5])
    emission = np.array([2.0e6, 0.0])
    tangent, brightness = optically_thin_limb_brightness_rayleigh(
        altitude, emission, np.array([100.5])
    )
    radius = (3388.25 + 100.5) * 1000.0
    outer = (3388.25 + 101.0) * 1000.0
    expected_path = 2.0 * np.sqrt(outer**2 - radius**2)
    np.testing.assert_allclose(tangent, [100.5])
    np.testing.assert_allclose(
        brightness, [2.0e6 * expected_path / 1.0e10]
    )

# mpc_boilerplate/tests/test_refraction.py
# import pytest

# Third-party imports
import numpy as np

# Import the specific package/module/function we are testing
from mpc_boilerplate.refraction import snell


def test_perpendicular():
    """
    For any indexes, a ray normal to the surface should not bend.
    We'll try a couple different combinations of indexes....
    """
    actual = snell(0, 2.00, 3.00)
    expected = 0
    assert actual == expected

    actual = snell(0, 3.00, 2.00)
    expected = 0
    assert actual == expected


def test_air_water():
    """
    Test a specific example
    """
    n_air, n_water = 1.00, 1.33
    actual = snell(np.deg2rad(45), n_air, n_water)
    expected = 0.5605584137424605
    assert np.allclose(actual, expected)

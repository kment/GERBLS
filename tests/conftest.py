import gerbls
import numpy as np
import pytest


def generate_LD_light_curve():
    """
    Generates a limb-darkened light curve for testing and returns it.
    """
    # Orbital parameters
    b = 0.
    mag0 = 1.
    P = 1.37
    r = 0.005**0.5
    t0 = 0.38

    # Stellar and limb darkening parameters
    target = gerbls.pyTarget()
    target.u1 = 0.4
    target.u2 = 0.3

    model = gerbls.LDModel(b=b, mag0=mag0, P=P, r=r, t0=t0, target=target)
    time = np.arange(0., 30., 2. / 60 / 24)
    err = np.random.randn(len(time)) * 0.0003 + 0.003
    mag = np.random.randn(len(time)) * err + model.eval(time)

    phot = gerbls.pyDataContainer()
    phot.store(time, mag, err)

    return phot


@pytest.fixture
def phot_test():
    """
    Returns a generated limb-darkened light curve for testing.
    """
    phot = generate_LD_light_curve()

    assert phot.size == 21600
    return phot


@pytest.fixture
def phot_test_from_file():
    """
    Loads in a light curve from a file for testing and returns it.
    """
    data = np.loadtxt("tests/phottest.dat")

    phot = gerbls.pyDataContainer()
    phot.store(data[:, 0], data[:, 1], data[:, 2])

    assert phot.size == 21600
    return phot

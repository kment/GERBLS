import gerbls
import numpy as np
import pytest

@pytest.fixture
def phot_test():
    """
    Loads in a light curve for testing and returns it.
    """

    data = np.loadtxt("tests/phottest.dat")

    phot = gerbls.pyDataContainer()
    phot.store(data[:,0], data[:,1], data[:,2])

    assert phot.size == 21600
    return phot
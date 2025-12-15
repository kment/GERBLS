from __future__ import annotations
import gerbls
import numpy as np
import numpy.typing as npt
from typing import Union


def run_bls(time: npt.ArrayLike,
            mag: npt.ArrayLike,
            err: npt.ArrayLike,
            min_period: float,
            max_period: float,
            durations: Union[list, float],
            t_samp: float = 0.):
    """
    A basic convenience function to generate a BLS spectrum.
    The input data will be resampled to a uniform time cadence to run the BLS,
    use ``t_samp`` to explicitly provide the cadence for resampling.
    Automatic downsampling will be used to maintain an optimal period spacing.

    Parameters
    ----------
    time : ArrayLike
        Array of observation timestamps.
    mag : ArrayLike
        Array of observed fluxes.
    err : ArrayLike
        Array of flux uncertainties for each observation.
    min_period : float
        Minimum BLS period to search.
    max_period : float
        Maximum BLS period to search.
    durations : list or float
        List of transit durations to test at each period.
    t_samp : float, optional
        Time sampling to bin the data before running the BLS.
        If 0 (default), the median time difference between observations is used.

    Returns
    -------
    dict
        A dictionary with BLS results:

        ========= ===================================================
        Key       Value
        ========= ===================================================
        ``P``     list of tested periods
        ``dchi2`` BLS statistic (:math:`\Delta\chi^2`) at each period
        ``t0``    best-fit transit mid-point at each period
        ``dur``   best-fit duration at each period
        ``mag0``  best-fit flux baseline at each period
        ``dmag``  best-fit transit depth at each period
        ``snr``   estimated SNR at each period
        ========= ===================================================
    """
    # Input checks
    assert len(time) == len(mag) == len(err), "time, mag, and err must have the same length."

    # Make sure the data is time-sorted and formatted as Numpy arrays
    if np.all(np.diff(time) >= 0):
        time = np.array(time)
        mag = np.array(mag)
        err = np.array(err)
    else:
        order = np.argsort(time)
        time = np.array(time)[order]
        mag = np.array(mag)[order]
        err = np.array(err)[order]
    
    if not hasattr(durations, "__len__"):
        durations = [durations]

    # Create a GERBLS data container
    phot = gerbls.pyDataContainer()
    phot.store(time, mag, err, convert_to_flux=False)

    # Set up and run the BLS
    bls = gerbls.pyFastBLS()
    bls.setup(phot,
              min_period,
              max_period,
              t_samp=t_samp,
              duration_mode='constant',
              durations=durations,
              downsample=True)
    bls.run(verbose=True)

    # Return the BLS spectrum
    blsa = gerbls.pyBLSAnalyzer(bls)
    return {'P': np.copy(blsa.P),
            'dchi2': np.copy(blsa.dchi2),
            't0': np.copy(blsa.t0),
            'dur': np.copy(blsa.dur),
            'mag0': np.copy(blsa.mag0),
            'dmag': np.copy(blsa.dmag),
            'snr': np.copy(blsa.snr)}

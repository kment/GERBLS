import gerbls
import numpy as np
import numpy.typing as npt


def run_bls(time: npt.ArrayLike,
            mag: npt.ArrayLike,
            err: npt.ArrayLike,
            min_period: float,
            max_period: float,
            durations: list = [],
            t_samp: float = 0.):
    """
    A basic convenience function to generate a BLS spectrum.
    The data must be evenly sampled in time to run the BLS,
    use `t_samp` to specify the cadence for any resampling.

    Parameters
    ----------
    time : npt.ArrayLike
        Array of observation timestamps.
    mag : npt.ArrayLike
        Array of observed fluxes.
    err : npt.ArrayLike
        Array of flux uncertainties for each observation.
    min_period : float
        Minimum BLS period to search.
    max_period : float
        Maximum BLS period to search.
    durations : list
        List of transit durations to test at each period.
    t_samp : float, optional
        Time sampling to bin the data before running the BLS.
        If 0 (default), the median time difference between observations is used.

    Returns
    -------
    dict
        Dictionary with BLS results:
        * `P` is the list of tested periods
        * `dchi2` is the BLS statistic (Delta chi squared) at each period
        * `t0` is the best-fit transit mid-point at each period
        * `dur` is the best-fit duration at each period
        * `mag0` is the best-fit flux baseline at each period
        * `dmag` is the best-fit transit depth at each period
    """

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
              durations=durations)
    bls.run(verbose=True)

    # Return the BLS spectrum
    blsa = gerbls.pyBLSAnalyzer(bls)
    return {'P': np.copy(blsa.P),
            'dchi2': np.copy(-blsa.dchi2),
            't0': np.copy(blsa.t0),
            'dur': np.copy(blsa.dur),
            'mag0': np.copy(blsa.mag0),
            'dmag': np.copy(blsa.dmag)}

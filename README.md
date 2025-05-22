# GERBLS

**GERBLS** (**G**reatly **E**xpedited **R**obust **B**ox **L**east **S**quares) is a lightweight fast-folding implementation of the BLS (Box Least Squares) algorithm. It is designed to facilitate transiting planet searches in photometric data via an easy setup and fast runtimes.

`GERBLS` can outperform popular brute-force BLS implementations such as `astropy.timeseries.BoxLeastSquares` by **over 10-20x** in runtime speed.

## Installation

Currently, `GERBLS` requires a Python version of 3.9 or above. Additional dependencies are `numpy` and a build-time dependency on `Cython`. These will be checked and/or installed automatically.

To install the latest version of `GERBLS`, run the following code:
```
pip install gerbls
```

If you encounter any issues while installing or using GERBLS, or would like to request a feature to be added to the code, please do not hesitate to [contact me](mailto:kxm821@psu.edu).

## Basic usage

A detrended light curve is required to run the BLS. You may use any of your favorite detrending algorithms; `scipy.signal.savgol_filter` is a relatively good option for long-term variability. A convenience function has been implemented for easy generation of the BLS spectrum, with the following parameters:
```
def run_bls(time: npt.ArrayLike, 
            mag: npt.ArrayLike, 
            err: npt.ArrayLike,
            min_period: float,
            max_period: float,
            t_samp: float = 0.):
    """
    A basic convenience function to generate a BLS spectrum.
    The data must be evenly sampled in time to run the BLS,
    use t_samp to specify the cadence for any resampling.

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
    t_samp : float, optional
        Time sampling to bin the data before running the BLS.
        If 0 (default), the median time difference between observations is used.

    Returns
    -------
    np.ndarray
        Array of tested periods.
    np.ndarray
        BLS statistic (Delta chi-squared) for each tested period.
    """
```

For example, running the following Python script generates a BLS spectrum for orbital periods between 0.4 and 10 days, where the light curve has been stored in arrays `time`, `mag`, and `err`:
```
from gerbls import run_bls
P, dchi2 = run_bls(time, mag, err, 0.4, 10)
```

The function `gerbls.run_bls` returns two arrays: a list of tested orbital periods `P` (in days) and a list of BLS statistics `dchi2` for each tested period. The searched periods are evenly spaced in frequency, and the spacing is set by the time sampling of the data. The BLS statistic `dchi2` ($\Delta\chi^2$) is the difference between the total $\chi^2$ parameters of a box-shaped model and a constant flux model fit to the data. In the case of pure Gaussian white noise, the signal-to-noise ratio of the fitted transit can be estimated as $\sqrt{\Delta\chi^2}$.

## Documentation

A full "Read The Docs" documentation page is currently in the works.

## Features in development

There are multiple additional features that are currently in various stages of development but need to be tested more thoroughly before they can be released publicly. These include:
- Various light curve detrending methods (Savitsky-Golay filter, Gaussian Process, etc.)
- Post-BLS limb-darkened transit model fitting
- Period-dependent bootstrap FAP calculation, which allows the significance of any potential transit to be evaluated (or alternative, an S/R threshold to be set) as a function of orbital period
- Additional tools to implement fake transit injection and recovery searches

## Acknowledgements

`GERBLS` includes some C code from the publicly available pulsar-searching [riptide](https://github.com/v-morello/riptide) package to implement the fast-folding algorithm.
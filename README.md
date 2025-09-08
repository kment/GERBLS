# GERBLS

**GERBLS** (**G**reatly **E**xpedited **R**obust **B**ox **L**east **S**quares) is a lightweight fast-folding implementation of the BLS (Box Least Squares) algorithm. It is designed to facilitate transiting planet searches in photometric data via an easy setup and fast runtimes.

`GERBLS` can outperform popular brute-force BLS implementations such as `astropy.timeseries.BoxLeastSquares` by **over 10-20x** in runtime speed.

## Installation

Refer to the [documentation](https://gerbls.readthedocs.io/en/stable/install.html) for installation instructions.

If you encounter any issues while installing or using GERBLS, or would like to request a feature to be added to the code, please do not hesitate to [contact me](mailto:kxm821@psu.edu).

## Documentation

A "Read The Docs" documentation page can be accessed [here](https://gerbls.readthedocs.io/).

## Features in development

There are multiple additional features that are currently in various stages of development but need to be tested more thoroughly before they can be released publicly. These include:
- Various light curve detrending methods (Savitsky-Golay filter, Gaussian Process, etc.)
- Post-BLS limb-darkened transit model fitting
- Period-dependent bootstrap FAP calculation, which allows the significance of any potential transit to be evaluated (or alternative, an S/R threshold to be set) as a function of orbital period
- Additional tools to implement fake transit injection and recovery searches

## Acknowledgements

`GERBLS` includes some C code from the publicly available pulsar-searching [riptide](https://github.com/v-morello/riptide) package to implement the fast-folding algorithm.

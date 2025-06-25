import gerbls
import numpy as np
from .exofunc import divide_into_chunks
from scipy.signal import savgol_filter
from typing import Union


def clean_savgol(phot: gerbls.pyDataContainer,
                 N_flares: int = 3,
                 sigma_clip: float = 5.,
                 window_length: float = 1.,
                 chunk_length: float = 0.,
                 verbose: bool = True):
    """
    Clean the data using a Savitsky-Golay filter.
    Also performs sigma clipping and flare rejection.
    Cubic splines will be used to interpolate over any masked data.

    Parameters
    ----------
    phot : gerbls.pyDataContainer
        Input data to be cleaned.
    N_flares : int, optional
        Number of iterations to detect flares, by default 3
    sigma_clip: float, optional
        Whether to remove outliers more than X sigma from the initial baseline, by default 5. To turn off, enter 0.
    window_length: float, optional
        Window length in time units (days) for the Savitsky-Golay filter, by default 1.
    chunk_length: float, optional
        If given, divide the data into chunks first where the time gaps in data are longer than this length,
        and apply the filter independently to each chunk, by default 0.
    verbose: bool, optional
        Whether to print any text/warnings, by default True

    Returns
    -------
    gerbls.pyDataContainer
        Cleaned data.
    np.array
        Mask corresponding to valid cleaned data.
    """

    # Mask corresponding to valid cleaned data
    cmask = np.ones(phot.size, dtype=bool)

    # Determine window size (=1 day) based on median time cadence
    filter_width = int(window_length / np.median(np.diff(phot.rjd)))
    if filter_width % 2 == 0:
        filter_width += 1

    # Subroutine that fits the filter to masked data
    def fit(mask):

        mag_data = np.copy(phot.mag)
        mag_fit = np.zeros(phot.size)

        # Use cubic splines to interpolate over unmasked data
        mag_data[~mask] = np.interp(phot.rjd[~mask], phot.rjd[mask], phot.mag[mask])

        # Divide data into chunks with gaps over 0.2 days
        for a, b in divide_into_chunks(phot.rjd, chunk_length):
            if b - a < filter_width:
                filter_width_ = b - a - (1 if (b - a) % 2 == 0 else 2)
                if verbose:
                    print(
                        f"Warning: Savgol window length shortened {filter_width} => {filter_width_}")
            else:
                filter_width_ = filter_width
            mag_fit[a:b] = savgol_filter(mag_data[a:b], filter_width_, 3)

        return mag_fit[mask]

    if sigma_clip > 0:
        mag0 = fit(cmask)
        cmask = (abs(phot.mag - mag0) <= phot.err * sigma_clip)

    for _ in range(N_flares):
        phot_ = phot.mask(cmask)
        cmask[cmask] = ~find_flares(phot_, fit(cmask))

    cphot = phot.mask(cmask)
    cphot.remove_planet(fit(cmask))

    return cphot, cmask


def find_flares(phot: gerbls.pyDataContainer,
                mag0: Union[np.ndarray, float] = 1):
    """
    Returns a mask with detected flares.
    A flare is at least 4 consecutive positive 2-sigma deviants.

    Parameters
    ----------
    phot : cinject.pyDataContainer
        Input data.
    mag0 : float or np.ndarray, optional
        Flux baseline to be subtracted before flare detection, by default 1.

    Returns
    -------
    np.ndarray[bool]
        Mask denoting flare data.
    """

    mag = phot.mag - mag0
    mask_ = (mag / phot.err >= 2)   # Mask for single points
    cs = np.cumsum(mask_, dtype=int)

    # Indices for first points of flares
    i_start = np.where(cs[4:] - cs[:-4] >= 4)[0] + 1

    # Mask for flares
    mask = np.zeros(phot.size, dtype=bool)

    j = 0
    for i in i_start:
        if i < j:
            continue
        j = i + 4
        while i > 0 and mag[i-1] > 0:
            i -= 1
        while j < phot.size and mag[j] > 0:
            j += 1
        mask[i:j] = True

    return mask

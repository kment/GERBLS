# cython: language_level = 3
# Various utility functions to be included in gerbls.pyx
# See gerbls.pyx for module imports

# Invert a boolean mask (array)
cdef bool_t[:] invert_mask(bool_t[:] mask):
    cdef Py_ssize_t i, n = len(mask)
    cdef bool_t[:] out = np.empty(n, dtype=np.bool_)
    
    for i in range(n):
        out[i] = not mask[i]
    
    return out

# General function to raise an ImportError about a missing dependency
def raise_import_error(str source_function, str missing_dep):
    raise ImportError(
        f"{source_function} requires {missing_dep}, which is an optional dependency. Please check "
        f"that {missing_dep.split('.')[0]} has been properly installed.")

def resample(pyDataContainer data, double t_samp):
    """
    Returns a copy of the input data, resampled to the specified time cadence.

    Parameters
    ----------
    data : gerbls.pyDataContainer
        Input data.
    t_samp : float
        Desired time cadence.
    
    Returns
    -------
    gerbls.pyDataContainer
    """
    out = pyDataContainer.from_ptr(resample_uniform(data.cPtr[0], t_samp).release(), True)
    return out
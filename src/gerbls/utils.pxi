# cython: language_level = 3
# Various utility functions to be included in gerbls.pyx
# See gerbls.pyx for module imports

# Loss (negative log-likelihood) and gradient vector for a sample drawn from a chi-squared
# distribution.
cdef tuple chi2_loss_and_grad_exp(double[:] x, double[:] sample):

    cdef Py_ssize_t i, N = sample.shape[0]
    cdef double lnx_scaled, x_scaled

    cdef double x0 = exp(x[0])
    cdef double x1 = exp(x[1])
    cdef double ix2 = exp(-x[2])

    cdef double loss = (x[2] + x0 * LN2 * 0.5 + lgamma(x0 * 0.5)) * N
    cdef double[:] grad = np.zeros(3, dtype=np.double)

    grad[0] = N * (LN2 + digamma(x0 * 0.5))
    grad[1] = -N * 0.5
    grad[2] = N * x0 * 0.5
    for i in range(N):
        x_scaled = (sample[i] - x1) * ix2
        if x_scaled <= 0:
            return (10.0**9, grad)
        lnx_scaled = log(x_scaled)
        loss += x_scaled * 0.5 - (x0 * 0.5 - 1) * lnx_scaled
        grad[0] -= lnx_scaled
        grad[1] += (x0 * 0.5 - 1) / x_scaled
        grad[2] -= x_scaled * 0.5
    grad[0] *= x0 * 0.5
    grad[1] *= x1 * ix2

    return (loss, grad)

def fit_chi2_dist(double[:] sample):
    """
    A wrapper to fit a chi-squared distribution to a sample using :func:`scipy.optimize.minimize`.
    Returns a :class:`scipy.optimize.OptimizeResult` object.
    To retrieve the three distribution parameters [DoF, loc, scale], use ``numpy.exp(res.x)``.

    Parameters
    ----------
    sample : ArrayLike
        Sample of values drawn from the chi-squared distribution.
    
    Returns
    -------
    res : scipy.optimize.OptimizeResult
        The optimization result (refer to SciPy docs for more specific documentation).
    """
    cdef double[:] x0 = np.array([log(np.mean(sample)), 0.01, 0.01])
    cdef res = minimize(chi2_loss_and_grad_exp, x0, method='BFGS', jac=True, args=(sample,))
    
    return res

# Invert a boolean mask (array)
cdef bool_t[:] invert_mask(bool_t[:] mask):
    cdef Py_ssize_t i, n = len(mask)
    cdef bool_t[:] out = np.empty(n, dtype=np.bool_)
    
    for i in range(n):
        out[i] = not mask[i]
    
    return out

def raise_import_error(str source_function, str missing_dep):
    """
    Raises a Python ImportError about a missing dependency.

    Parameters
    ----------
    source_function : str
        Name of the function that failed to run.
    missing_dep : str
        Name of the missing dependency module.
    """
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
    out : gerbls.pyDataContainer
        Resampled data.
    """
    out = pyDataContainer.from_ptr(resample_uniform(data.cPtr[0], t_samp).release(), True)
    return out
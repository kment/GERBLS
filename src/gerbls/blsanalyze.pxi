# cython: language_level = 3
# BLS analyzer to be included in gerbls.pyx
# See core.pyx for module imports

cdef class pyBLSAnalyzer:
    """
    pyBLSAnalyzer(model, noise=None, allow_noise_interp=True)

    BLS results analyzer.
    
    Parameters
    ----------
    model : gerbls.pyBLSModel
        BLS model generator that was used to generate the BLS spectrum.
    noise : gerbls.pyNoiseBLS, optional
        BLS noise model (if desired), by default None.
    allow_noise_interp : bool, optional
        Whether to interpolate the noise spectrum over the orbital period range. If False, `model`
        and `noise` must have the exact same number of tested periods. Has no effect if `noise` is
        None. By default True.
    
    
    .. property:: dchi2
        :type: numpy.ndarray

        Get the array of best-fit :math:`\Delta\chi^2` values for each tested period.
    
    .. property:: dmag
        :type: numpy.ndarray

        Get the array of best-fit transit depths for each tested period.

    .. property:: dur
        :type: numpy.ndarray

        Get the array of best-fit transit durations for each tested period.
    
    .. property:: f
        :type: numpy.ndarray

        Get the array of tested orbital frequencies (`= 1/period`).
    
    .. property:: mag0
        :type: numpy.ndarray

        Get the array of best-fit out-of-transit flux baselines for each tested period.
    
    .. property:: model
        :type: gerbls.pyBLSModel

        Get the BLS model tied to this object.
    
    .. property:: N_bins
        :type: numpy.ndarray

        Get the number of data points (bins) in the phase-folded light curve for each tested period.
    
    .. property:: noise
        :type: gerbls.pyNoiseBLS

        Get the noise BLS instance tied to this object.

    .. property:: P
        :type: numpy.ndarray

        Get the array of tested orbital periods.
    
    .. property:: snr
        :type: numpy.ndarray

        Get an estimated SNR at each period from :math:`\\textrm{SNR} \\approx \sqrt{\Delta\chi^2}`.
    
    .. property:: t0
        :type: numpy.ndarray

        Get the array of best-fit transit midpoint times for each tested period.
    """
    cdef readonly pyBLSModel model
    cdef readonly pyNoiseBLS noise
    cdef bool_t [:] _mask
    cdef bool_t noise_interp
    
    def __cinit__(self,
                  pyBLSModel model not None,
                  pyNoiseBLS noise = None,
                  bool_t allow_noise_interp = True):
        
        self.model = model
        self.noise = noise
        self.noise_interp = allow_noise_interp
        self.initialize_mask()
    
    @property
    def dchi2(self):
        self.model.assert_setup()
        return np.asarray(self.model.view_dchi2())
    
    @property
    def dmag(self):
        self.model.assert_setup()
        return np.asarray(self.model.view_dmag())
    
    #@property
    #def dmag_err(self):
    #    self.model.assert_setup()
    #    return np.asarray(self.model.view_dmag_err())
    #
    @property
    def dur(self):
        self.model.assert_setup()
        return np.asarray(self.model.view_dur())
    
    @property
    def f(self):
        self.model.assert_setup()
        return np.asarray(self.model.view_freq())
        
    @property
    def mag0(self):
        self.model.assert_setup()
        return np.asarray(self.model.view_mag0())
    
    @property
    def mask(self):
        return np.asarray(self._mask)

    @property
    def N_bins(self):
        self.model.assert_setup()
        return np.asarray(self.model.view_bins())
    
    @property
    def P(self):
        return self.f**-1
    
    @property
    def snr(self):
        return self.dchi2**0.5
    
    @property
    def t0(self):
        self.model.assert_setup()
        return np.asarray(self.model.view_t0())
    
    def fit_bls_trend(self, size_t window_length = 1001):
        """
        Fit a trendline to the SNR values using a median filter.

        Parameters
        ----------
        window_length : int
            Window length for the median filter.
        
        Returns
        -------
        np.ndarray
            Fitted SNR trend at each tested period.
        """
        return median_filter(self.snr, window_length, mode='reflect')
    
    def generate_models(self, N_models, double unmaskf = 0.005, bool_t use_SDE = False, **kwargs):
        """
        Identify the top BLS models (periods) in terms of highest :math:`\Delta\chi^2` values.

        Parameters
        ----------
        N_models : int
            Number of models to generate.
        unmaskf : float, optional
            The frequencies of any generated models must differ by at least this amount, by default
            0.005.
        use_SDE : bool, optional
            Whether to use the Signal Detection Efficiency (SDE) to identify peaks instead of the
            :math:`\Delta\chi^2` values, by default False.
        **kwargs
            Any keyword arguments are passed to :meth:`get_SDE`. Has no effect if `use_SDE` is
            False.
        
        Returns
        -------
        list
            List of :class:`gerbls.pyBLSResult` corresponding to the identified models.
        """
        # Perform input checks
        assert unmaskf > 0, "unmaskf must be positive."

        self.initialize_mask()
        return [self.generate_next_model(unmaskf, use_SDE, **kwargs) for _ in range(N_models)]
    
    def generate_next_model(self, double unmaskf = 0.005, bool_t use_SDE = False, **kwargs):
        """:meta private:"""
        cdef size_t index, mask_index
        
        if not self.mask.any():
            return None
        
        if use_SDE:
            mask_index = np.argmax(self.get_SDE(**kwargs)[self.mask])
        elif self.noise is None:
            mask_index = np.argmax(self.dchi2[self.mask])
        else:
            if self.noise_interp:
                noise_dchi2 = self.noise.dchi2(self.P)
                mask_index = np.argmax(self.dchi2[self.mask] - noise_dchi2[self.mask])
            else:
                assert self.model.N_freq == len(self.noise.dchi2_arr), "Noise model has the wrong number of periods."
                mask_index = np.argmax(self.dchi2[self.mask] - self.noise.dchi2_arr[self.mask])
        index = np.where(self.mask)[0][mask_index]
        
        # Returned frequencies must be some range apart
        self.unmask_freq(self.f[index], unmaskf)
        
        return pyBLSResult(self, index)
    
    def get_SDE(self, **kwargs):
        """
        Calculate the Signal Detection Efficiency (SDE) at each tested period.

        Parameters
        ----------
        **kwargs
            Passed to BLS trend calculation.
        
        Returns
        -------
        np.ndarray
            Array of SDE values at each tested period.
        """
        bls_trend = self.fit_bls_trend(**kwargs)
        bls_scatter = np.std(self.snr - bls_trend)
        return (self.snr - bls_trend) / bls_scatter
    
    cdef void initialize_mask(self):
        self.model.assert_setup()
        self._mask = np.ones(self.model.N_freq, dtype = np.bool_)
        # Ignore anti-transits
        self._mask *= (self.dmag > 0)
    
    # Mask out BLS frequencies less than df away from f_
    cpdef void unmask_freq(self, double f_, double df):
        """:meta private:"""
        self._mask *= (np.abs(self.f - f_) >= df)

cdef class pyBLSResult:
    """
    pyBLSResult(blsa, index)

    Fitted BLS model at a specific orbital period.

    Parameters
    ----------
    blsa : gerbls.pyBLSAnalyzer
        BLS analyzer object.
    index : int
        Index of the orbital period stored in the BLS analyzer.
    
    
    .. property:: dchi2
        :type: float
        
        Get the :math:`\Delta\chi^2` of the fitted model.

    .. property:: dmag
        :type: float
        
        Get the transit depth.

    .. property:: dur
        :type: float
        
        Get the transit duration.

    .. property:: mag0
        :type: float
        
        Get the out-of-transit flux baseline.

    .. property:: P
        :type: float
        
        Get the orbital period.
    
    .. property:: r
        :type: float
        
        Calculate the planet-to-star radius ratio.

    .. property:: snr
        :type: float
        
        Get an estimated SNR from :math:`\\textrm{SNR} \\approx \sqrt{\Delta\chi^2}`.
    
    .. property:: snr_from_dchi2
        :type: float
        
        Alias of :attr:`snr`.

    .. property:: t0
        :type: float
        
        Get the transit midpoint time.
    """
    cdef readonly double P
    cdef readonly double dchi2
    cdef readonly double mag0
    cdef readonly double dmag
    cdef readonly double t0
    cdef readonly double dur

    def __cinit__(self, pyBLSAnalyzer blsa not None, size_t index):

        assert blsa.model.N_freq > 0, "BLS model has no generated frequencies."
        assert 0 <= index < blsa.model.N_freq, "index out of bounds for the BLS model."

        self.P = blsa.P[index]
        self.dchi2 = blsa.dchi2[index]
        self.mag0 = blsa.mag0[index]
        self.dmag = blsa.dmag[index]
        self.t0 = (blsa.t0[index] + blsa.dur[index] / 2) % blsa.P[index]
        self.dur = blsa.dur[index]
    
    def __str__(self):
        return (
            f"pyBLSResult(P={self.P}, dchi2={self.dchi2}, mag0={self.mag0}, dmag={self.dmag}, "
            f"t0={self.t0}, dur={self.dur}, snr={self.snr_from_dchi2})"
        )
    
    @property
    def r(self):
        return (self.dmag / self.mag0)**0.5
    
    @property
    def snr(self):
        return self.dchi2**0.5
    
    @property
    def snr_from_dchi2(self):
        return self.snr
    
    def get_dmag_err(self, pyDataContainer phot not None):
        """
        Calculate the uncertainty in :attr:`dmag` (transit depth).

        Parameters
        ----------
        phot : gerbls.pyDataContainer
            Fitted data.
        
        Returns
        -------
        float
        """
        # Perform input checks
        assert phot.size > 0, "Data container cannot be empty."

        cdef bool_t[:] mask = self.get_transit_mask(phot.rjd)
        cdef double err_in = np.sum(phot.err[mask]**-2)**-0.5
        cdef double err_out = np.sum(phot.err[invert_mask(mask)]**-2)**-0.5
        return (err_in**2 + err_out**2)**0.5

    def get_SNR(self, pyDataContainer phot not None):
        """
        Calculate the transit SNR from uncertainty in :attr:`dmag`.

        Parameters
        ----------
        phot : gerbls.pyDataContainer
            Fitted data.
        
        Returns
        -------
        float
        """
        # Perform input checks
        assert phot.size > 0, "Data container cannot be empty."

        return self.dmag / self.get_dmag_err(phot)

    def get_transit_mask(self, double[:] t):
        """
        Determine which of the given input times are in-transit.

        Parameters
        ----------
        t : ArrayLike
            Array of observation times.
        
        Returns
        -------
        numpy.ndarray
            Boolean array with True values corresponding to in-transit data points.
        """
        # Perform input checks
        assert len(t) > 0, "No observation times given."

        return (abs((np.array(t) - self.t0 + self.P / 2) % self.P - self.P / 2) < self.dur / 2)
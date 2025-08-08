# cython: language_level = 3
# BLS model and analyzer to be included in gerbls.pyx
# See gerbls.pyx for module imports

cdef class pyBLSModel:
    """
    Base class for BLS model generators. Should not be created directly.
    """
    cdef BLSModel* cPtr
    cdef bool_t alloc           # Whether responsible for memory allocation
    
    def __cinit__(self):
        if type(self) is pyBLSModel:
            self.alloc = False
    
    def __dealloc__(self):
        if self.alloc and type(self) is pyBLSModel:
            del self.cPtr
    
    @property
    def freq(self):
        """np.ndarray: Array of tested frequencies."""
        return np.asarray(self.view_freq())
    
    #def get_max_duration(self, double P):
    #    return self.cPtr.get_max_duration(P)
    
    @property
    def N_freq(self):
        """int: Number of tested frequencies."""
        return self.cPtr.N_freq()
    
    def run(self, bool_t verbose = True):
        """
        Run the BLS generator.

        Parameters
        ----------
        verbose : bool
            Whether to print output to the console, by default True.
        
        Returns
        -------
        None
        """
        self.cPtr.run(verbose)

    cdef size_t [::1] view_bins(self):
        return <size_t [:self.N_freq]>self.cPtr.N_bins.data()
        
    cdef double [::1] view_dchi2(self):
        return <double [:self.N_freq]>self.cPtr.dchi2.data()
    
    cdef double [::1] view_dmag(self):
        return <double [:self.N_freq]>self.cPtr.chi2_dmag.data()
    
    cdef double [::1] view_dur(self):
        return <double [:self.N_freq]>self.cPtr.chi2_dt.data()
    
    cdef double [::1] view_freq(self):
        return <double [:self.N_freq]>self.cPtr.freq.data()
    
    cdef double [::1] view_mag0(self):
        return <double [:self.N_freq]>self.cPtr.chi2_mag0.data()
    
    cdef double [::1] view_t0(self):
        return <double [:self.N_freq]>self.cPtr.chi2_t0.data()

cdef class pyBruteForceBLS(pyBLSModel):
    """
    Brute-force (slow) BLS generator.
    :meth:`setup` should be used before :meth:`run`.
    """
    cdef BLSModel_bf* dPtr
    
    def __cinit__(self):
        self.alloc = False
    
    def __dealloc__(self):
        if self.alloc:
            del self.dPtr
    
    def setup(self,
              pyDataContainer data not None,
              double min_period,
              double max_period,
              pyTarget target = None,
              double dt_per_step = 0.,
              double t_bins = 0.,
              size_t N_bins_min = 0,
              str duration_mode = "",
              double min_duration_factor = 0.,
              double max_duration_factor = 0.):
        """
        Set up the BLS generation.

        Parameters
        ----------
        data : gerbls.pyDataContainer
            Input data.
        min_period : float
            Minimum searched orbital period.
        max_period : float
            Maximum searched orbital period.
        target : gerbls.pyTarget, optional
            Stellar parameters, by default None.
        dt_per_step : float, optional
            Period spacing will be calculated such that over the course of the entire time baseline
            of the data, any transit midtime will not be expected to shift by more than this value
            due to finite period spacing, by default 0.003.
        t_bins : float, optional
            Time cadence that phase-folded light curves will be binned to, by default 0.007.
        N_bins_min : int, optional
            Regardless of the value specified by `t_bins`, phase-folded light curves at each period
            are guaranteed to have at least this many bins in total, by default 100.
        duration_mode : {'constant', 'fractional', 'physical'}, optional
            Affects how the maximum tested transit duration is determined at each period, by default
            'fractional'.
        min_duration_factor : float, optional
            Affects the minimum searched transit duration at each period, by default 0.
        max_duration_factor : float, optional
            Affects the maximum searched transit duration at each period, by default 0.1.
        
        Returns
        -------
        None
        """
        cdef Target* targetPtr = (<Target *>NULL if target == None else target.cPtr)
        self.dPtr = new BLSModel_bf(data.cPtr[0],
                                    1/max_period,
                                    1/min_period,
                                    targetPtr, 
                                    dt_per_step,
                                    t_bins,
                                    N_bins_min,
                                    convert_duration_mode(duration_mode),
                                    min_duration_factor,
                                    max_duration_factor)
        self.cPtr = self.dPtr
        self.alloc = True
    
    # Setup with a pre-defined frequency array
    def setup_from_freq(self,
                        pyDataContainer data not None,
                        double[:] freq_,
                        pyTarget target = None,
                        double t_bins = 0.,
                        size_t N_bins_min = 0,
                        str duration_mode = "",
                        double min_duration_factor = 0.,
                        double max_duration_factor = 0.):
        """
        Set up the BLS generation with a predefined array of orbital frequencies.

        Parameters
        ----------
        data : gerbls.pyDataContainer
            Input data.
        freq : ArrayLike
            Array of orbital frequencies (`= 1/period`) to test.
        target : gerbls.pyTarget, optional
            Stellar parameters, by default None.
        t_bins : float, optional
            Time cadence that phase-folded light curves will be binned to, by default 0.007.
        N_bins_min : int, optional
            Regardless of the value specified by `t_bins`, phase-folded light curves at each period
            are guaranteed to have at least this many bins in total, by default 100.
        duration_mode : {'constant', 'fractional', 'physical'}, optional
            Affects how the maximum tested transit duration is determined at each period, by default
            'fractional'.
        min_duration_factor : float, optional
            Affects the minimum searched transit duration at each period, by default 0.
        max_duration_factor : float, optional
            Affects the maximum searched transit duration at each period, by default 0.1.
        
        Returns
        -------
        None
        """
        cdef Target* targetPtr = (<Target *>NULL if target == None else target.cPtr)
        self.dPtr = new BLSModel_bf(data.cPtr[0],
                                    list(freq_),
                                    targetPtr,
                                    t_bins,
                                    N_bins_min,
                                    convert_duration_mode(duration_mode),
                                    min_duration_factor,
                                    max_duration_factor)
        self.cPtr = self.dPtr
        self.alloc = True

cdef class pyFastBLS(pyBLSModel):
    """
    Fast-folding BLS generator.
    :meth:`setup` should be used before :meth:`run`.
    """
    cdef BLSModel_FFA* dPtr
    
    def __cinit__(self):
        self.alloc = False
    
    def __dealloc__(self):
        if self.alloc:
            del self.dPtr
    
    @property
    def rdata(self):
        """
        gerbls.pyDataContainer: Resampled data generated by :meth:`run`, with a time sampling given
        by :attr:`t_samp`.
        """
        return pyDataContainer.from_ptr(self.dPtr.rdata.get(), False)
    
    def run_double(self, bool_t verbose = True):
        """
        Run the BLS generator with all output results in `double` precision.

        Parameters
        ----------
        verbose : bool
            Whether to print output to the console, by default True.
        
        Returns
        -------
        None
        """
        self.dPtr.run_double(verbose)
    
    def setup(self,
              pyDataContainer data,
              double min_period,
              double max_period,
              pyTarget target = None,
              double t_samp = 0.,
              bool_t verbose = True,
              str duration_mode = "",
              vector[double] durations = [],
              double min_duration_factor = 0.,
              double max_duration_factor = 0.,
              bool_t downsample = False,
              double downsample_invpower = 3.,
              double downsample_threshold = 1.1):
        """
        Set up the BLS generation.

        Parameters
        ----------
        data : gerbls.pyDataContainer
            Input data.
        min_period : float
            Minimum searched orbital period.
        max_period : float
            Maximum searched orbital period.
        target : gerbls.pyTarget, optional
            Stellar parameters, by default None.
        t_samp : float, optional
            Desired initial time sampling of the data, by default 0. Overwrites the value in
            :meth:`t_samp`. If 0, the median time cadence of the input data will be used instead.
        verbose : bool, optional
            Whether to print output to the console, by default True.
        duration_mode : {'constant', 'fractional', 'physical'}, optional
            Affects how the maximum tested transit duration is determined at each period, by default
            'fractional'.
        durations : list, optional
            If given, use a specific list of duration factors instead of a range.
        min_duration_factor : float, optional
            Affects the minimum searched transit duration at each period, by default 0. Has no
            effect if `durations` is given.
        max_duration_factor : float, optional
            Affects the maximum searched transit duration at each period, by default 0.1. Has no
            effect if `durations` is given.
        downsample : bool, optional
            Whether to automatically downsample the data at longer periods, by default False.
        downsample_invpower : float, optional
            Affects the rate of downsampling, by default 3.
        downsample_threshold : float, optional
            Affects the threshold that triggers downsampling, by default 1.1.
        
        Returns
        -------
        None
        """
        cdef Target* targetPtr = (<Target *>NULL if target == None else target.cPtr)
        if t_samp == 0:
            t_samp = np.median(np.diff(data.rjd))
            if verbose:
                print(
                    f"BLS time sampling set to the median cadence of input data: "
                    f"{t_samp*24*60:.2f} minutes.",
                    flush=True)
        self.dPtr = new BLSModel_FFA(data.cPtr[0],
                                     1./max_period,
                                     1./min_period,
                                     targetPtr,
                                     convert_duration_mode(duration_mode),
                                     (&durations if durations.size() else NULL),
                                     min_duration_factor,
                                     max_duration_factor,
                                     t_samp,
                                     downsample,
                                     downsample_invpower,
                                     downsample_threshold)
        self.cPtr = self.dPtr
        self.alloc = True
    
    @property
    def t_samp(self):
        """
        float: Desired (initial) time sampling during BLS generation.
        Value can be set manually.
        """
        return self.dPtr.t_samp
    @t_samp.setter
    def t_samp(self, double value):
        self.dPtr.t_samp = value
    
    @property
    def time_spent(self):
        """np.ndarray: The runtime spent at each orbital period during :meth:`run`."""
        return np.asarray(<double [:self.dPtr.time_spent.size()]>self.dPtr.time_spent.data())

cdef class pyBLSAnalyzer:
    """
    BLS results analyzer.
    
    Parameters
    ----------
    model : gerbls.pyBLSModel
        BLS model generator.
    """
    cdef size_t [:] _bins
    cdef double [:] _dchi2
    cdef double [:] _dmag
    cdef double [:] _dur
    cdef double [:] _freq
    cdef double [:] _mag0
    cdef bool_t [:] _mask
    cdef double [:] _t0
    cdef int N_freq
    cdef double t_samp
    
    def __cinit__(self, pyBLSModel model):
        self._bins = model.view_bins()
        self._dchi2 = model.view_dchi2()
        self._dmag = model.view_dmag()
        self._dur = model.view_dur()
        self._freq = model.view_freq()
        self._mag0 = model.view_mag0()
        self._t0 = model.view_t0()
        self.N_freq = model.N_freq
        self.t_samp = (model.t_samp if hasattr(model, "t_samp") else 0)
        self.initialize_mask()
    
    @property
    def dchi2(self):
        """np.ndarray: Array of best-fit :math:`\Delta\chi^2` values for each tested period."""
        return np.asarray(self._dchi2)
    
    @property
    def dmag(self):
        """np.ndarray: Array of best-fit transit depths for each tested period."""
        return np.asarray(self._dmag)
    
    @property
    def dur(self):
        """np.ndarray: Array of best-fit transit durations for each tested period."""
        return np.asarray(self._dur)
    
    @property
    def f(self):
        """np.ndarray: Array of tested orbital frequencies (`= 1/period`)."""
        return np.asarray(self._freq)
    
    cdef void initialize_mask(self):
        self._mask = np.ones(self.N_freq, dtype = np.bool_)
        # Ignore anti-transits
        self._mask *= (self.dmag > 0)
        
    @property
    def mag0(self):
        """np.ndarray: Array of best-fit out-of-transit flux baselines for each tested period."""
        return np.asarray(self._mag0)
    
    @property
    def mask(self):
        return np.asarray(self._mask)

    @property
    def N_bins(self):
        return np.asarray(self._bins)
    
    @property
    def P(self):
        """np.ndarray: Array of tested orbital periods."""
        return self.f**-1
    
    @property
    def t0(self):
        """np.ndarray: Array of best-fit transit midpoint times for each tested period."""
        return np.asarray(self._t0)
    
    def generate_models(self, N_models, double unmaskf = 0.005):
        """
        Identify the top BLS models (periods) in terms of highest :math:`\Delta\chi^2` values.

        Parameters
        ----------
        N_models : int
            Number of models to generate.
        unmaskf : float, optional
            The frequencies of any generated models must differ by at least this amount, by default
            0.005.
        
        Returns
        -------
        List of `gerbls.pyBLSResult`
        """
        self.initialize_mask()
        return [self.generate_next_model(unmaskf) for _ in range(N_models)]
    
    def generate_next_model(self, double unmaskf = 0.005):
        
        if not self.mask.any():
            return None
        
        cdef size_t mask_index = np.argmax(-self.dchi2[self.mask])
        cdef size_t index = np.where(self.mask)[0][mask_index]
        
        # Returned frequencies must be some range apart
        self.unmask_freq(self._freq[index], unmaskf)
        
        return pyBLSResult(self, index)
    
    # Mask out BLS frequencies less than df away from f_
    cpdef void unmask_freq(self, double f_, double df):
        self._mask *= (np.abs(self.f - f_) >= df)

cdef class pyBLSResult:
    """
    Fitted BLS model at a specific orbital period.

    Parameters
    ----------
    blsa : gerbls.pyBLSAnalyzer
        BLS analyzer object.
    index : int
        Index of the orbital period stored in the BLS analyzer.
    
    Attributes
    ----------
    dchi2 : float
        :math:`\Delta\chi^2` of the fitted model.
    dmag : float
        Transit depth.
    dur : float
        Transit duration.
    mag0 : float
        Out-of-transit flux baseline.
    P : float
        Orbital period.
    t0 : float
        Transit midpoint time.
    """
    cdef readonly double P
    cdef readonly double dchi2
    cdef readonly double mag0
    cdef readonly double dmag
    cdef readonly double t0
    cdef readonly double dur

    def __cinit__(self, pyBLSAnalyzer blsa, size_t index):
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
        """float: Calculate the planet-to-star radius ratio."""
        return (self.dmag / self.mag0)**0.5
    
    @property
    def snr_from_dchi2(self):
        """float: An initial estimate of the SNR from :math:`\\textrm{SNR} \\approx \sqrt{\Delta\chi^2}`."""
        return ((-self.dchi2)**0.5 if self.dchi2 <= 0 else -np.inf)
    
    def get_dmag_err(self, pyDataContainer phot):
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
        cdef bool_t[:] mask = self.get_transit_mask(phot.rjd)
        cdef double err_in = np.sum(phot.err[mask]**-2)**-0.5
        cdef double err_out = np.sum(phot.err[invert_mask(mask)]**-2)**-0.5
        return (err_in**2 + err_out**2)**0.5

    def get_SNR(self, pyDataContainer phot):
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
        np.ndarray
            Boolean array with True values corresponding to in-transit data points.
        """
        return (abs((np.array(t) - self.t0 + self.P / 2) % self.P - self.P / 2) < self.dur / 2)

cdef int convert_duration_mode(str duration_mode):
    """
    Converts a string representation of a duration mode to its integer counterpart.
    """
    cdef dict allowed_duration_modes = {'': 0,
                                        'constant': 1,
                                        'fractional': 2,
                                        'physical': 3}
    assert (
        duration_mode in allowed_duration_modes
        ), f"duration_mode must be one of: {allowed_duration_modes.keys()}"
    
    return allowed_duration_modes[duration_mode]
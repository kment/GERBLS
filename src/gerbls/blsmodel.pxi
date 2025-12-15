# cython: language_level = 3
# BLS models to be included in gerbls.pyx
# See core.pyx for module imports

# Allowed duration modes for BLS models
cdef dict allowed_duration_modes = {'': DurationMode.None,
                                    'constant': DurationMode.Constant,
                                    'fractional': DurationMode.Fractional,
                                    'physical': DurationMode.Physical}

cdef class pyBLSModel:
    """
    Base class for BLS model generators. Should not be created directly.

    .. property:: duration_mode
        :type: str

        Get the current `duration_mode` value, which affects how the maximum tested transit duration
        is determined at each orbital period.
    
    .. property:: durations
        :type: list

        Get the list of tested transit durations, if specified during model setup.

    .. property:: freq
        :type: numpy.ndarray

        Get the array of tested frequencies.
    
    .. property:: max_duration_factor
        :type: float

        Get the `max_duration_factor` that affects the maximum tested transit duration.
    
    .. property:: max_period
        :type: float

        Get the maximum tested orbital period.
    
    .. property:: min_duration_factor
        :type: float

        Get the `min_duration_factor` that affects the minimum tested transit duration.
    
    .. property:: min_period
        :type: float

        Get the minimum tested orbital period.
    
    .. property:: N_freq
        :type: int

        Get the number of tested frequencies.
    
    .. property:: target
        :type: pyTarget | None

        Get a non-mutable reference to the stellar parameters object, if specified.
    """
    cdef BLSModel* cPtr
    cdef bool_t alloc           # Whether responsible for memory allocation
    
    def __cinit__(self):
        if type(self) is pyBLSModel:
            self.alloc = False
    
    def __dealloc__(self):
        if self.alloc and type(self) is pyBLSModel:
            del self.cPtr
    
    cdef void assert_setup(self):
        assert self.cPtr is not NULL, "Model needs to be set up first."
    
    def calculate_N_freq(self):
        """
        Calculate the number of tested frequencies without running the model.

        Returns
        -------
        int
        """
        self.assert_setup()
        return self.cPtr.calculate_N_freq()

    @property
    def duration_mode(self):
        self.assert_setup()
        return next((k for k, v in allowed_duration_modes.items() if v == self.cPtr.duration_mode), 
                    "")

    cdef DurationMode duration_mode_enum(self, str duration_mode):
        """Convert a string representation of a duration mode to its enum counterpart."""
        assert (
            duration_mode in allowed_duration_modes
            ), f"duration_mode must be one of: {allowed_duration_modes.keys()}"
        
        return allowed_duration_modes[duration_mode]

    @property
    def durations(self):
        self.assert_setup()
        return list(<double [:self.cPtr.durations.size()]>self.cPtr.durations.data())

    @property
    def freq(self):
        self.assert_setup()
        return np.asarray(self.view_freq())
    
    @property
    def max_duration_factor(self):
        self.assert_setup()
        return 1./self.cPtr.max_duration_factor
    
    @property
    def max_period(self):
        self.assert_setup()
        return 1./self.cPtr.f_min
    
    @property
    def min_duration_factor(self):
        self.assert_setup()
        return 1./self.cPtr.min_duration_factor

    @property
    def min_period(self):
        self.assert_setup()
        return 1./self.cPtr.f_max
    
    @property
    def N_freq(self):
        self.assert_setup()
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
        self.assert_setup()
        self.cPtr.run(verbose, True)

    @property
    def target(self):
        self.assert_setup()
        if self.cPtr.target is NULL:
            return None
        else:
            return pyTarget.from_const_ptr(self.cPtr.target)

    cdef size_t [::1] view_bins(self):
        return <size_t [:self.N_freq]>self.cPtr.N_bins.data()
        
    cdef double [::1] view_dchi2(self):
        return <double [:self.N_freq]>self.cPtr.dchi2.data()
    
    cdef double [::1] view_dmag(self):
        return <double [:self.N_freq]>self.cPtr.chi2_dmag.data()
    
    #cdef double [::1] view_dmag_err(self):
    #    return <double [:self.N_freq]>self.cPtr.chi2_dmag_err.data()
    
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
    :meth:`setup` should be used before :meth:`~pyBLSModel.run`.
    Refer to the base class for additional properties and methods.
    """
    cdef BLSModel_bf* dPtr
    
    def __cinit__(self):
        self.alloc = False
    
    def __dealloc__(self):
        if self.alloc:
            del self.dPtr
    
    cdef pyBruteForceBLS duplicate(self):
        cdef pyBruteForceBLS new_model = pyBruteForceBLS.__new__(pyBruteForceBLS)
        new_model.alloc = True
        new_model.dPtr = <BLSModel_bf*>self.cPtr.duplicate().release()
        new_model.cPtr = new_model.dPtr
        return new_model

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

        .. caution:: Only references to ``data`` and ``target`` are stored. Crashes or unexpected
            results may occur if the referenced objects get deallocated.

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
        # Perform input checks
        assert data.size > 0, "data cannot be empty."
        assert max_period > min_period > 0, "Invalid min and/or max period."
        assert dt_per_step >= 0, "dt_per_step cannot be negative."
        assert t_bins >= 0, "t_bins cannot be negative."
        assert N_bins_min >= 0, "N_bins_min cannot be negative."
        assert min_duration_factor >= 0, "min_duration_factor cannot be negative."
        assert max_duration_factor >= 0, "max_duration_factor cannot be negative."

        cdef const Target* targetPtr = (<const Target *>NULL if target == None else target.cptr)
        self.dPtr = new BLSModel_bf(data.cPtr[0],
                                    1/max_period,
                                    1/min_period,
                                    targetPtr, 
                                    dt_per_step,
                                    t_bins,
                                    N_bins_min,
                                    self.duration_mode_enum(duration_mode),
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

        .. caution:: Only references to ``data`` and ``target`` are stored. Crashes or unexpected
            results may occur if the referenced objects get deallocated.

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
        # Perform input checks
        assert data.size > 0, "data cannot be empty."
        assert len(freq_) > 0, "List of frequencies cannot be empty."
        assert t_bins >= 0, "t_bins cannot be negative."
        assert N_bins_min >= 0, "N_bins_min cannot be negative."
        assert min_duration_factor >= 0, "min_duration_factor cannot be negative."
        assert max_duration_factor >= 0, "max_duration_factor cannot be negative."

        cdef const Target* targetPtr = (<const Target *>NULL if target == None else target.cptr)
        self.dPtr = new BLSModel_bf(data.cPtr[0],
                                    list(freq_),
                                    targetPtr,
                                    t_bins,
                                    N_bins_min,
                                    self.duration_mode_enum(duration_mode),
                                    min_duration_factor,
                                    max_duration_factor)
        self.cPtr = self.dPtr
        self.alloc = True

cdef class pyFastBLS(pyBLSModel):
    """
    Fast-folding BLS generator.
    :meth:`setup` should be used before :meth:`~pyBLSModel.run`.
    Refer to the base class for additional properties and methods.

    .. property:: rdata
        :type: gerbls.pyDataContainer

        Get the resampled data generated by :meth:`~pyBLSModel.run`, with a time sampling given by
        :attr:`t_samp`.
    
    .. property:: t_samp
        :type: float

        Get or set the desired (initial) time sampling during BLS generation.
    
    .. property:: time_spent
        :type: numpy.ndarray

        Get the runtime spent at each orbital period during :meth:`~pyBLSModel.run`.
    """
    cdef BLSModel_FFA* dPtr
    
    def __cinit__(self):
        self.alloc = False
    
    def __dealloc__(self):
        if self.alloc:
            del self.dPtr
    
    cpdef pyFastBLS duplicate(self):
        cdef pyFastBLS new_model = pyFastBLS.__new__(pyFastBLS)
        new_model.alloc = True
        new_model.dPtr = <BLSModel_FFA*>self.cPtr.duplicate().release()
        new_model.cPtr = new_model.dPtr
        return new_model
    
    @property
    def rdata(self):
        self.assert_setup()
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
        self.dPtr.run_double(verbose, True)
    
    def setup(self,
              pyDataContainer data not None,
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

        .. caution:: Only references to ``data`` and ``target`` are stored. Crashes or unexpected
            results may occur if the referenced objects get deallocated.

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
            :attr:`t_samp`. If 0, the median time cadence of the input data will be used instead.
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
        # Perform input checks
        assert data.size > 0, "data cannot be empty."
        assert max_period > min_period > 0, "Invalid min and/or max period."
        assert t_samp >= 0, "t_samp cannot be negative."
        assert min_duration_factor >= 0, "min_duration_factor cannot be negative."
        assert max_duration_factor >= 0, "max_duration_factor cannot be negative."

        cdef const Target* targetPtr = (<const Target *>NULL if target == None else target.cptr)
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
                                     self.duration_mode_enum(duration_mode),
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
        self.assert_setup()
        return self.dPtr.t_samp
    @t_samp.setter
    def t_samp(self, double value):
        self.assert_setup()
        assert value > 0, "t_samp must be positive."
        self.dPtr.t_samp = value
    
    @property
    def time_spent(self):
        self.assert_setup()
        return np.asarray(<double [:self.dPtr.time_spent.size()]>self.dPtr.time_spent.data())
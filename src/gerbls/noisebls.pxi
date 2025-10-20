# cython: language_level = 3
# Noise BLS model to be included in gerbls.pyx
# See core.pyx for module imports

# Allowed selection modes for noise BLS
cdef dict allowed_selection_modes = {'': NoiseMode.None,
                                     'fit': NoiseMode.FittedChi2Dist,
                                     'max': NoiseMode.MaximumDChi2}

cdef class pyNoiseBLS:
    """
    Noise BLS model.

    Parameters
    ----------
    model : pyBLSModel
        BLS model used to generate the noise BLS. An internal copy will be created.
    

    .. property:: dchi2_arr
        :type: float
        
        Get the :math:`\Delta\chi^2` array of the generated noise BLS model.

    .. property:: f
        :type: numpy.ndarray

        Get the array of tested orbital frequencies (`= 1/period`).
    
    .. property:: N_freq
        :type: int

        Get the number of tested frequencies.

    .. property:: N_sim
        :type: int

        Get the number of simulations that were used to generate the noise BLS.
    
    .. property:: P
        :type: numpy.ndarray

        Get the array of tested orbital periods.

    .. property:: selection_mode
        :type: str

        Get the current `selection_mode` value, which affects how the noise BLS value is calculated
        from simulated BLS spectra.
    """
    cdef NoiseBLS* cPtr
    cdef bool_t alloc           # Whether responsible for memory allocation
    
    def __cinit__(self):
        self.alloc = False
    
    def __dealloc__(self):
        if self.alloc:
            del self.cPtr
    
    def __init__(self, pyBLSModel model not None):
        model.assert_setup()
        self.cPtr = new NoiseBLS(model.cPtr[0])
        self.alloc = True
    
    cdef void assert_generated(self):
        assert self.cPtr.dchi2.size() > 0, "Noise BLS has not been generated."

    def dchi2(self, object P):
        """
        Evaluate the :math:`\Delta\chi^2` of the generated noise BLS at the given period(s). Uses
        linear interpolation over the tested period range.

        Parameters
        ----------
        P : float or ArrayLike
            The period(s) at which to evaluate :math:`\Delta\chi^2`.
        
        Returns
        -------
        float or numpy.ndarray
        """
        return np.interp(P, self.P, self.dchi2_arr)

    @property
    def dchi2_arr(self):
        self.assert_generated()
        return np.asarray(self.view_dchi2())
    
    @property
    def f(self):
        return np.asarray(self.view_freq())
    
    def generate(self, size_t N_sim, str selection_mode = "", bool_t verbose = True):
        """
        Generate the noise BLS spectrum.

        Parameters
        ----------
        N_sim : int
            Number of simulations.
        selection_mode : {'fit', 'max'}, optional
            Affects how the noise BLS value is calculated from simulated BLS spectra, by default
            'max'.
        verbose : bool, optional
            Whether to print output to the console, by default True.
        
        Returns
        -------
        None
        """
        # Perform input checks
        assert N_sim > 0, "N_sim must be positive."

        cdef vector[double] dchi2_ = self.cPtr.generate(N_sim,
                                                        self.selection_mode_enum(selection_mode),
                                                        verbose)
        cdef double* dchi2_ptr = dchi2_.data()
        cdef size_t N_freq = dchi2_.size() // N_sim  # Number of tested periods
        cdef size_t i, j
        cdef double[:] chi2_param, memview
        for i in range(N_freq):
            if self.selection_mode == 'max':
                self.cPtr.dchi2[i] = dchi2_[i * N_sim]
                for j in range(N_sim):
                    if dchi2_[i * N_sim + j] > self.cPtr.dchi2[i]:
                        self.cPtr.dchi2[i] = dchi2_[i * N_sim + j]
            elif self.selection_mode == 'fit':
                memview = <double[:N_sim]> (dchi2_ptr + i * N_sim)
                chi2_param = np.exp(fit_chi2_dist(memview).x)
                self.cPtr.dchi2[i] = chi2.isf(1./N_freq, *chi2_param)
    
    def generate_binned(self,
                        size_t N_sim,
                        double bin_width,
                        str selection_mode = "",
                        double bin_width_prec = 0.1,
                        bool_t verbose = True):
        """
        Generate the noise BLS spectrum.

        Parameters
        ----------
        N_sim : int
            Number of simulations.
        selection_mode : {'fit', 'max'}, optional
            Affects how the noise BLS value is calculated from simulated BLS spectra, by default
            'max'.
        verbose : bool, optional
            Whether to print output to the console, by default True.
        
        Returns
        -------
        None
        """
        # Perform input checks
        assert N_sim > 0, "N_sim must be positive."
        assert bin_width > 0, "bin_width must be positive."
        
        cdef size_t i
        cdef double bin_width_
        cdef object rng = np.random.default_rng()
        cdef size_t N_freq = self.cPtr.model.get().calculate_N_freq()

        dchi2_arr = np.zeros((N_sim, N_freq), dtype=np.double)
        cdef double[:, ::1] dchi2_view = dchi2_arr
        cdef double[::1, :] dchi2_T = dchi2_view.T
        cdef double[::1] dchi2_model
        cdef double[:] chi2_param

        data = pyDataContainer.from_ptr(self.cPtr.model.get().data)
        self.cPtr.dchi2.resize(N_freq)

        rjd_ = np.zeros(data.size, dtype=np.double)
        mag_ = np.zeros(data.size, dtype=np.double)
        cdef double[::1] rjd_view = rjd_
        cdef double[::1] mag_view = mag_
        sdata = pyDataContainer()
        sdata.assign(data.view_rjd(), mag_view, data.view_err())
        self.cPtr.model.get().data = sdata.cPtr

        for i in range(N_sim):

            bin_width_ = ((rng.random() * 2 - 1) * bin_width_prec + 1) * bin_width

            bin_ids = (data.rjd / bin_width_ + rng.random()).astype(int)
            bin_ids_unique = np.unique(bin_ids)
            bins_order = rng.permutation(bin_ids_unique)
            rjd_view[:] = data.view_rjd()

            for j, j_bin in enumerate(bin_ids_unique):
                mask_ = (bin_ids == j_bin)
                rjd_[mask_] += (bins_order[j] - j_bin) * bin_width
            
            mag_[:] = data.mag[rjd_.argsort()]

            self.cPtr.model.get().run(False, False)

            assert self.cPtr.model.get().N_freq() == N_freq, "Wrong number of model periods."

            dchi2_model = <double[:N_freq]> self.cPtr.model.get().dchi2.data()
            dchi2_view[i, :] = dchi2_model
        
        if verbose:
            print("Calculating noise spectrum...")

        for i in range(N_freq):
            if selection_mode == 'max':
                self.cPtr.dchi2[i] = np.max(dchi2_T[i, :])
            elif selection_mode == 'fit':
                chi2_param = np.exp(fit_chi2_dist(dchi2_T[i, :]).x)
                self.cPtr.dchi2[i] = chi2.isf(1./N_freq, *chi2_param)

    @property
    def N_freq(self):
        return self.cPtr.model.get().N_freq()

    @property
    def N_sim(self):
        return self.cPtr.N_sim
    
    @property
    def P(self):
        return self.f**-1

    @property
    def selection_mode(self):
        return next(
            (k for k, v in allowed_selection_modes.items() if v == self.cPtr.selection_mode), "")

    cdef NoiseMode selection_mode_enum(self, str selection_mode):
        """Convert a string representation of a selection mode to its enum counterpart."""
        assert (
            selection_mode in allowed_selection_modes
            ), f"selection_mode must be one of: {allowed_selection_modes.keys()}"
        
        return allowed_selection_modes[selection_mode]

    cdef double [::1] view_dchi2(self):
        return <double [:self.cPtr.dchi2.size()]>self.cPtr.dchi2.data()
    
    cdef double [::1] view_freq(self):
        return <double [:self.N_freq]>self.cPtr.model.get().freq.data()
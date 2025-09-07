# cython: language_level = 3
# Noise BLS model to be included in gerbls.pyx
# See gerbls.pyx for module imports

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
    

    .. property:: dchi2
        :type: float
        
        Get the :math:`\Delta\chi^2` array of the generated noise BLS model.
    
    .. property:: N_sim
        :type: int

        Get the number of simulations that were used to generate the noise BLS.
    
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
    
    @property
    def dchi2(self):
        return np.asarray(self.view_dchi2())
    
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
        cdef size_t i
        cdef size_t j
        cdef double[:] memview, res
        for i in range(N_freq):
            if self.selection_mode == 'max':
                for j in range(N_sim):
                    if dchi2_[i * N_sim + j] < self.cPtr.dchi2[i]:
                        self.cPtr.dchi2[i] = dchi2_[i * N_sim + j]
            elif self.selection_mode == 'fit':
                memview = <double[:N_sim]> (dchi2_ptr + i * N_sim)
                self.cPtr.dchi2[i] = chi2.isf(1./N_freq, *(fit_chi2_dist(memview)))

    @property
    def N_sim(self):
        return self.cPtr.N_sim

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
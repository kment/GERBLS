# cython: language_level = 3
# BLS model and analyzer to be included in gerbls.pyx

cdef class pyBLSModel:
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
        return np.asarray(self.view_freq())
    
    @property
    def N_freq(self):
        return self.cPtr.N_freq()
    
    def run(self, bool_t verbose=False):
        self.cPtr.run(verbose)
        
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

cdef class pyFastBLS(pyBLSModel):
    cdef BLSModel_FFA* dPtr
    
    def __cinit__(self):
        self.alloc = False
    
    def __dealloc__(self):
        if self.alloc:
            del self.dPtr
    
    def create(self, pyDataContainer data, double f_min, double f_max, pyTarget target):
        self.dPtr = new BLSModel_FFA(data.cPtr[0], f_min, f_max, target.cPtr)
        self.cPtr = self.dPtr
        self.alloc = True
    
    @property
    def dchi2(self):
        snr = np.asarray(<double [:self.dPtr.snr.size()]>self.dPtr.snr.data())
        N_widths = self.dPtr.widths.size()
        return snr.reshape((int(len(snr) / N_widths), N_widths))
    
    @property
    def foldbins(self):
        return np.asarray(<size_t [:self.dPtr.foldbins.size()]>self.dPtr.foldbins.data())
    
    @property
    def periods(self):
        return np.asarray(<double [:self.dPtr.periods.size()]>self.dPtr.periods.data())
    
    @property
    def rdata(self):
        return pyDataContainer.from_ptr(self.dPtr.rdata.get(), False)
    
    def run_double(self, bool_t verbose=True):
        self.dPtr.run_double(verbose)
    
    @property
    def t_samp(self):
        return self.dPtr.t_samp
    @t_samp.setter
    def t_samp(self, double value):
        self.dPtr.t_samp = value
        
    @property
    def t_widths(self):
        return np.asarray(<size_t [:self.dPtr.widths.size()]>self.dPtr.widths.data())
    
    @property
    def t0(self):
        t0 = np.asarray(<size_t [:self.dPtr.t0.size()]>self.dPtr.t0.data())
        N_widths = self.dPtr.widths.size()
        return t0.reshape((int(len(t0) / N_widths), N_widths))

cdef class pyBLSAnalyzer:
    cdef double [:] _dchi2
    cdef double [:] _dmag
    cdef double [:] _dur
    cdef double [:] _freq
    cdef double [:] _mag0
    cdef bool_t [:] _mask
    cdef double [:] _t0
    cdef int N_freq
    cdef double t_samp
    cdef readonly pyDataContainer data
    
    def __cinit__(self, pyBLSModel model, pyDataContainer data_):
        self._dchi2 = model.view_dchi2()
        self._dmag = model.view_dmag()
        self._dur = model.view_dur()
        self._freq = model.view_freq()
        self._mag0 = model.view_mag0()
        self._t0 = model.view_t0()
        self.N_freq = model.N_freq
        self.t_samp = (model.t_samp if hasattr(model, "t_samp") else 0)
        self.data = data_
        self.initialize_mask()
    
    @property
    def dchi2(self):
        return np.asarray(self._dchi2)
    
    @property
    def dmag(self):
        return np.asarray(self._dmag)
    
    @property
    def dur(self):
        return np.asarray(self._dur)
    
    @property
    def f(self):
        return np.asarray(self._freq)
    
    cdef void initialize_mask(self):
        self._mask = np.ones(self.N_freq, dtype=np.bool_)
        # Ignore anti-transits
        self._mask *= (self.dmag > 0)
        
    @property
    def mag0(self):
        return np.asarray(self._mag0)
    
    @property
    def mask(self):
        return np.asarray(self._mask)
    
    @property
    def P(self):
        return self.f**-1
    
    @property
    def t0(self):
        return np.asarray(self._t0)
    
    cpdef (double, double) generate_next_model_quick(self, double unmaskf=0.005):
        
        if not self.mask.any():
            return 0, 0
        
        cdef size_t mask_index = np.argmax(-self.dchi2[self.mask])
        cdef size_t index = np.where(self.mask)[0][mask_index]
        cdef double P = 1 / self._freq[index]
        
        # Returned frequencies must be some range apart
        self.unmask_freq(self._freq[index], unmaskf)
        
        return P, -self._dchi2[index]
    
    # Mask out BLS frequencies less than df away from f_
    cpdef void unmask_freq(self, double f_, double df):
        self._mask *= (np.abs(self.f - f_) >= df)
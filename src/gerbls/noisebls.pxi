# cython: language_level = 3
# Noise BLS model to be included in gerbls.pyx
# See gerbls.pyx for module imports

from scipy.stats import chi2

cdef class pyNoiseBLS:
    cdef NoiseBLS* cPtr
    cdef bool_t alloc           # Whether responsible for memory allocation
    
    def __cinit__(self):
        self.alloc = False
    
    def __dealloc__(self):
        if self.alloc:
            del self.cPtr
    
    def __init__(self, pyBLSModel model not None):
        self.cPtr = new NoiseBLS(model.cPtr[0])
        self.alloc = True
    
    @property
    def dchi2(self):
        return np.asarray(self.view_dchi2())
    
    def generate(self, size_t N_sim, int selection_mode = 0, bool_t verbose = True):
        cdef vector[double] dchi2_ = self.cPtr.generate(N_sim, selection_mode, verbose)
        cdef double* dchi2_ptr = dchi2_.data()
        cdef size_t N_freq = dchi2_.size() // N_sim  # Number of tested periods
        cdef size_t i
        cdef size_t j
        cdef double[:] memview, res
        for i in range(N_freq):
            if selection_mode == 0:
                for j in range(N_sim):
                    if dchi2_[i * N_sim + j] < self.cPtr.dchi2[i]:
                        self.cPtr.dchi2[i] = dchi2_[i * N_sim + j]
            elif selection_mode == 1:
                memview = <double[:N_sim]> (dchi2_ptr + i * N_sim)
                self.cPtr.dchi2[i] = chi2.isf(1./N_freq, *(fit_chi2_dist(memview)))

    cdef double [::1] view_dchi2(self):
        return <double [:self.cPtr.dchi2.size()]>self.cPtr.dchi2.data()
# cython: language_level = 3str
# distutils: language = c++

from gerbls.cgerbls cimport *
from libc.stdlib cimport srand
from libc.time cimport time as ctime
from libcpp cimport bool as bool_t
import numpy as np
cimport numpy as np

include "blsmodel.pxi"

# Initialize random number generator
srand(ctime(NULL))

cdef class pyDataContainer:
    cdef DataContainer* cPtr
    cdef bool_t alloc           # Whether responsible for cPtr memory allocation
    
    def __cinit__(self):
        self.alloc = False
    
    def __dealloc__(self):
        if self.alloc:
            del self.cPtr
    
    # Allocate memory for and take ownership of cPtr
    cdef void allocate(self):
        self.cPtr = new DataContainer()
        self.alloc = True
    
    # Assign data without copying
    cpdef void assign(self, double[::1] rjd, double[::1] mag, double[::1] err):
        cdef Py_ssize_t i
        if self.cPtr is NULL:
            self.allocate()
        self.cPtr.set(&rjd[0], &mag[0], &err[0], rjd.shape[0])
    
    def clean(self, double P_rot = 0, int N_flares = 3):
        cdef bool_t[::1] mask = np.zeros(self.cPtr.size, dtype = bool)
        data = pyDataContainer()
        data.cPtr = self.cPtr.clean(P_rot, &mask[0], N_flares).release()
        data.alloc = True
        return data, np.asarray(mask)
    
    def clean_hw(self, double hw, int N_flares = 3):
        cdef bool_t[::1] mask = np.zeros(self.cPtr.size, dtype = bool)
        data = pyDataContainer()
        data.cPtr = self.cPtr.clean_hw(hw, &mask[0], N_flares).release()
        data.alloc = True
        return data, np.asarray(mask)
    
    @property
    def err(self):
        return np.asarray(self.view_err())
    
    def find_flares(self, double[:] mag0 = None):
        from clean import find_flares
        return find_flares(self, mag0)
    
    @staticmethod
    cdef pyDataContainer from_ptr(DataContainer* ptr, bool_t _alloc = False):
        cdef pyDataContainer data = pyDataContainer()
        data.cPtr = ptr
        data.alloc = _alloc
        return data
    
    @property
    def mag(self):
        return np.asarray(self.view_mag())
    
    def mask(self, bool_t[:] mask):
        data = pyDataContainer()
        data.store_sec(self.sec[mask], 
                       self.rjd[mask], 
                       self.mag[mask], 
                       self.err[mask],
                       convert_to_flux = False)
        return data
    
    def phase_folded(self, double P_rot, double t_extend):
        data = pyDataContainer()
        data.cPtr = self.cPtr.phase_folded(P_rot, t_extend).release()
        data.alloc = True
        return data
    
    # Divide out a planetary signal
    def remove_planet(self, double[:] lc_model):
        cdef double [::1] lc_model_ = np.ascontiguousarray(lc_model)
        cdef size_t i
        for i in range(self.cPtr.size):
            self.cPtr.mag[i] /= lc_model[i]
    
    # Rescale error bars in each sector such that std(mag)=median(err)
    def rescale_err(self):
        cdef size_t i
        cdef int[::1] sec = self.view_sec()
        for sec_ in self.sectors:
            mask = (self.sec == sec_)
            factor = np.std(self.mag[mask]) / np.median(self.err[mask])
            for i in range(self.cPtr.size):
                if sec[i] == sec_:
                    self.cPtr.err[i] *= factor
    
    @property
    def rjd(self):
        return np.asarray(self.view_rjd())
    
    def running_median(self, double hwidth):
        return np.asarray(self.cPtr.running_median(hwidth))
    
    def running_median_eval(self, double hwidth, double[:] t_eval):
        cdef double[::1] t_ = np.ascontiguousarray(t_eval)
        return np.asarray(self.cPtr.running_median_eval(hwidth, &t_[0], t_.shape[0]))
    
    def running_median_per(self, double hwidth, double P_rot):
        return np.asarray(self.cPtr.running_median_per(hwidth, P_rot))
    
    @property
    def sec(self):
        return np.asarray(self.view_sec())
    
    @property
    def sectors(self):
        return np.unique(self.view_sec())
    
    @property
    def size(self):
        return self.cPtr.size
    
    #def splfit(self, double P_rot, int M=50):
    #    return np.asarray(self.cPtr.splfit(P_rot, M))
    
    #def splfit_eval(self, double[:] t_eval, int M=50):
    #    cdef double[::1] t_ = np.ascontiguousarray(t_eval)
    #    return np.asarray(self.cPtr.splfit_eval(M, &t_[0], t_.shape[0]))
    
    def split_by_sector(self):
        cdef dict data = {}
        for sector in self.sectors:
            mask = (self.sec == sector)
            data[sector] = pyDataContainer()
            data[sector].store_sec(self.sec[mask], self.rjd[mask], self.mag[mask],
                                   self.err[mask], False)
        return data
    
    # Store data by making a copy
    def store(self, double[:] rjd_, double[:] mag_, double[:] err_, bool_t convert_to_flux = False):
        cdef Py_ssize_t i
        cdef double[::1] rjd = np.ascontiguousarray(rjd_)
        cdef double[::1] mag = np.ascontiguousarray(mag_)
        cdef double[::1] err = np.ascontiguousarray(err_)
        if self.cPtr is NULL:
            self.allocate()
        self.cPtr.store(&rjd[0], &mag[0], &err[0], rjd.shape[0])
        if convert_to_flux:
            for i in range(rjd.shape[0]):
                self.cPtr.mag[i] = 10.0**(-0.4 * self.cPtr.mag[i])
                self.cPtr.err[i] = 0.4 * np.log(10.0) * self.cPtr.mag[i] * self.cPtr.err[i]
    
    def store_sec(self,
                  int[:] sec_,
                  double[:] rjd_,
                  double[:] mag_,
                  double[:] err_, 
                  bool_t convert_to_flux = False):
        cdef Py_ssize_t i
        self.store(rjd_, mag_, err_, convert_to_flux)
        for i in range(len(sec_)):
            self.cPtr.sec[i] = sec_[i]
            
    def store_sec_d(self,
                    double[:] sec_,
                    double[:] rjd_,
                    double[:] mag_,
                    double[:] err_, 
                    bool_t convert_to_flux = False):
        self.store_sec(np.asarray(sec_, dtype=np.int32), rjd_, mag_, err_, convert_to_flux)
    
    cdef double [::1] view_err(self):
        return <double [:self.cPtr.size]>self.cPtr.err
    
    cdef double [::1] view_mag(self):
        return <double [:self.cPtr.size]>self.cPtr.mag
    
    cdef double [::1] view_rjd(self):
        return <double [:self.cPtr.size]>self.cPtr.rjd
    
    cdef int[::1] view_sec(self):
        if self.cPtr.sec is NULL:
            return None
        else:
            return <int [:self.cPtr.size]>self.cPtr.sec

cdef class pyTarget:
    cdef Target* cPtr
    
    def __cinit__(self):
        self.cPtr = new Target()
    
    def __dealloc__(self):
        del self.cPtr
    
    @property
    def L(self):
        return self.cPtr.L
    @L.setter
    def L(self, double L_):
        self.cPtr.L = L_
    
    @property
    def L_comp(self):
        return self.cPtr.L_comp
    @L_comp.setter
    def L_comp(self, double L_):
        self.cPtr.L_comp = L_
    
    @property
    def logg(self):
        return self.cPtr.logg()
    
    @property
    def M(self):
        return self.cPtr.M
    @M.setter
    def M(self, double M_):
        self.cPtr.M = M_
    
    @property
    def Prot(self):
        return self.cPtr.P_rot
    @Prot.setter
    def Prot(self, double Prot_):
        self.cPtr.P_rot = Prot_
    
    @property
    def Prot2(self):
        return self.cPtr.P_rot2
    @Prot2.setter
    def Prot2(self, double Prot2_):
        self.cPtr.P_rot2 = Prot2_
    
    @property
    def R(self):
        return self.cPtr.R
    @R.setter
    def R(self, double R_):
        self.cPtr.R = R_
        
    @property
    def u(self):
        return [self.cPtr.u1, self.cPtr.u2]
    
    @property
    def Teff(self):
        return self.cPtr.Teff()

def resample(pyDataContainer data, double t_samp):
    out = pyDataContainer.from_ptr(resample_uniform(data.cPtr[0], t_samp).release(), True)
    return out
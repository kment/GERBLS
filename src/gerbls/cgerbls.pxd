# distutils: language = c++
# Cython bindings for imearth

from libcpp cimport bool as bool_t
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
#from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

cdef extern from "cpp/ffafunc.hpp":
    unique_ptr[DataContainer] resample_uniform(DataContainer, double)

cdef extern from "cpp/model.hpp":
    cdef cppclass BLSModel:
        int max_duration_mode
        double max_duration_factor
        vector[double] chi2_dmag
        vector[double] chi2_dt
        vector[double] chi2_mag0
        vector[double] chi2_t0
        vector[double] dchi2
        vector[double] freq
        vector[size_t] N_bins
        
        size_t N_freq()
        void run(bool_t)
    
    cdef cppclass BLSModel_bf(BLSModel):
        BLSModel_bf(
            DataContainer, double, double, Target*, double, double, size_t, int, double, double)
        BLSModel_bf(DataContainer, vector[double], Target*, double, size_t, int, double, double)
    
    cdef cppclass BLSModel_FFA(BLSModel):
        double t_samp
        unique_ptr[DataContainer] rdata
        
        vector[double] periods
        vector[size_t] widths
        vector[size_t] foldbins
        vector[double] snr
        vector[size_t] t0
        BLSModel_FFA(DataContainer, double, double, Target*, int, double, double)
        
        void run_double(bool_t)

cdef extern from "cpp/physfunc.hpp":
    double gdraw(double, double, double)
    double get_aR_ratio(double, double, double)
    double get_inc(double, double, double, double)
    void get_phase_range(double, double*, double*)
    double get_transit_dur(double, double, double, double)
    double grand(double, double)

cdef extern from "cpp/structure.hpp":
    cdef cppclass DataContainer:
        int* sec
        double* rjd
        double* mag
        double* err
        size_t size
        
        unique_ptr[DataContainer] clean(double, bool_t*, int)
        unique_ptr[DataContainer] clean_hw(double, bool_t*, int)
        unique_ptr[bool_t] find_flares(const double*)
        unique_ptr[bool_t] find_flares()
        void imprint(double*, size_t)
        unique_ptr[DataContainer] phase_folded(double, double)
        vector[double] running_median(double)
        vector[double] running_median_eval(double, double*, size_t)
        vector[double] running_median_per(double, double)
        void set(double*, double*, double*, size_t)
        #vector[double] splfit(double, int)
        #vector[double] splfit_eval(int, double*, size_t)
        #unordered_map[int, unique_ptr[DataContainer]] split_by_sector()
        void store(double*, double*, double*, size_t)

    cdef cppclass Target:
        double M
        double R
        double L
        double u1
        double u2
        double L_comp
        double P_rot
        double P_rot2
        
        double logg()
        double Teff()
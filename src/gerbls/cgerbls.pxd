# distutils: language = c++
# Cython bindings for imearth

from libcpp cimport bool as bool_t
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
#from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

cdef extern from "cpp/ffafunc.hpp":
    size_t periodogram_length(size_t, double, double, double, bool, double, double)
    unique_ptr[DataContainer] resample_uniform(DataContainer, double, double)

cdef extern from "cpp/model.hpp":
    cdef enum class DurationMode:
        None
        Constant
        Fractional
        Physical
    
    cdef enum class NoiseMode:
        None
        FittedChi2Dist
        MaximumDChi2

    cdef cppclass BLSModel:
        DataContainer* data
        const Target* target
        double f_min
        double f_max
        double min_duration_factor
        double max_duration_factor
        DurationMode duration_mode
        vector[double] chi2_dmag
        #vector[double] chi2_dmag_err
        vector[double] chi2_dt
        vector[double] chi2_mag0
        vector[double] chi2_t0
        vector[double] dchi2
        vector[double] durations
        vector[double] freq
        vector[size_t] N_bins
        
        size_t calculate_N_freq()
        unique_ptr[BLSModel] duplicate() const
        size_t N_freq()
        void run(bool_t, bool_t)
    
    cdef cppclass BLSModel_bf(BLSModel):
        double t_bins
        BLSModel_bf(DataContainer,
                    double,
                    double,
                    const Target*,
                    double,
                    double,
                    size_t,
                    DurationMode,
                    double,
                    double)
        BLSModel_bf(DataContainer,
                    vector[double],
                    const Target*,
                    double,
                    size_t,
                    DurationMode,
                    double,
                    double)
    
    cdef cppclass BLSModel_FFA(BLSModel):
        double t_samp
        unique_ptr[DataContainer] rdata
        vector[double] time_spent
        
        BLSModel_FFA(DataContainer,
                     double,
                     double,
                     const Target*,
                     DurationMode,
                     const vector[double]*,
                     double,
                     double,
                     double,
                     bool_t,
                     double,
                     double)
        
        void run_double(bool_t, bool_t)
    
    cdef cppclass NoiseBLS:
        NoiseMode selection_mode
        size_t N_sim
        unique_ptr[BLSModel] model
        vector[double] dchi2

        NoiseBLS(BLSModel)

        vector[double] generate(size_t, NoiseMode, bool_t)

cdef extern from "cpp/physfunc.hpp":
    double get_aR_ratio(double, double, double)
    double get_inc(double, double, double, double)
    void get_phase_range(double, double*, double*)
    double get_transit_dur(double, double, double, double)

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
        
        double logg() const
        double Teff() const
    
    cdef cppclass Vector2D_double "Vector2D<double>":
        vector[double] data
        size_t cols
        size_t rows
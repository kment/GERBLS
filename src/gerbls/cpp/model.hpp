/*
 * model.hpp
 *
 *  Created on: Aug 20, 2017
 *      Author: Kristo Ment
 */

#ifndef MODEL_HPP_
#define MODEL_HPP_

#include "structure.hpp"
#include <fstream>
#include <tuple>
#include <unordered_map>

// Forward declarations
template <typename T> struct BLSResult;

// BLS model (base class)
struct BLSModel {

    // Settings
    double f_min = 0.025;             // Minimum search frequency
    double f_max = 5;                 // Maximum search frequency
    int duration_mode = 2;            // Affects tested transit durations
    double min_duration_factor = 0;   // Affects get_min_duration()
    double max_duration_factor = 0.1; // Affects get_max_duration()
    std::vector<double> durations;    // List of transit durations to test

    // Pointer to associated data
    DataContainer *data = nullptr;

    // Pointer to associated target
    const Target *target = nullptr;

    // Array to store tested frequencies
    std::vector<double> freq;

    // Constructor and destructor
    BLSModel(DataContainer &data_ref,
             double f_min = 0.,
             double f_max = 0.,
             const Target *targetPtr = nullptr,
             int duration_mode = 0,
             const std::vector<double> *durations = nullptr,
             double min_duration_factor = 0.,
             double max_duration_factor = 0.);
    virtual ~BLSModel() = default;

    bool explicit_durations() const; // Whether to use explicit durations (as opposed to range)
    std::tuple<double, double>
        get_duration_limits(double P) const; // Min & max transit duration to test at a given period
    size_t N_freq();                         // Get number of frequencies
    void set_widths(double P,
                    double tau,
                    std::vector<size_t> &widths) const; // Set transit widths for period P

    // Virtual functions to be overwritten
    virtual void run(bool verbose);

    // Required results for each tested frequency
    std::vector<double> dchi2, chi2_mag0, chi2_dmag, chi2_t0, chi2_dt;

    // Number of phase-folded data bins at each tested frequency
    std::vector<size_t> N_bins;
};

// BLS model (brute force)
struct BLSModel_bf : public BLSModel {

    // Grid search ranges and steps
    double dt_per_step = 0.003; // Maximum orbital shift between frequencies in days
    double t_bins = 0.007;      // Time bin width in days
    size_t N_bins_min = 100;    // Minimum number of bins

    // Arrays to store best chi2 values for each tested frequency
    std::vector<double> chi2, chi2r;

    // Constructors
    BLSModel_bf(DataContainer &data_ref,
                double f_min = 0.,
                double f_max = 0.,
                const Target *targetPtr = nullptr,
                double dt_per_step = 0.,
                double t_bins = 0.,
                size_t N_bins_min = 0,
                int duration_mode = 0,
                double min_duration_factor = 0.,
                double max_duration_factor = 0.);
    BLSModel_bf(DataContainer &data_ref,
                const std::vector<double> &freq,
                const Target *targetPtr = nullptr,
                double t_bins = 0.,
                size_t N_bins_min = 0,
                int duration_mode = 0,
                double min_duration_factor = 0.,
                double max_duration_factor = 0.);

    // Methods to overwrite parent virtual functions
    void run(bool verbose = true);

    // Private methods
private:
    void initialize(double t_bins, size_t N_bins_min);
};

// BLS model (FFA)
struct BLSModel_FFA : public BLSModel {

    // Settings
    bool downsample = false; // Automatic downsampling for shorter periods
    double ds_invpower = 3.;
    double ds_threshold = 1.1;     // Downsample when the max transit duration
                                   // drops by this fraction
    size_t N_bins_transit_min = 1; // Minimum number of bins per transit
    double t_samp = 2. / 60 / 24;  // Uniform cadence to resample data to

    // Pointer to the resampled data
    std::unique_ptr<DataContainer> rdata;

    // Time spent evaluating dchi2 at each tested period
    std::vector<double> time_spent;

    // Constructor
    BLSModel_FFA(DataContainer &data_ref,
                 double f_min = 0.,
                 double f_max = 0.,
                 const Target *targetPtr = nullptr,
                 int duration_mode = 0,
                 const std::vector<double> *durations = nullptr,
                 double min_duration_factor = 0.,
                 double max_duration_factor = 0.,
                 double t_samp = 0.,
                 bool downsample = false,
                 double ds_invpower = 0.,
                 double ds_threshold = 0.,
                 size_t N_bins_transit_min = 0);

    // Methods
    template <typename T> void process_results(std::vector<BLSResult<T>> &results);
    void run(bool verbose = true);
    void run_double(bool verbose = true);
    template <typename T> void run_prec(bool verbose = true);
};

#endif /* MODEL_HPP_ */

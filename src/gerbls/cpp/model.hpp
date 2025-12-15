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

// Scoped enumerations
enum class DurationMode { None, Constant, Fractional, Physical };
enum class NoiseMode { None, FittedChi2Dist, MaximumDChi2 };

// BLS model (base class)
struct BLSModel {

    // Settings
    double f_min = 0.025;                                  // Minimum search frequency
    double f_max = 5;                                      // Maximum search frequency
    DurationMode duration_mode = DurationMode::Fractional; // Affects tested transit durations
    double min_duration_factor = 0;                        // Affects get_min_duration()
    double max_duration_factor = 0.1;                      // Affects get_max_duration()
    std::vector<double> durations;                         // List of transit durations to test

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
             DurationMode duration_mode = DurationMode::None,
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
    virtual size_t calculate_N_freq(); // Calculate N_freq() without running the model
    virtual std::unique_ptr<BLSModel> duplicate() const; // Create a new similar object
    virtual void run(bool verbose, bool full_results);

    // Required results for each tested frequency
    std::vector<double> dchi2, chi2_mag0, chi2_dmag, chi2_t0, chi2_dt;
    //std::vector<double> chi2_dmag_err;

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
                DurationMode duration_mode = DurationMode::None,
                double min_duration_factor = 0.,
                double max_duration_factor = 0.);
    BLSModel_bf(DataContainer &data_ref,
                const std::vector<double> &freq,
                const Target *targetPtr = nullptr,
                double t_bins = 0.,
                size_t N_bins_min = 0,
                DurationMode duration_mode = DurationMode::None,
                double min_duration_factor = 0.,
                double max_duration_factor = 0.);

    // Methods to overwrite parent virtual functions
    size_t calculate_N_freq() override;
    std::unique_ptr<BLSModel> duplicate() const override;
    void run(bool verbose = true, bool full_results = true) override;

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
                 DurationMode duration_mode = DurationMode::None,
                 const std::vector<double> *durations = nullptr,
                 double min_duration_factor = 0.,
                 double max_duration_factor = 0.,
                 double t_samp = 0.,
                 bool downsample = false,
                 double ds_invpower = 0.,
                 double ds_threshold = 0.,
                 size_t N_bins_transit_min = 0);

    // Methods to overwrite parent virtual functions
    size_t calculate_N_freq() override;
    std::unique_ptr<BLSModel> duplicate() const override;
    void run(bool verbose = true, bool full_results = true) override;

    // Additional methods
    template <typename T>
    void process_results(std::vector<BLSResult<T>> &results, bool full_results = true);
    void run_double(bool verbose = true, bool full_results = true);
    template <typename T> void run_prec(bool verbose = true, bool full_results = true);
};

// Noise BLS
struct NoiseBLS {

    // Settings
    size_t N_sim = 0;
    NoiseMode selection_mode = NoiseMode::MaximumDChi2;

    // Pointer to the BLS model
    std::unique_ptr<BLSModel> model;

    // Noise BLS spectrum
    std::vector<double> dchi2;

    // Constructor
    NoiseBLS(BLSModel &model);

    // Other methods
    std::vector<double> generate(size_t N_sim, NoiseMode selection_mode, bool verbose = true);
};

#endif /* MODEL_HPP_ */
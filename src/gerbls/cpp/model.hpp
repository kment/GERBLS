/*
 * model.hpp
 *
 *  Created on: Aug 20, 2017
 *      Author: Kristo Ment
 */

#ifndef MODEL_HPP_
#define MODEL_HPP_

#include "ffafunc.hpp"
#include "structure.hpp"
#include <fstream>
#include <tuple>
#include <unordered_map>

// BLS model (base class)
struct BLSModel {
    // Settings
    double f_min = 0.025;             // Minimum search frequency
    double f_max = 5;                 // Maximum search frequency
    int duration_mode = 2;            // Affects tested transit durations
    double min_duration_factor = 0;   // Affects get_min_duration()
    double max_duration_factor = 0.1; // Affects get_max_duration()

    // Pointer to associated data
    DataContainer *data = nullptr;

    // Pointer to associated target
    const Target *target = nullptr;

    // Array to store tested frequencies
    std::vector<double> freq;

    // Constructor and destructor
    BLSModel(DataContainer &data_ref,
             double f_min,
             double f_max,
             const Target *targetPtr,
             int duration_mode,
             double min_duration_factor,
             double max_duration_factor);
    virtual ~BLSModel() = default;

    std::tuple<double, double>
        get_duration_limits(double P); // Min and max transit duration to test at a given period
    size_t N_freq();                   // Get number of frequencies

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
                double f_min,
                double f_max,
                const Target *targetPtr,
                double dt_per_step,
                double t_bins,
                size_t N_bins_min,
                int duration_mode,
                double min_duration_factor,
                double max_duration_factor);
    BLSModel_bf(DataContainer &data_ref,
                const std::vector<double> &freq,
                const Target *targetPtr,
                double t_bins,
                size_t N_bins_min,
                int duration_mode,
                double min_duration_factor,
                double max_duration_factor);

    // Methods to overwrite parent virtual functions
    void run(bool verbose = true);

    // Private methods
private:
    void initialize(double t_bins, size_t N_bins_min);
};

// BLS model (FFA)
struct BLSModel_FFA : public BLSModel {
    // Settings
    double t_samp = 2. / 60 / 24; // Uniform cadence to resample data to

    // Pointer to the resampled data
    std::unique_ptr<DataContainer> rdata;

    // TEMPORARY FFA results
    std::vector<double> periods;
    std::vector<size_t> widths;
    std::vector<size_t> foldbins;
    std::vector<double> snr;
    std::vector<size_t> t0;

    // Inherit constructor from parent
    using BLSModel::BLSModel;

    // Methods
    template <typename T> void process_results(std::vector<BLSResult<T>> &results);
    void run(bool verbose = true);
    void run_double(bool verbose = true);
    template <typename T> void run_prec(bool verbose = true);
};

#endif /* MODEL_HPP_ */

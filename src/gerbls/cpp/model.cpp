/*
 * model.cpp
 *
 *  Created on: Aug 20, 2017
 *      Author: Kristo Ment
 */

#include "ffafunc.hpp"
// #include "interpolation.h" 	// ALGLIB dependency
#include "model.hpp"
#include "physfunc.hpp"
#include "structure.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>

// Constructor
// If any numeric value is 0 then the default value is used
BLSModel::BLSModel(DataContainer &data_ref,
                   double f_min,
                   double f_max,
                   const Target *targetPtr,
                   DurationMode duration_mode,
                   const std::vector<double> *durations,
                   double min_duration_factor,
                   double max_duration_factor)
{
    data = &data_ref;
    target = targetPtr;

    if (f_min > 0)
        this->f_min = f_min;
    if (f_max > 0)
        this->f_max = f_max;
    if (duration_mode != DurationMode::None)
        this->duration_mode = duration_mode;
    if (durations != nullptr)
        this->durations = *durations;
    if (min_duration_factor > 0)
        this->min_duration_factor = min_duration_factor;
    if (max_duration_factor > 0)
        this->max_duration_factor = max_duration_factor;
}

size_t BLSModel::calculate_N_freq()
{
    throw std::logic_error("calculate_N_freq() is not defined for an object of type " +
                           std::string(typeid(*this).name()));
}

std::unique_ptr<BLSModel> BLSModel::duplicate() const
{
    throw std::logic_error("duplicate() is not defined for an object of type " +
                           std::string(typeid(*this).name()));
}

// Whether to use transit durations defined by durations, as opposed to a range defined by
// min_duration_factor and max_duration_factor
bool BLSModel::explicit_durations() const
{
    return !durations.empty();
}

// Get the minimum and maximum tested transit duration at a given period P
std::tuple<double, double> BLSModel::get_duration_limits(double P) const
{
    double min_duration, max_duration;

    switch (duration_mode) {
    // Constant duration limits
    case DurationMode::Constant:
        min_duration = min_duration_factor;
        max_duration = max_duration_factor;
        break;

    // Duration limits proportional to the orbital period
    case DurationMode::Fractional:
        min_duration = min_duration_factor * P;
        max_duration = max_duration_factor * P;
        break;

    // Duration limits proportional to the predicted physical transit duration
    case DurationMode::Physical:
        if (target == nullptr) {
            throw std::runtime_error("Target must not be null with duration_mode == Physical.");
            return std::make_tuple(0, 0);
        }
        else {
            double transit_dur = get_transit_dur(P, target->M, target->R, 0);
            min_duration = min_duration_factor * transit_dur;
            max_duration = max_duration_factor * transit_dur;
        }
        break;

    // Invalid duration code
    default:
        throw std::runtime_error("BLSModel::get_duration_limits() called with invalid duration_mode"
                                 " (int value " +
                                 std::to_string(static_cast<int>(duration_mode)) +
                                 ")");
        return std::make_tuple(0, 0);
    }

    return std::make_tuple(min_duration, max_duration);
}

// Get number of frequencies
size_t BLSModel::N_freq()
{
    return freq.size();
}

void BLSModel::run(bool verbose, bool full_results)
{
    throw std::logic_error("run() is not defined for an object of type " +
                           std::string(typeid(*this).name()));
}

// Set the searched transit widths (durations in bins) for a given period P
// tau is the time sampling of data bins
// CAUTION! Assumes that widths has the same size as this->durations (not checked)
void BLSModel::set_widths(double P, double tau, std::vector<size_t> &widths) const
{

    switch (duration_mode) {
    // Constant duration limits
    case DurationMode::Constant:
        for (size_t i = 0; i < durations.size(); i++) {
            widths[i] = round(durations[i] / tau);
        }
        break;

    // Duration limits proportional to the orbital period
    case DurationMode::Fractional:
        for (size_t i = 0; i < durations.size(); i++) {
            widths[i] = round(durations[i] * P / tau);
        }
        break;

    // Duration limits proportional to the predicted physical transit duration
    case DurationMode::Physical:
        if (target == nullptr) {
            throw std::runtime_error("Target must not be null with duration_mode == Physical.");
        }
        else {
            const double transit_dur = get_transit_dur(P, target->M, target->R, 0);
            for (size_t i = 0; i < durations.size(); i++) {
                widths[i] = round(durations[i] * transit_dur / tau);
            }
        }
        break;

    // Invalid duration code
    default:
        throw std::runtime_error("BLSModel::set_widths() called with invalid duration_mode"
                                 " (int value " +
                                 std::to_string(static_cast<int>(duration_mode)) +
                                 ")");
    }
}

// Constructor
// If any numeric value is 0 then the default value is used
// If no target is given, use default values
BLSModel_bf::BLSModel_bf(DataContainer &data_ref,
                         double f_min,
                         double f_max,
                         const Target *targetPtr,
                         double dt_per_step,
                         double t_bins,
                         size_t N_bins_min,
                         DurationMode duration_mode,
                         double min_duration_factor,
                         double max_duration_factor) :
    BLSModel(data_ref,
             f_min,
             f_max,
             targetPtr,
             duration_mode,
             nullptr,
             min_duration_factor,
             max_duration_factor)
{
    // Override numeric values if given
    if (dt_per_step > 0)
        this->dt_per_step = dt_per_step;

    // Multiplicative frequency step (round to 8 decimal places)
    double df = this->dt_per_step / data_ref.get_time_range();
    df = round(df * 1e8) / 1e8;

    // Generate frequencies
    size_t f_steps = (int)(log(f_max / f_min) / log(1 + df) + 1);
    freq.resize(f_steps);
    freq[0] = f_min;
    for (size_t i = 1; i < f_steps; i++) freq[i] = freq[i - 1] * (1 + df);

    initialize(t_bins, N_bins_min);
}

// Constructor with a fixed array of search frequencies
BLSModel_bf::BLSModel_bf(DataContainer &data_ref,
                         const std::vector<double> &freq,
                         const Target *targetPtr,
                         double t_bins,
                         size_t N_bins_min,
                         DurationMode duration_mode,
                         double min_duration_factor,
                         double max_duration_factor) :
    BLSModel(
        data_ref, 0, 0, targetPtr, duration_mode, nullptr, min_duration_factor, max_duration_factor)
{
    // Set min and max frequencies
    f_min = *std::min_element(freq.begin(), freq.end());
    f_max = *std::max_element(freq.begin(), freq.end());

    dt_per_step = 0;
    this->freq = freq;

    initialize(t_bins, N_bins_min);
}

// Calculate the number of tested frequencies without running the model
size_t BLSModel_bf::calculate_N_freq()
{
    return freq.size();
}

// Make a duplicate (newly initialized object with the same settings)
std::unique_ptr<BLSModel> BLSModel_bf::duplicate() const
{
    return std::unique_ptr<BLSModel>(new BLSModel_bf(*data,
                                                     freq,
                                                     target,
                                                     t_bins,
                                                     N_bins_min,
                                                     duration_mode,
                                                     min_duration_factor,
                                                     max_duration_factor));
}

// Initial operations (called by constructors)
void BLSModel_bf::initialize(double t_bins, size_t N_bins_min)
{
    if (t_bins > 0)
        this->t_bins = t_bins;
    if (N_bins_min > 0)
        this->N_bins_min = N_bins_min;

    // Expand chi2 vectors
    chi2.resize(N_freq());
    dchi2.resize(N_freq());
    chi2r.resize(N_freq());
    chi2_mag0.resize(N_freq());
    chi2_dmag.resize(N_freq());
    //chi2_dmag_err.assign(N_freq(), 0);
    chi2_t0.resize(N_freq());
    chi2_dt.resize(N_freq());
    N_bins.resize(N_freq());
}

void BLSModel_bf::run(bool verbose, bool full_results)
{
    double P, Z, Zi, mi, dchi2_, dchi2_max, m0_best, dmag_best;
    double dt_min_P, dt_max_P;
    size_t N_bins_, N_bins_real, t_start_best, dt_best, dt_min, dt_max;

    // Arrays for binned magnitudes
    size_t N_bins_max = std::max(N_bins_min, (size_t)(1 / f_min / t_bins));
    double mag[N_bins_max], mag_err[N_bins_max];

    // Process data
    data->calculate_mag_frac();

    for (size_t i = 0; i < N_freq(); i++) {
        if ((verbose) and (i % std::max(1, (int)(N_freq() / 100)) == 0))
            std::cout << "BLS     NFREQ: " << N_freq()
                      << "     STATUS: " << (int)(100 * i / N_freq()) << "%          \r"
                      << std::flush;

        P = 1. / freq[i];

        dchi2_max = 0;
        m0_best = 0;
        dmag_best = 0;
        t_start_best = 0;
        dt_best = 0;

        // Calculate binned magnitudes
        N_bins_ = std::max(N_bins_min, (size_t)(P / t_bins));
        bin(P, N_bins_, data, mag, mag_err, &N_bins_real);

        // Estimate the range of transit durations
        std::tie(dt_min_P, dt_max_P) = get_duration_limits(P);
        dt_min = std::max((size_t)(1), (size_t)(N_bins_ * dt_min_P / P));
        dt_max = std::max(dt_min, (size_t)(N_bins_ * dt_max_P / P));

        // Obtain the sum of 1 / mag_err^2
        Z = 0;
        for (size_t j = 0; j < N_bins_; j++) Z += 1 / SQ(mag_err[j]);

        // Loop over transit starts (in bins)
        for (size_t t_start = 0; t_start < N_bins_; t_start++) {
            mi = 0;
            Zi = 0;

            // Loop over transit durations (in bins)
            for (size_t dt = 0; dt < dt_max; dt++) {
                mi += mag[(t_start + dt) % N_bins_] / SQ(mag_err[(t_start + dt) % N_bins_]);
                Zi += 1 / SQ(mag_err[(t_start + dt) % N_bins_]);
                dchi2_ = (Zi == Z ? 0 : SQ(mi) / Zi / (1 - Zi / Z));

                if ((dt >= dt_min - 1) && (dchi2_ > dchi2_max)) {
                    dchi2_max = dchi2_;
                    if (full_results) {
                        m0_best = -mi / Z / (1 - Zi / Z);
                        dmag_best = -mi / Zi / (1 - Zi / Z);
                        t_start_best = t_start;
                        dt_best = dt;
                    }
                }
            }
        }

        // Calculate chi2 for the best combination (= chi2_const - dchi2)
        dchi2[i] = dchi2_max;
        if (full_results) {
            chi2[i] = -dchi2_max;
            for (size_t j = 0; j < N_bins_; j++) chi2[i] += SQ(mag[j] / mag_err[j]);
            chi2r[i] = chi2[i] / (N_bins_real - 1);
            chi2_mag0[i] = data->mag_avg * (m0_best + 1);
            chi2_dmag[i] = dmag_best * chi2_mag0[i];
            chi2_t0[i] = P * t_start_best / N_bins_;
            chi2_dt[i] = P * (dt_best + 1) / N_bins_;
            N_bins[i] = N_bins_;
        }
    }

    if (verbose)
        std::cout << "BLS     NFREQ: " << N_freq() << "     STATUS: 100%\n";
}

// Constructor
// If any numeric value is 0 then the default value is used
BLSModel_FFA::BLSModel_FFA(DataContainer &data_ref,
                           double f_min,
                           double f_max,
                           const Target *targetPtr,
                           DurationMode duration_mode,
                           const std::vector<double> *durations,
                           double min_duration_factor,
                           double max_duration_factor,
                           double t_samp,
                           bool downsample,
                           double ds_invpower,
                           double ds_threshold,
                           size_t N_bins_transit_min) :
    BLSModel(data_ref,
             f_min,
             f_max,
             targetPtr,
             duration_mode,
             durations,
             min_duration_factor,
             max_duration_factor)
{
    // Override numeric values if given
    if (t_samp > 0)
        this->t_samp = t_samp;
    if (ds_invpower > 0)
        this->ds_invpower = ds_invpower;
    if (ds_threshold > 0)
        this->ds_threshold = ds_threshold;
    if (N_bins_transit_min > 0)
        this->N_bins_transit_min = N_bins_transit_min;

    this->downsample = downsample;
}

// Calculate the number of tested frequencies without running the model
size_t BLSModel_FFA::calculate_N_freq()
{
    size_t data_size = resample_uniform_size(*data, t_samp);
    return periodogram_length(
        data_size, t_samp, 1. / f_max, 1. / f_min, downsample, ds_invpower, ds_threshold);
}

// Make a duplicate (newly initialized object with the same settings)
std::unique_ptr<BLSModel> BLSModel_FFA::duplicate() const
{
    return std::unique_ptr<BLSModel>(new BLSModel_FFA(*data,
                                                      f_min,
                                                      f_max,
                                                      target,
                                                      duration_mode,
                                                      &durations,
                                                      min_duration_factor,
                                                      max_duration_factor,
                                                      t_samp,
                                                      downsample,
                                                      ds_invpower,
                                                      ds_threshold,
                                                      N_bins_transit_min));
}

// Generate required results
template <typename T>
void BLSModel_FFA::process_results(std::vector<BLSResult<T>> &results, bool full_results)
{
    const size_t N_freq = results.size();
    BLSResult<T> *pres = results.data();

    freq.resize(N_freq);
    dchi2.assign(N_freq, 0);
    if (full_results) {
        chi2_mag0.assign(N_freq, 0);
        chi2_dmag.assign(N_freq, 0);
        //chi2_dmag_err.assign(N_freq, 0);
        chi2_t0.assign(N_freq, 0);
        chi2_dt.assign(N_freq, 0);
        N_bins.assign(N_freq, 0);
        time_spent.assign(N_freq, 0);
    }

    for (size_t i = 0; i < N_freq; i++) {
        freq[i] = 1 / pres->P;
        dchi2[i] = pres->dchi2;
        if (full_results) {
            chi2_mag0[i] = pres->mag0;
            chi2_dmag[i] = pres->dmag;
            //chi2_dmag_err[i] = pres->dmag_err;
            chi2_t0[i] = fmod(rdata->rjd[0] + t_samp * (pres->t0 - 0.5), pres->P);
            chi2_dt[i] = t_samp * pres->dur;
            N_bins[i] = pres->N_bins;
            time_spent[i] = pres->time_spent;
        }
        pres++;
    }
}

void BLSModel_FFA::run(bool verbose, bool full_results)
{
    run_prec<float>(verbose, full_results);
}

void BLSModel_FFA::run_double(bool verbose, bool full_results)
{
    run_prec<double>(verbose, full_results);
}

// Data will be resampled uniformly to cadence tsamp
template <typename T> void BLSModel_FFA::run_prec(bool verbose, bool full_results)
{
    if (verbose)
        std::cout << "Starting FFA...\n";

    // Resample to desired tsamp
    rdata = resample_uniform(*data, t_samp);
    std::vector<T> mag(rdata->size, 0); // Magnitudes
    std::vector<T> wts(rdata->size, 0); // Weights (1/err^2)
    for (size_t i = 0; i < rdata->size; i++) {
        if (rdata->valid_mask[i]) {
            mag[i] = rdata->mag[i];
            wts[i] = 1. / rdata->err[i] / rdata->err[i];
        }
    }

    // Function wrapper to return the maximum tested transit duration at each period
    // auto get_duration_limits_ =
    //    std::bind(&BLSModel::get_duration_limits, this, std::placeholders::_1);

    auto t_start = std::chrono::high_resolution_clock::now();
    std::vector<BLSResult<T>> pgram =
        std::move(periodogram<T>(mag.data(), wts.data(), mag.size(), *this, verbose));
    auto t_end = std::chrono::high_resolution_clock::now();

    if (verbose) {
        std::chrono::duration<double> rtime = t_end - t_start;
        std::cout << "Number of tested periods: " << pgram.size() << "\n";
        std::cout << "BLS runtime: " << rtime.count() << " sec\n";
    }
    process_results(pgram, full_results);
}

// Generate the noise spectrum
// Currently only generates the necessary BLS spectra but does not analyze them to generate the
// noise spectrum (this would need additional C++ libraries)
std::vector<double> NoiseBLS::generate(size_t N_sim, NoiseMode selection_mode, bool verbose)
{

    // Make sure a model has been set up
    if (!model)
        throw std::runtime_error("Cannot generate a noise spectrum due to an uninitialized model.");

    // Override default selection_mode if given
    if (selection_mode != NoiseMode::None)
        this->selection_mode = selection_mode;

    // Set up information and arrays
    this->N_sim = N_sim;
    size_t N_freq = model->calculate_N_freq();
    Vector2D<double> dchi2_arr(N_sim, N_freq);
    dchi2.resize(N_freq);

    // Hold pointer to original data and allocate a temporary simulated data set
    DataContainer *data_orig = model->data;
    std::unique_ptr<DataContainer> data_sim(data_orig->duplicate());
    model->data = data_sim.get();

    // Set up a random number generator
    std::random_device rd;
    std::mt19937 random_gen(rd());

    auto t_start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < N_sim; i++) {

        if (verbose)
            std::cout << "\rNoise BLS iteration: " << (i + 1) << " / " << N_sim << std::flush;

        // Shuffle the data
        std::shuffle(model->data->mag, model->data->mag + model->data->size, random_gen);

        model->run(false, false);

        // Ensure the generated spectrum has the right number of frequencies
        if (model->N_freq() != N_freq)
            throw std::runtime_error("Generated noise spectrum has " +
                                     std::to_string(model->N_freq()) +
                                     " frequencies (Expected: " +
                                     std::to_string(N_freq) +
                                     ").");

        // Save the dchi2 values
        std::copy(model->dchi2.begin(), model->dchi2.end(), dchi2_arr[i].ptr);
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> rtime = t_end - t_start;

    if (verbose) {
        std::cout << " DONE" << std::endl;
        std::cout << "Total noise BLS runtime: " << rtime.count() << " sec" << std::endl;
    }

    // Make model point to the original data again
    model->data = data_orig;

    // Transpose to a [N_freq x N_sim] vector
    Vector2D<double> dchi2_arr_T(dchi2_arr.transpose());

    // Returns a flattened [N_freq x N_sim] vector for Cython processing (see noisebls.pxi)
    return dchi2_arr_T.data;
}

// Constructor
NoiseBLS::NoiseBLS(BLSModel &model)
{
    this->model = std::move(model.duplicate());
}
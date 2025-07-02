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
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>

// Constructor
// If any numeric value is 0 then the default value is used
BLSModel::BLSModel(DataContainer &data_ref,
                   double f_min,
                   double f_max,
                   const Target *targetPtr,
                   int duration_mode,
                   double min_duration_factor,
                   double max_duration_factor)
{
    data = &data_ref;
    target = targetPtr;

    if (f_min > 0)
        this->f_min = f_min;
    if (f_max > 0)
        this->f_max = f_max;
    if (duration_mode > 0)
        this->duration_mode = duration_mode;
    if (min_duration_factor > 0)
        this->min_duration_factor = min_duration_factor;
    if (max_duration_factor > 0)
        this->max_duration_factor = max_duration_factor;
}

// Get the minimum and maximum tested transit duration at a given period P
std::tuple<double, double> BLSModel::get_duration_limits(double P)
{
    double min_duration, max_duration;

    switch (duration_mode) {
    // Constant duration limits
    case 1:
        min_duration = min_duration_factor;
        max_duration = max_duration_factor;
        break;

    // Duration limits proportional to the orbital period
    case 2:
        min_duration = min_duration_factor * P;
        max_duration = max_duration_factor * P;
        break;

    // Duration limits proportional to the predicted physical transit duration
    case 3:
        if (target == nullptr) {
            throw std::runtime_error("Target must not be null with max_duration_mode == 3.");
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
        throw std::runtime_error("BLSModel::get_max_duration() called with invalid "
                                 "duration_mode = " +
                                 std::to_string(duration_mode));
        return std::make_tuple(0, 0);
    }

    return std::make_tuple(min_duration, max_duration);
}

// Get number of frequencies
size_t BLSModel::N_freq()
{
    return freq.size();
}

void BLSModel::run(bool verbose)
{
    std::cout << "run() is not defined for an object of type " << typeid(*this).name();
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
                         int duration_mode,
                         double min_duration_factor,
                         double max_duration_factor) :
    BLSModel(
        data_ref, f_min, f_max, targetPtr, duration_mode, min_duration_factor, max_duration_factor)
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
                         int duration_mode,
                         double min_duration_factor,
                         double max_duration_factor) :
    BLSModel(data_ref, 0, 0, targetPtr, duration_mode, min_duration_factor, max_duration_factor)
{
    // Set min and max frequencies
    f_min = *std::min_element(freq.begin(), freq.end());
    f_max = *std::max_element(freq.begin(), freq.end());

    dt_per_step = 0;
    this->freq = freq;

    initialize(t_bins, N_bins_min);
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
    chi2_t0.resize(N_freq());
    chi2_dt.resize(N_freq());
    N_bins.resize(N_freq());
}

void BLSModel_bf::run(bool verbose)
{
    double P, Z, Zi, mi, dchi2_, dchi2_min, m0_best, dmag_best;
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

        P = 1 / freq[i];

        dchi2_min = 0;
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
                dchi2_ = -(Zi == Z ? 0 : SQ(mi) / Zi / (1 - Zi / Z));

                if ((dt >= dt_min - 1) and (dchi2_ < dchi2_min)) {
                    m0_best = -mi / Z / (1 - Zi / Z);
                    dmag_best = -mi / Zi / (1 - Zi / Z);
                    t_start_best = t_start;
                    dt_best = dt;
                    dchi2_min = dchi2_;
                }
            }
        }

        // Calculate chi2 for the best combination (= dchi2 + chi2_const)
        dchi2[i] = dchi2_min;
        chi2[i] = dchi2_min;
        for (size_t j = 0; j < N_bins_; j++) chi2[i] += SQ(mag[j] / mag_err[j]);
        chi2r[i] = chi2[i] / (N_bins_real - 1);
        chi2_mag0[i] = data->mag_avg * (m0_best + 1);
        chi2_dmag[i] = dmag_best * chi2_mag0[i];
        chi2_t0[i] = P * t_start_best / N_bins_;
        chi2_dt[i] = P * (dt_best + 1) / N_bins_;
        N_bins[i] = N_bins_;
    }

    if (verbose)
        std::cout << "BLS     NFREQ: " << N_freq() << "     STATUS: 100%\n";
}

// Generate required results
template <typename T> void BLSModel_FFA::process_results(std::vector<BLSResult<T>> &results)
{
    const size_t N_freq = results.size();
    BLSResult<T> *pres = results.data();

    freq.resize(N_freq);
    dchi2.assign(N_freq, 0);
    chi2_mag0.assign(N_freq, 0);
    chi2_dmag.assign(N_freq, 0);
    chi2_t0.assign(N_freq, 0);
    chi2_dt.assign(N_freq, 0);
    N_bins.assign(N_freq, 0);

    for (size_t i = 0; i < N_freq; i++) {
        freq[i] = 1 / pres->P;
        dchi2[i] = -(pres->dchi2);
        chi2_mag0[i] = pres->mag0;
        chi2_dmag[i] = pres->dmag;
        chi2_t0[i] = fmod(rdata->rjd[0] + pres->P * (pres->t0 - 0.5) / pres->N_bins, pres->P);
        chi2_dt[i] = t_samp * pres->dur;
        N_bins[i] = pres->N_bins;
        pres++;
    }
}

void BLSModel_FFA::run(bool verbose)
{
    run_prec<float>(verbose);
}

void BLSModel_FFA::run_double(bool verbose)
{
    run_prec<double>(verbose);
}

// Data will be resampled uniformly to cadence tsamp
template <typename T> void BLSModel_FFA::run_prec(bool verbose)
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
            wts[i] = 1 / rdata->err[i] / rdata->err[i];
        }
    }

    /*size_t length = periodogram_length(mag.size(), t_samp, 1/f_max, 1/f_min);
    if (verbose)
        std::cout << "Number of tested periods: " << length << "\n";
    freq.assign(length, 0);
    foldbins.assign(length, 0);
    dchi2.assign(length, 0);
    chi2_t0.assign(length, 0);*/

    // Function wrapper to return the maximum tested transit duration at each period
    auto get_duration_limits_ =
        std::bind(&BLSModel::get_duration_limits, this, std::placeholders::_1);

    auto t_start = std::chrono::high_resolution_clock::now();
    std::vector<BLSResult<T>> pgram = std::move(periodogram<T>(
        mag.data(), wts.data(), mag.size(), t_samp, get_duration_limits_, 1 / f_max, 1 / f_min));
    auto t_end = std::chrono::high_resolution_clock::now();

    if (verbose) {
        std::chrono::duration<double> rtime = t_end - t_start;
        std::cout << "Number of tested periods: " << pgram.size() << "\n";
        std::cout << "BLS runtime: " << rtime.count() << " sec\n";
    }
    process_results(pgram);
}
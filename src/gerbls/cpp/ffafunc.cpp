/*
    This code uses portions of the riptide package, with modifications - see riptide/LICENSE.
    The original riptide code is included in riptide/periodogram.hpp.
    Portions under Copyright (c) 2017-2021 Vincent Morello
    Portions under Copyright (c) 2025 Kristo Ment
*/

#include "ffafunc.hpp"
#include "model.hpp"
#include "physfunc.hpp"
#include "riptide/block.hpp"
#include "riptide/kernels.hpp"
#include "riptide/periodogram.hpp"
#include "riptide/snr.hpp"
#include "riptide/transforms.hpp"
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric> // std::accumulate
#include <tuple>
#include <vector>

// out = x - y
template <typename T>
void array_diff(const T *__restrict__ x,
                const T *__restrict__ y,
                const size_t size,
                T *__restrict__ out)
{
    for (size_t i = 0; i < size; i++) out[i] = x[i] - y[i];
}

// Find maximum dchi2 given sum arrays of in-transit (mag*wts) and (wts)
// Returns: BLSResult corresponding to dchi2_max
// wtotal = sum(wts)
template <typename T>
void array_dchi2_max(const T *__restrict__ prod,
                     const T *__restrict__ wts,
                     const size_t size,
                     const T wtotal,
                     BLSResult<T> &result,
                     const size_t width)
{
    size_t i_max = 0;
    T dchi2, dchi2_max = 0;

    for (size_t i = 0; i < size; i++) {
        dchi2 = prod[i] * prod[i] / wts[i] / (1 - wts[i] / wtotal);
        if (dchi2 > dchi2_max) {
            i_max = i;
            dchi2_max = dchi2;
        }
    }

    if (dchi2_max > result.dchi2) {
        result.t0 = ((i_max < size - 1) ? i_max + 1 : 0);
        result.dur = width;
        result.mag0 = -prod[i_max] / wtotal / (1 - wts[i_max] / wtotal);
        result.dmag = -prod[i_max] / wts[i_max] / (1 - wts[i_max] / wtotal);
        result.dchi2 = dchi2_max;
    }
}

// Output: results (must be array of size wprod.rows)
template <typename T>
void chisq_2d(riptide::ConstBlock<T> wprod,
              riptide::ConstBlock<T> weights,
              const std::vector<size_t> &widths,
              BLSResult<T> *results)
{
    for (size_t i = 0; i < wprod.rows; ++i) {
        chisq_row<T>(wprod.rowptr(i), weights.rowptr(i), wprod.cols, widths, *results);
        results++;
    }
}

// Output: results (must be array of size wprod.rows)
template <typename T>
void chisq_2d(riptide::ConstBlock<T> wprod,
              riptide::ConstBlock<T> weights,
              const size_t min_width,
              const size_t max_width,
              BLSResult<T> *results)
{
    for (size_t i = 0; i < wprod.rows; ++i) {
        chisq_row<T>(
            wprod.rowptr(i), weights.rowptr(i), wprod.cols, min_width, max_width, *results);
        results++;
    }
}

// Evaluate dchi2 for a fixed period (FFA row)
// Output: BLSResult corresponding to the highest dchi2
template <typename T>
void chisq_row(const T *__restrict__ wprod,
               const T *__restrict__ wts,
               const size_t size,
               const std::vector<size_t> &widths,
               BLSResult<T> &result)
{
    size_t previous_width = 0;
    const size_t max_width = *std::max_element(widths.begin(), widths.end());
    T cpfsum1[size + max_width];
    T cpfsum2[size + max_width];
    T inmag[size];
    T inwts[size];
    riptide::circular_prefix_sum<T>(wprod, size, size + max_width, cpfsum1);
    riptide::circular_prefix_sum<T>(wts, size, size + max_width, cpfsum2);
    const T wtotal = cpfsum2[size - 1]; // sum of weights

    for (size_t width : widths) {
        if (width < 1 || width == previous_width)
            continue;
        array_diff<T>(cpfsum1 + width, cpfsum1, size, inmag);
        array_diff<T>(cpfsum2 + width, cpfsum2, size, inwts);
        array_dchi2_max<T>(inmag, inwts, size, wtotal, result, width);
        previous_width = width;
    }
}

// Evaluate dchi2 for a fixed period (FFA row)
// Output: BLSResult corresponding to the highest dchi2
template <typename T>
void chisq_row(const T *__restrict__ wprod,
               const T *__restrict__ wts,
               const size_t size,
               const size_t min_width,
               const size_t max_width,
               BLSResult<T> &result)
{
    T cpfsum1[size + max_width];
    T cpfsum2[size + max_width];
    T inmag[size];
    T inwts[size];
    riptide::circular_prefix_sum<T>(wprod, size, size + max_width, cpfsum1);
    riptide::circular_prefix_sum<T>(wts, size, size + max_width, cpfsum2);
    const T wtotal = cpfsum2[size - 1]; // sum of weights

    for (size_t width = min_width; width <= max_width; width++) {
        array_diff<T>(cpfsum1 + width, cpfsum1, size, inmag);
        array_diff<T>(cpfsum2 + width, cpfsum2, size, inwts);
        array_dchi2_max<T>(inmag, inwts, size, wtotal, result, width);
    }
}

// Compute the periodogram of a time series that has been normalised to zero mean and unit variance.
// P Output: vector of BLSResult<T> containing the best fit for each tested orbital period
template <typename T>
std::vector<BLSResult<T>> periodogram(const T *__restrict__ mag,
                                      const T *__restrict__ wts,
                                      size_t size,
                                      const BLSModel_FFA &model,
                                      bool verbose)
{
    // periodogram_check_arguments(size, tsamp, period_min, period_max, bins_min, bins_max);

    // Model parameters
    const double period_min = 1. / model.f_max;
    const double period_max = 1. / model.f_min;

    // Temporary variables
    double min_width_P, max_width_P;
    std::vector<size_t> widths(model.durations.size()); // Transit duration widths to search
    size_t bstart, bstop;
    T *wprod_, *weights_;

    // Maximum possible size of the data
    const size_t n_max = size + (size_t)(period_max / model.t_samp);

    // Calculate periodogram length and allocate memory for output
    const size_t length = periodogram_length(size,
                                             model.t_samp,
                                             period_min,
                                             period_max,
                                             model.downsample,
                                             model.ds_invpower,
                                             model.ds_threshold);
    std::vector<BLSResult<T>> results(length);
    BLSResult<T> *presult = results.data();

    // Geometric growth factor for the downsampling parameter
    double ds_geo = model.ds_threshold;

    // Minimum and maximum number of bins in each downsampling loop iteration
    size_t bins_min = period_min / model.t_samp;
    // size_t bins_max = downsample ? bins_min * ds_threshold + 1 : period_max / tsamp;

    // Number of required downsampling cycles
    const size_t num_downsamplings =
        model.downsample ? ceil(log(period_max / period_min) / log(ds_geo) / model.ds_invpower) : 1;

    if (verbose) {
        if (model.downsample)
            std::cout << "Downsampling: ON     Number of downsamplings: " << num_downsamplings
                      << "\n";
        else
            std::cout << "Downsampling: OFF\n";
    }

    // Allocate buffers
    const size_t bufsize = n_max;
    std::unique_ptr<T[]> ds_wprod(new T[bufsize]);
    std::unique_ptr<T[]> ds_weights(new T[bufsize]);
    std::unique_ptr<T[]> ffabuf_mem(new T[bufsize]);
    std::unique_ptr<T[]> ffaout_mem1(new T[bufsize]);
    std::unique_ptr<T[]> ffaout_mem2(new T[bufsize]);
    T *ffabuf = ffabuf_mem.get();
    T *ffamag = ffaout_mem1.get();
    T *ffawts = ffaout_mem2.get();

    // Calculate weights and weighted mags
    std::unique_ptr<T[]> weights(new T[n_max]);
    std::unique_ptr<T[]> wprod(new T[n_max]); // Product of mag and weights
    double ftotal = 0;                        // Weighted sum of fluxes
    double wtotal = 0;                        // Sum of weights
    for (size_t i = 0; i < size; i++) {
        ftotal += mag[i] * wts[i];
        wtotal += wts[i];
    }
    const double favg = ftotal / wtotal; // Weighted mean flux
    for (size_t i = 0; i < size; i++) {
        weights[i] = favg * favg * wts[i];
        wprod[i] = (mag[i] / favg - 1) * weights[i];
    }
    // Pad data with zeros to ensure the last FFA row never gets cropped
    for (size_t i = size; i < n_max; i++) {
        wprod[i] = 0;
        weights[i] = 0;
    }

    /* Downsampling loop */
    for (size_t ids = 0; ids < num_downsamplings; ++ids) {

        const double f = pow(ds_geo, ids);  // current downsampling factor
        const double tau = f * model.t_samp; // current sampling time
        const double period_max_samples = period_max / tau;
        const size_t n = riptide::downsampled_size(size, f); // current number of real input samples
        const size_t n_padded = n + period_max_samples;      // number of samples after zero-padding

        // Initial downsampling factor is 1
        if (ids == 0) {
            wprod_ = wprod.get();
            weights_ = weights.get();
        }
        else {
            wprod_ = ds_wprod.get();
            weights_ = ds_weights.get();

            // Downsample the input data
            riptide::downsample<T>(wprod.get(), size, f, wprod_);
            riptide::downsample<T>(weights.get(), size, f, weights_);

            // Pad data with zeros to ensure the last FFA row never gets cropped
            for (size_t i = n; i < n_padded; i++) {
                wprod_[i] = 0;
                weights_[i] = 0;
            }
        }

        // Min and max number of bins with which to FFA transform
        // NOTE: bstop is INclusive
        if (model.downsample) {
            bstart = bins_min * pow(model.ds_threshold, (model.ds_invpower - 1) * ids);
            bstop =
                std::min(bstart * pow(model.ds_threshold, model.ds_invpower), period_max_samples);
        }
        else {
            bstart = bins_min;
            bstop = period_max_samples;
        }

        // Pre-calculate tested transit widths if they are constant and explicitly defined
        if (model.explicit_durations() && model.duration_mode == 1) {
            model.set_widths(0., tau, widths);
        }

        /* FFA transform loop */
        for (size_t bins = bstart; bins <= bstop; ++bins) {

            auto t_start = std::chrono::high_resolution_clock::now();

            const size_t rows = n_padded / bins;
            const double period_ceil = std::min(period_max_samples, bins + 1.);
            const size_t rows_eval = std::min(rows, riptide::ceilshift(rows, bins, period_ceil));

            riptide::transform<T>(wprod_, rows, bins, ffabuf, ffamag);
            riptide::transform<T>(weights_, rows, bins, ffabuf, ffawts);

            auto block1 = riptide::ConstBlock<T>(ffamag, rows_eval, bins);
            auto block2 = riptide::ConstBlock<T>(ffawts, rows_eval, bins);

            if (model.explicit_durations()) {
                if (model.duration_mode != 1) {
                    model.set_widths((bins + 1) * tau, tau, widths);
                }
                chisq_2d<T>(block1, block2, widths, presult);
            }
            else {
                std::tie(min_width_P, max_width_P) = model.get_duration_limits((bins + 1) * tau);
                const size_t min_width = std::max((size_t)(1), (size_t)(min_width_P / tau));
                const size_t max_width = std::max(min_width, (size_t)(max_width_P / tau));
                chisq_2d<T>(block1, block2, min_width, max_width, presult);
            }

            auto t_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> rtime = t_end - t_start;

            for (size_t s = 0; s < rows_eval; ++s) {
                presult[s].P = tau * bins * bins / (bins - s / (rows - 1.0));
                presult[s].mag0 = favg * (presult[s].mag0 + 1);
                presult[s].dmag *= presult[s].mag0;
                presult[s].N_bins = bins;
                presult[s].time_spent = rtime.count() / rows_eval;
            }

            presult += rows_eval;
        }
    }

    return results;
}

// Explicit instantiations for float and double
template std::vector<BLSResult<float>> periodogram(
    const float *__restrict__, const float *__restrict__, size_t, const BLSModel_FFA &, bool);
template std::vector<BLSResult<double>> periodogram(
    const double *__restrict__, const double *__restrict__, size_t, const BLSModel_FFA &, bool);

/*
Returns the total number of trial periods in a periodogram
*/
size_t periodogram_length(size_t size,
                          double tsamp,
                          double period_min,
                          double period_max,
                          bool downsample,
                          double ds_invpower,
                          double ds_threshold)
{
    // periodogram_check_arguments(size, tsamp, period_min, period_max, bins_min, bins_max);

    size_t bstart, bstop;

    // Geometric growth factor for the downsampling factor
    double ds_geo = ds_threshold;

    // Minimum and maximum number of bins in each downsampling loop iteration
    size_t bins_min = period_min / tsamp;
    // size_t bins_max = downsample ? bins_min * ds_threshold + 1 : period_max / tsamp;

    // Number of required downsampling cycles
    const size_t num_downsamplings =
        downsample ? ceil(log(period_max / period_min) / log(ds_geo) / ds_invpower) : 1;
    size_t length = 0; // total number of period trials, to be calculated

    /* Downsampling loop */
    for (size_t ids = 0; ids < num_downsamplings; ++ids) {

        const double f = pow(ds_geo, ids); // current downsampling factor
        const double tau = f * tsamp;      // current sampling time
        const double period_max_samples = period_max / tau;
        const size_t n = riptide::downsampled_size(size, f); // current number of real input samples
        const size_t n_padded = n + period_max_samples;      // number of samples after zero-padding

        if (downsample) {
            bstart = bins_min * pow(ds_threshold, (ds_invpower - 1) * ids);
            bstop = std::min(bstart * pow(ds_threshold, ds_invpower), period_max_samples);
        }
        else {
            bstart = bins_min;
            bstop = period_max_samples;
        }

        /* FFA transform loop */
        for (size_t bins = bstart; bins <= bstop; ++bins) {
            const size_t rows = n_padded / bins;
            const double period_ceil = std::min(period_max_samples, bins + 1.);
            const size_t rows_eval = std::min(rows, riptide::ceilshift(rows, bins, period_ceil));
            length += rows_eval;
        }
    }

    return length;
}

// Resample the light curve with a uniform sampling interval
// Assumes the data are already time-sorted
std::unique_ptr<DataContainer> resample_uniform(const DataContainer &data, double tsamp)
{
    std::unique_ptr<DataContainer> out(new DataContainer);
    size_t N_sampled = ceil((data.rjd[data.size - 1] - data.rjd[0]) / tsamp + 0.5);
    out->allocate(N_sampled);
    out->valid_mask.reset(new bool[N_sampled]);
    size_t i = 0;

    for (size_t j = 0; j < N_sampled; j++) {
        out->rjd[j] = data.rjd[0] + j * tsamp;
        out->mag[j] = 0;
        out->err[j] = 0;
        while ((i < data.size) && (data.rjd[i] < data.rjd[0] + (j + 0.5) * tsamp)) {
            out->mag[j] += data.mag[i] / data.err[i] / data.err[i];
            out->err[j] += 1 / data.err[i] / data.err[i];
            i++;
        }
    }

    for (size_t j = 0; j < N_sampled; j++) {
        if (out->err[j] > 0) {
            out->mag[j] /= out->err[j];
            out->err[j] = sqrt(1 / out->err[j]);
            out->valid_mask[j] = true;
        }
        else {
            out->err[j] = 1e10;
            out->valid_mask[j] = false;
        }
    }

    return out;
}

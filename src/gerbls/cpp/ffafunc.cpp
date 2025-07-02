/*
    This code uses portions of the riptide package, with modifications - see riptide/LICENSE.
    The original riptide code is included in riptide/periodogram.hpp.
    Portions under Copyright (c) 2017-2021 Vincent Morello
    Portions under Copyright (c) 2025 Kristo Ment
*/

#include "ffafunc.hpp"
#include "physfunc.hpp"
#include "riptide/block.hpp"
#include "riptide/kernels.hpp"
#include "riptide/periodogram.hpp"
#include "riptide/snr.hpp"
#include "riptide/transforms.hpp"
#include <cmath>
#include <functional>
#include <iostream>
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
              const size_t min_width,
              const size_t max_width,
              BLSResult<T> *results)
{
    for (size_t i = 0; i < wprod.rows; ++i) {
        // snr1(block.rowptr(i), block.cols, widths, stdnoise, out);
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
// get_duration_limits(P) should return the min and max transit duration at the given orbital period
// P Output: vector of BLSResult<T> containing the best fit for each tested orbital period
template <typename T>
std::vector<BLSResult<T>>
    periodogram(const T *__restrict__ mag,
                const T *__restrict__ wts,
                size_t size,
                double tsamp,
                // const std::vector<size_t>& widths,
                std::function<std::tuple<double, double>(double)> get_duration_limits,
                double period_min,
                double period_max)
{
    // periodogram_check_arguments(size, tsamp, period_min, period_max, bins_min, bins_max);

    // Temporary variables
    double min_width_P, max_width_P;

    // Calculate periodogram length and allocate memory for output
    const size_t length = periodogram_length(size, tsamp, period_min, period_max);
    std::vector<BLSResult<T>> results(length);
    BLSResult<T> *presult = results.data();

    // Pad data with zeros to ensure the last FFA row never gets cropped
    const double period_max_samples = period_max / tsamp;
    const size_t n = size + period_max_samples; // Total size of the data

    // Allocate buffers
    const size_t bufsize = n;
    std::unique_ptr<T[]> ffabuf_mem(new T[bufsize]);
    std::unique_ptr<T[]> ffaout_mem1(new T[bufsize]);
    std::unique_ptr<T[]> ffaout_mem2(new T[bufsize]);
    T *ffabuf = ffabuf_mem.get();
    T *ffamag = ffaout_mem1.get();
    T *ffawts = ffaout_mem2.get();

    /* Downsampling loop */
    // const double tau = tsamp; // current sampling time
    // const size_t n = size; // current number of input samples
    // const size_t num_widths = widths.size();

    // Calculate weights and weighted mags
    std::unique_ptr<T[]> weights(new T[n]);
    std::unique_ptr<T[]> wprod(new T[n]); // Product of wmag and weights
    double ftotal = 0;                    // Weighted sum of fluxes
    double wtotal = 0;                    // Sum of weights
    for (size_t i = 0; i < size; i++) {
        ftotal += mag[i] * wts[i];
        wtotal += wts[i];
    }
    const double favg = ftotal / wtotal; // Weighted mean flux
    for (size_t i = 0; i < size; i++) {
        weights[i] = favg * favg * wts[i];
        wprod[i] = (mag[i] / favg - 1) * weights[i];
    }
    for (size_t i = size; i < n; i++) {
        weights[i] = 0;
        wprod[i] = 0;
    }

    // Min and max number of bins with which to FFA transform in order to
    // cover all trial periods between period_min and period_max.
    // NOTE: bstop is INclusive
    // Also, we MUST enforce bstop <= n, to avoid doing an FFA transform with 0 rows
    const size_t bstart = period_min / tsamp;
    const size_t bstop = std::min({n, size_t(period_max / tsamp)});

    /* FFA transform loop */
    for (size_t bins = bstart; bins <= bstop; ++bins) {
        const size_t rows = n / bins;
        // const float stdnoise = 1.0; //sqrt(rows * downsampled_variance(size, f));
        const double period_ceil = std::min(period_max_samples, bins + 1.0);
        const size_t rows_eval = std::min(rows, riptide::ceilshift(rows, bins, period_ceil));
        std::tie(min_width_P, max_width_P) = get_duration_limits((bins + 1) * tsamp);
        const size_t min_width = std::max((size_t)(1), (size_t)(min_width_P / tsamp));
        const size_t max_width = std::max(min_width, (size_t)(max_width_P / tsamp));

        riptide::transform<T>(wprod.get(), rows, bins, ffabuf, ffamag);
        riptide::transform<T>(weights.get(), rows, bins, ffabuf, ffawts);

        auto block1 = riptide::ConstBlock<T>(ffamag, rows_eval, bins);
        auto block2 = riptide::ConstBlock<T>(ffawts, rows_eval, bins);
        chisq_2d<T>(block1, block2, min_width, max_width, presult);

        for (size_t s = 0; s < rows_eval; ++s) {
            presult[s].P = tsamp * bins * bins / (bins - s / (rows - 1.0));
            presult[s].mag0 = favg * (presult[s].mag0 + 1);
            presult[s].dmag *= presult[s].mag0;
            presult[s].N_bins = bins;
        }

        presult += rows_eval;
    }

    return results;
}

// Explicit instantiations for float and double
template std::vector<BLSResult<float>>
    periodogram(const float *__restrict__,
                const float *__restrict__,
                size_t,
                double,
                std::function<std::tuple<double, double>(double)>,
                double,
                double);
template std::vector<BLSResult<double>>
    periodogram(const double *__restrict__,
                const double *__restrict__,
                size_t,
                double,
                std::function<std::tuple<double, double>(double)>,
                double,
                double);

/*
Returns the total number of trial periods in a periodogram
*/
size_t periodogram_length(size_t size, double tsamp, double period_min, double period_max)
// size_t bins_min,
// size_t bins_max)
{
    // periodogram_check_arguments(size, tsamp, period_min, period_max, bins_min, bins_max);

    size_t length = 0; // total number of period trials, to be calculated

    /* Downsampling loop */
    // const double tau = f * tsamp; // current sampling time
    const double period_max_samples = period_max / tsamp;
    // const size_t n = size; // current number of input samples

    // Pad data with zeros to ensure the last FFA row never gets cropped
    const size_t n = size + period_max_samples; // Total size of the data

    const size_t bstart = period_min / tsamp;
    const size_t bstop = period_max / tsamp;

    /* FFA transform loop */
    for (size_t bins = bstart; bins <= bstop; ++bins) {
        const size_t rows = n / bins;
        const double period_ceil = std::min(period_max_samples, bins + 1.0);
        const size_t rows_eval = std::min(rows, riptide::ceilshift(rows, bins, period_ceil));
        length += rows_eval;
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

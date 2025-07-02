#ifndef FFAFUNC_HPP_
#define FFAFUNC_HPP_

#include "structure.hpp"
#include <functional>
#include <tuple>
#include <vector>

// Forward declarations to avoid compiler errors
namespace riptide {
    template <typename T> class BlockTemplate;
    template <typename T> using ConstBlock = BlockTemplate<const T>;
}

template <typename T> struct BLSResult {
    T P;
    size_t t0;
    size_t dur;
    T mag0;
    T dmag;
    T dchi2;
    size_t N_bins;
};

template <typename T>
void array_diff(const T *__restrict__ x,
                const T *__restrict__ y,
                const size_t size,
                T *__restrict__ out);
template <typename T>
void array_dchi2_max(const T *__restrict__ prod,
                     const T *__restrict__ wts,
                     const size_t size,
                     const T wtotal,
                     BLSResult<T> &result,
                     const size_t width);
template <typename T>
void chisq_2d(riptide::ConstBlock<T> wprod,
              riptide::ConstBlock<T> weights,
              const size_t min_width,
              const size_t max_width,
              BLSResult<T> *results);
template <typename T>
void chisq_row(const T *__restrict__ wprod,
               const T *__restrict__ wts,
               const size_t size,
               const size_t min_width,
               const size_t max_width,
               BLSResult<T> &result);

template <typename T>
std::vector<BLSResult<T>>
    periodogram(const T *__restrict__ mag,
                const T *__restrict__ wts,
                size_t size,
                double tsamp,
                // const std::vector<size_t>& widths,
                std::function<std::tuple<double, double>(double)> get_max_duration,
                double period_min,
                double period_max);

size_t periodogram_length(size_t size, double tsamp, double period_min, double period_max);

// Resample the light curve with a uniform sampling interval
std::unique_ptr<DataContainer> resample_uniform(const DataContainer &data, double tsamp);

#endif /* FFAFUNC_HPP_ */

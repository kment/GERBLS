#ifndef FFAFUNC_HPP_
#define FFAFUNC_HPP_

#include "structure.hpp"
#include <tuple>
#include <vector>

// Forward declarations to avoid compiler errors
namespace riptide {
    template <typename T> class BlockTemplate;
    template <typename T> using ConstBlock = BlockTemplate<const T>;
}

template <typename T>
struct BLSResult {
    T P;
    size_t t0;
    size_t dur;
    T mag0;
    T dmag;
    T dchi2;
    size_t N_bins;
} ;

template <typename T> void array_diff(const T*, const T*, const size_t, T*);
template <typename T> void array_dchi2_max(const T*, const T*, const size_t, const T, BLSResult<T>&, const size_t);
template <typename T> void chisq_2d(riptide::ConstBlock<T>, riptide::ConstBlock<T>, const size_t, BLSResult<T>*);
template <typename T> void chisq_row(const T* __restrict__, const T* __restrict__, const size_t, const size_t, BLSResult<T>&);

template <typename T>
std::vector<BLSResult<T>> periodogram(
    const T* __restrict__ mag,
    const T* __restrict__ wts,
    size_t size,
    double tsamp,
    //const std::vector<size_t>& widths,
    const Target&,
    double period_min,
    double period_max);

size_t periodogram_length(size_t, double, double, double);

// Resample the light curve with a uniform sampling interval
std::unique_ptr<DataContainer> resample_uniform(const DataContainer&, double);

#endif /* FFAFUNC_HPP_ */

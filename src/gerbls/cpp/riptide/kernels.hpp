/*
    This code originates from the riptide package - see LICENSE.
    Copyright (c) 2017-2021 Vincent Morello
    Minor modifications by Kristo Ment (5/16/25).
*/

#ifndef KERNELS_HPP
#define KERNELS_HPP

#include <cstddef> // size_t
#include <cstring> // memcpy
#include <algorithm> // std::min

namespace riptide {

/* Add x with y, and store the result in z */
template <typename T>
void add(const T* __restrict__ x, const T* __restrict__ y, size_t size, T* __restrict__ z)
    {
    for (size_t i = 0; i < size; ++i)
        z[i] = x[i] + y[i];
    }

/* 
Add x with y rolled by shift elements to the LEFT, and store the result in z. shift must be positive.
In numpy that would equivalent to: z = x + roll(y, -shift)
*/
template <typename T>
void fused_rollback_add(const T* __restrict__ x, const T* __restrict__ y, size_t size, size_t shift, T* __restrict__ z)
    {
    const size_t p = shift % size;
    const size_t q = size - p;
    add<T>(x, y + p, q, z);
    add<T>(x + q, y, p, z + q);
    }


/*
Rotate x backwards by shift elements and stored the result in out. shift must be positive.
In numpy that would equivalent to: out = roll(x, -shift)
*/
void rollback(const float* __restrict__ x, size_t size, size_t shift, float* __restrict__ out)
    {
    const size_t p = shift % size;
    const size_t q = size - p;
    memcpy(out, x + p, q * sizeof(float));
    memcpy(out + q, x, p * sizeof(float));
    }


// out = x + a
template <typename T>
void add_scalar(const T* __restrict__ x, size_t size, const T a, T* __restrict__ out)
    {
    for (size_t i = 0; i < size; ++i)
        out[i] = x[i] + a;
    }


// max(x[i] - y[i])
float diff_max(const float* __restrict__ x, const float* __restrict__ y, size_t size)
    {
    float dmax = x[0] - y[0];
    for (size_t i = 1; i < size; ++i)
        {
        const float d = x[i] - y[i];
        if (d > dmax)
            dmax = d;
        }
    return dmax;
    }

/*
size is the size of the input array x, while nsum is the length of the prefix summation.
out will be filled with the prefix sum of x as if its elements were circularly repeating beyond its size.

out[0] = x[0]
out[1] = x[0] + x[1]
...
out[size-1] = x[0] + x[1] + ... + x[size-1]
out[size] = x[0] + x[1] + ... + x[size-1] + x[0]
out[size+1] = x[0] + x[1] + ... + x[size-1] + x[0] + x[1]
...
*/
template <typename T>
void circular_prefix_sum(const T* __restrict__ x, size_t size, size_t nsum, T* __restrict__ out)
    {
    // IMPORTANT: use a double precision accumulator
    double acc = 0;
    const size_t jmax = std::min(size, nsum);

    for (size_t j = 0; j < jmax; ++j)
        {
        acc += x[j];
        out[j] = acc;
        }

    if (nsum <= size)
        return;

    // We can now cast back the input sum to float, which results in faster 
    // code and costs almost no numerical precision
    const T sumx = acc;
    const size_t q = nsum / size; // number of full summation wraps 
    const size_t r = nsum % size; // additional elements to sum beyond last full wrap

    // We have already done a full summation wrap at this point,
    // hence starting at i = 1
    for (size_t i = 1; i < q; ++i)
        add_scalar<T>(out, size, i * sumx, out + i * size);

    // Last incomplete wrap
    add_scalar<T>(out, r, q * sumx, out + q * size);
    }

} // namespace riptide

#endif // KERNELS_HPP
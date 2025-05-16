/*
    This code originates from the riptide package - see LICENSE.
    Copyright (c) 2017-2021 Vincent Morello
    Minor modifications by Kristo Ment (5/16/25).
*/

#ifndef TRANSFORMS_HPP
#define TRANSFORMS_HPP

#include <cstddef> // size_t
#include <cstring> // memcpy()

#include "kernels.hpp"
#include "block.hpp"


namespace riptide {

template <typename T>
void merge(ConstBlock<T> thead, ConstBlock<T> ttail, Block<T> out)
    {
    const size_t m = out.rows;
    const size_t p = out.cols;
    const float kh = (thead.rows - 1.0) / (m - 1.0);
    const float kt = (ttail.rows - 1.0) / (m - 1.0);

    for (size_t s = 0; s < m; ++s)
        {
        const size_t h = kh * s + 0.5;
        const size_t t = kt * s + 0.5;
        const size_t b = s - (h + t);
        fused_rollback_add<T>(thead.rowptr(h), ttail.rowptr(t), p, h+b, out.rowptr(s));
        }
    }


template <typename T>
void transform(ConstBlock<T> input, Block<T> temp, Block<T> out)
    {
    const size_t m = input.rows;
    const size_t p = input.cols;

    if (m == 2)
        {
        add<T>(input.rowptr(0), input.rowptr(1), p, out.rowptr(0));
        fused_rollback_add<T>(input.rowptr(0), input.rowptr(1), p, 1, out.rowptr(1));
        return;
        }
    else if (m == 1)
        {
        memcpy(out.data, input.data, p * sizeof(T));
        return;
        }

    transform<T>(input.head(), out.head(), temp.head());
    transform<T>(input.tail(), out.tail(), temp.tail());
    merge<T>(temp.head().as_const(), temp.tail().as_const(), out);
    }


template <typename T>
Block<T> transform(const T* input, size_t rows, size_t cols, T* temp, T* out)
    {
    transform<T>(
        ConstBlock<T>(input, rows, cols),
        Block<T>(temp, rows, cols),
        Block<T>(out, rows, cols)
        );
    return Block<T>(out, rows, cols);
    }


} // namespace riptide

#endif // TRANSFORMS_HPP
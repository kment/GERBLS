/*
 * physfunc.cpp
 *
 * 	A collection of useful physics/math functions
 *
 *  Created on: Aug 22, 2017
 *      Author: Kristo Ment
 */

#include "physfunc.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

// Bin data by period P
// Bin times will be [0,P/N_bins), ..., [P-P/N_bins,P)
void bin(double P,
         size_t N_bins,
         const DataContainer *data,
         double *bin_mags,
         double *bin_errors,
         size_t *N_bins_real)
{
    size_t i;
    *N_bins_real = 0;

    for (i = 0; i < N_bins; i++) {
        bin_mags[i] = 0;
        bin_errors[i] = 0;
    }

    for (size_t j = 0; j < data->size; j++) {
        i = (int)((fmod(data->rjd[j], P) / P) * N_bins);
        bin_mags[i] += data->mag_frac[j] / SQ(data->err_frac[j]);
        bin_errors[i] += 1 / SQ(data->err_frac[j]);
    }

    for (i = 0; i < N_bins; i++) {
        if (bin_errors[i]) {
            bin_mags[i] /= bin_errors[i];
            bin_errors[i] = sqrt(1 / bin_errors[i]);
            (*N_bins_real)++;
        }
        else {
            bin_errors[i] = 1000000000;
        }
    }
}

// Draw a sample from a distribution a <= x <= b such that p(x) = C * x^q
double gdraw(double a, double b, double q)
{
    double x;
    double C = 1 / std::max(pow(a, q), pow(b, q));
    do {
        x = grand(a, b);
    } while (C * pow(x, q) >= grand(0, 1));
    return x;
}

// Calculate the a/R ratio given P [days], M [Msun], R [Rsun]
double get_aR_ratio(double P, double M, double R)
{
    return pow(M, 0.3333) * pow(P / 365.25, 0.6667) / R / 0.00465;
}

// Solve Kepler's eqn M = E - e * sin(E) for E using Newton's method
// tol is the desired relative accuracy
double get_eccentric(double M, double e, double tol)
{
    double E0 = (e > 0.8 ? PI : M);
    double acc = 1;
    double E;
    int max_tries = 1000;

    for (int i = 0; i < max_tries; i++) {
        E = E0 - (E0 - e * std::sin(E0) - M) / (1 - e * std::cos(E0));
        acc = std::abs((E - E0) / (E != 0 ? E : 1));
        E0 = E;
        if (acc <= tol)
            return E;
    }

    std::cerr << "Required accuracy for E not achieved: " << acc << " > " << tol << "\n";
    std::cerr << "M = " << M << " e = " << e << "\n";
    return E;
}

// Calculate the orbital inclination (in deg) based on the impact parameter
double get_inc(double P, double M, double R, double b)
{
    return acos(b / get_aR_ratio(P, M, R)) * 180 / PI;
}

double get_inc(double aR_ratio, double b)
{
    return acos(b / aR_ratio) * 180 / PI;
}

// Estimate the min and max transit duration as a fraction of the orbital period P
void get_phase_range(double P, double *phase_min, double *phase_max)
{
    *phase_min = acos(1 - 0.002 / pow(P, 1.3333)) / TWOPI;
    *phase_max = acos(1 - 0.1 / pow(P, 1.3333)) / TWOPI;
}

// Estimate the transit duration given orbital period P (in days)
// Depends on the stellar mass M [Msun] and radius R [Rsun], and the impact parameter b
double get_transit_dur(double P, double M, double R, double b)
{
    double aR_ratio = get_aR_ratio(P, M, R);
    double inc = get_inc(aR_ratio, b) * PI / 180;
    return P * asin(sqrt(1 - SQ(b)) / aR_ratio / sin(inc)) / PI;
}

// Generate a random number between a and b
double grand(double a, double b)
{
    return a + (b - a) * rand() / RAND_MAX;
}

// Find the arithmetic mean of an array or a vector
double mean(const double *arr, size_t size)
{
    return sum(arr, size) / size;
}

double mean(const std::vector<double> &v)
{
    return mean(v.data(), v.size());
}

// Find the median of a vector
// Needs to make a copy because of std::nth_element()
double median(std::vector<double> v)
{
    // If size is even
    if (v.size() % 2 == 0) {
        std::nth_element(v.begin(), v.begin() + v.size() / 2 - 1, v.end());
        double a = v[v.size() / 2 - 1];
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        double b = v[v.size() / 2];
        return (a + b) / 2;
    }
    // If size is odd
    else {
        std::nth_element(v.begin(), v.begin() + (v.size() - 1) / 2, v.end());
        return v[(v.size() - 1) / 2];
    }
}

// Find the sum of an array or a vector
double sum(const double *arr, size_t size)
{
    double S = 0;
    for (size_t i = 0; i < size; i++) S += *(arr++);
    return S;
}

double sum(const std::vector<double> &v)
{
    return sum(v.data(), v.size());
}

// Squared trigonometric functions for convenience
double sin2(double x)
{
    return sin(x) * sin(x);
}

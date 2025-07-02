/*
 * physfunc.hpp
 *
 * 	A collection of useful physics/math functions
 *
 *  Created on: Aug 22, 2017
 *      Author: Kristo Ment
 */

#ifndef PHYSFUNC_HPP_
#define PHYSFUNC_HPP_

#include "structure.hpp"

#define PI 3.14159265359
#define TWOPI 6.28318530718
#define LOG10E 0.43429448190
#define SQ(x) ((x) * (x))

// Bin measurements by period
void bin(double P,
         size_t N_bins,
         const DataContainer *data,
         double *bin_mags,
         double *bin_errors,
         size_t *N_bins_real);

// Draw a random number from a power prior distribution
double gdraw(double a, double b, double q);

// Calculate the a/R ratio
double get_aR_ratio(double P, double M, double R);

// Solve Kepler's equation for the eccentric anomaly
double get_eccentric(double M, double e, double tol = 1e-7);

// Calculate the inclination in degrees
double get_inc(double P, double M, double R, double b);
double get_inc(double aR_ratio, double b);

// Estimate the min and max transit duration as a fraction of the orbital period P
void get_phase_range(double P, double *phase_min, double *phase_max);

// Estimate the transit duration given orbital period P (in days)
double get_transit_dur(double P, double M = 1, double R = 1, double b = 0);

// Generate a random number between a and b
double grand(double a, double b);

// Find the arithmetic mean of an array or a vector
double mean(const double *arr, size_t size);
double mean(const std::vector<double> &v);

// Find the median of a vector
double median(std::vector<double> v);

// Find the sum of an array or a vector
double sum(const double *arr, size_t size);
double sum(const std::vector<double> &v);

// Squared trigonometric functions for convenience
double sin2(double x);

#endif /* PHYSFUNC_HPP_ */

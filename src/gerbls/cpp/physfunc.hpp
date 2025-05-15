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
void bin(double, int, const DataContainer*, double*, double*, int*);

// Draw a random number from a power prior distribution
double gdraw(double, double, double);

// Calculate the a/R ratio
double get_aR_ratio(double, double, double);

// Solve Kepler's equation for the eccentric anomaly
double get_eccentric(double, double, double=1e-7);

// Calculate the inclination in degrees
double get_inc(double, double, double, double);
double get_inc(double, double);

// Estimate the min and max transit duration as a fraction of the orbital period P
void get_phase_range(double, double*, double*);

// Estimate the transit duration given orbital period P (in days)
double get_transit_dur(double, double=1, double=1, double=0);

// Generate a random number between a and b
double grand(double, double);

// Find the median of a vector
double median(std::vector<double>);

// Find the arithmetic mean of an array or a vector
double mean(const double*, size_t);
double mean(const std::vector<double>&);

// Find the sum of an array or a vector
double sum(const double*, size_t);
double sum(const std::vector<double>&);

// Squared trigonometric functions for convenience
double sin2(double);

#endif /* PHYSFUNC_HPP_ */

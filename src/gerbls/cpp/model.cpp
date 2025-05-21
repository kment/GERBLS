/*
 * model.cpp
 *
 *  Created on: Aug 20, 2017
 *      Author: Kristo Ment
 */

#include "ffafunc.hpp"
//#include "interpolation.h" 	// ALGLIB dependency
#include "model.hpp"
#include "physfunc.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>

// Constructor
// If f_min or f_max are 0 then the default values are used
// If no target is given, use default values
BLSModel::BLSModel(DataContainer& data_ref, double f_min_, double f_max_, const Target* targetPtr) {

	data = &data_ref;

	// Set min and max frequencies
	if (f_min_ > 0)
		f_min = f_min_;
	if (f_max_ > 0)
		f_max = f_max_;

    if (targetPtr == nullptr)
        target = new Target();
    else target = targetPtr;

}

// Get number of frequencies
size_t BLSModel::N_freq() {
	return freq.size();
}

void BLSModel::run(bool verbose) {
	std::cout << "run() is not defined for an object of type " << typeid(*this).name();
}

// Constructor
// If f_min or f_max are 0 then the default values are used
// If no target is given, use default values
BLSModel_FFA::BLSModel_FFA(DataContainer& data_ref, double f_min_, double f_max_, Target* targetPtr) :
						  BLSModel(data_ref, f_min_, f_max_, targetPtr) { };

// Generate required results
template <typename T>
void BLSModel_FFA::process_results(std::vector<BLSResult<T>>& results) {

	const size_t N_freq = results.size();
	BLSResult<T>* pres = results.data();

	freq.resize(N_freq);
	dchi2.assign(N_freq, 0);
	chi2_mag0.assign(N_freq, 0);
	chi2_dmag.assign(N_freq, 0);
	chi2_t0.assign(N_freq, 0);
	chi2_dt.assign(N_freq, 0);

	for (size_t i = 0; i < N_freq; i++) {
		freq[i] = 1 / pres->P;
		dchi2[i] = -(pres->dchi2);
		chi2_mag0[i] = pres->mag0;
		chi2_dmag[i] = pres->dmag;
		chi2_t0[i] = fmod(rdata->rjd[0] + pres->P * (pres->t0 - 0.5) / pres->N_bins, pres->P);
		chi2_dt[i] = t_samp * pres->dur;
		pres++;
	}

}

void BLSModel_FFA::run(bool verbose) {
	run_prec<float>(verbose);
}

void BLSModel_FFA::run_double(bool verbose) {
	run_prec<double>(verbose);
}

// Data will be resampled uniformly to cadence tsamp
template <typename T>
void BLSModel_FFA::run_prec(bool verbose) {

	if (verbose)
		std::cout << "Starting FFA...\n";

	// Resample to desired tsamp
	rdata = resample_uniform(*data, t_samp);
	std::vector<T> mag(rdata->size, 0);	// Magnitudes
	std::vector<T> wts(rdata->size, 0);	// Weights (1/err^2)
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

	auto t_start = std::chrono::high_resolution_clock::now();
	std::vector<BLSResult<T>> pgram = std::move(periodogram<T>(mag.data(), wts.data(), mag.size(), t_samp, *target, 1/f_max, 1/f_min));
	auto t_end = std::chrono::high_resolution_clock::now();

	if (verbose) {
		std::chrono::duration<double> rtime = t_end - t_start;
		std::cout << "Number of tested periods: " << pgram.size() << "\n";
		std::cout << "BLS runtime: " << rtime.count() << " sec\n";
	}
	process_results(pgram);

}
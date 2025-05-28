/*
 * model.hpp
 *
 *  Created on: Aug 20, 2017
 *      Author: Kristo Ment
 */

#ifndef MODEL_HPP_
#define MODEL_HPP_

#include "ffafunc.hpp"
#include "structure.hpp"
#include <fstream>
#include <unordered_map>

// BLS model (base class)
struct BLSModel {

	// Settings
	double f_min = 0.025;
	double f_max = 5;
	int max_duration_mode = 2;				// Affects get_max_duration()
	double max_duration_factor = 0.1;		// Affects get_max_duration()

	// Pointer to associated data
	DataContainer* data = nullptr;
	
	// Pointer to associated target
	const Target* target = nullptr;

	// Array to store tested frequencies
	std::vector<double> freq;

	// Constructor and destructor
	BLSModel(DataContainer&, double=0, double=0, const Target* =nullptr, int=0, double=0);
	virtual ~BLSModel() = default;

	double get_max_duration(double);		// Maximum transit duration to test at a given period
	size_t N_freq();	// Get number of frequencies

	// Virtual functions to be overwritten
	virtual void run(bool);

	// Required results for each tested frequency
	std::vector<double> dchi2, chi2_mag0, chi2_dmag, chi2_t0, chi2_dt;

};

// BLS model (brute force) removed

// BLS model (FFA)
struct BLSModel_FFA: public BLSModel {

	// Settings
	double t_samp = 2./60/24; 		// Uniform cadence to resample data to

	// Pointer to the resampled data
	std::unique_ptr<DataContainer> rdata;

	// TEMPORARY FFA results
	std::vector<double> periods;
	std::vector<size_t> widths;
	std::vector<size_t> foldbins;
	std::vector<double> snr;
	std::vector<size_t> t0;

	// Inherit constructor from parent
	using BLSModel::BLSModel;

	// Methods
	template <typename T> void process_results(std::vector<BLSResult<T>>&);
	void run(bool=true);
	void run_double(bool=true);
	template <typename T> void run_prec(bool=true);

};

#endif /* MODEL_HPP_ */

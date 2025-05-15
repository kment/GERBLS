/*
 * structure.hpp
 *
 *  Created on: Aug 16, 2017
 *      Author: Kristo Ment
 */

#ifndef STRUCTURE_HPP_
#define STRUCTURE_HPP_

#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

// New array-structure to hold data (dynamically allocated)
struct DataContainer {

    int* sec = nullptr;       // Data sector
    double* rjd;
    double* mag;
    double* err;
    size_t size = 0;
    bool alloc = false;
    
    // Optional extra variables
    double mag_avg = 0;
    std::unique_ptr<int[]> sectors;
    std::unique_ptr<double[]> mag_frac;
    std::unique_ptr<double[]> err_frac;
    std::unique_ptr<bool[]> valid_mask;   // Only valid data evaluates to true
    
    // Destructor
    ~DataContainer();
    
    void allocate(size_t);
    void calculate_mag_frac();
    std::unique_ptr<DataContainer> clean(double=0, bool* =nullptr, int=3);
    std::unique_ptr<DataContainer> clean_hw(double, bool* =nullptr, int=3);
    std::unique_ptr<bool[]> find_flares(const double*);
    std::unique_ptr<bool[]> find_flares();
    std::set<int> get_sectors();
    double get_time_range() const;
    void imprint(double*, size_t);
    std::unique_ptr<DataContainer> phase_folded(double, double=0, int* =nullptr);
    void read_from_file(std::string, const std::vector<int>* =nullptr);
    std::vector<double> running_median(double);
    std::vector<double> running_median_eval(double, double*, size_t);
    std::vector<double> running_median_per(double, double);
    void set(double*, double*, double*, size_t);
    void set(int*, double*, double*, double*, size_t);
    std::vector<double> splfit(double, int=50);
    std::vector<double> splfit_eval(int, double*, size_t);
    std::unordered_map<int, std::unique_ptr<DataContainer>> split_by_sector();
    void store(double*, double*, double*, size_t);
    void store(int*, double*, double*, double*, size_t);
    void store_unmasked(double*, double*, double*, bool*, size_t);
    void store_unmasked(int*, double*, double*, double*, bool*, size_t);

};

// Count lines in a data file
size_t file_count_lines(std::string, const std::vector<int>* =nullptr);

// Structure to keep information about the star
// Uses MRL values for LHS 1140 by default
struct Target {

    double M = 0.179;       // Mass in Msun
    double R = 0.209;       // Radius in Run
    double L = 0.00433;     // Luminosity in Lsun
    double u1 = 0;          // Limb darkening coefficients
    double u2 = 0;
    double L_comp = 0;      // Luminosity of binary companion as a fraction of L
    double P_rot = 0;       // Rotation period
    double P_rot2 = 0;      // Second rotation period (binary companion)
    
    double logg();
    double Teff();

};

#endif /* STRUCTURE_HPP_ */

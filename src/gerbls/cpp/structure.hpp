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
    int *sec = nullptr; // Data sector
    double *rjd;
    double *mag;
    double *err;
    size_t size = 0;
    bool alloc = false;

    // Optional extra variables
    double mag_avg = 0;
    std::unique_ptr<int[]> sectors;
    std::unique_ptr<double[]> mag_frac;
    std::unique_ptr<double[]> err_frac;
    std::unique_ptr<bool[]> valid_mask; // Only valid data evaluates to true

    // Destructor
    ~DataContainer();

    void allocate(size_t size);
    void calculate_mag_frac();
    std::unique_ptr<DataContainer> clean(double P_rot = 0, bool *mask = nullptr, int N_flares = 3);
    std::unique_ptr<DataContainer> clean_hw(double hw, bool *mask = nullptr, int N_flares = 3);
    std::unique_ptr<bool[]> find_flares(const double *mag0);
    std::unique_ptr<bool[]> find_flares();
    std::set<int> get_sectors();
    double get_time_range() const;
    void imprint(double *signal, size_t size);
    std::unique_ptr<DataContainer>
        phase_folded(double P_rot, double t_extend = 0, int *save_indices = nullptr);
    void read_from_file(std::string filename, const std::vector<int> *sectors = nullptr);
    std::vector<double> running_median(double hw);
    std::vector<double> running_median_eval(double hw, double *rjd_eval, size_t N_eval);
    std::vector<double> running_median_per(double hw, double P_rot);
    void set(double *rjd, double *mag, double *err, size_t size);
    void set(int *sec, double *rjd, double *mag, double *err, size_t size);
    // std::vector<double> splfit(double, int = 50);
    // std::vector<double> splfit_eval(int, double *, size_t);
    std::unordered_map<int, std::unique_ptr<DataContainer>> split_by_sector();
    void store(double *rjd, double *mag, double *err, size_t size);
    void store(int *sec, double *rjd, double *mag, double *err, size_t size);
    void store_unmasked(double *rjd, double *mag, double *err, bool *mask, size_t size);
    void store_unmasked(int *sec, double *rjd, double *mag, double *err, bool *mask, size_t size);
};

// Count lines in a data file
size_t file_count_lines(std::string filename, const std::vector<int> *sectors = nullptr);

// Structure to keep information about the star
// Uses MRL values for LHS 1140 by default
struct Target {
    double M = 0.179;   // Mass in Msun
    double R = 0.209;   // Radius in Run
    double L = 0.00433; // Luminosity in Lsun
    double u1 = 0;      // Limb darkening coefficients
    double u2 = 0;
    double L_comp = 0; // Luminosity of binary companion as a fraction of L
    double P_rot = 0;  // Rotation period
    double P_rot2 = 0; // Second rotation period (binary companion)

    double logg();
    double Teff();
};

#endif /* STRUCTURE_HPP_ */

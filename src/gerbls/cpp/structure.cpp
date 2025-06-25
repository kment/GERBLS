/*
 * structure.cpp
 *
 *  Created on: Aug 16, 2017
 *      Author: Kristo Ment
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>
// #include "interpolation.h"    // ALGLIB dependency
#include "physfunc.hpp"
#include "structure.hpp"

DataContainer::~DataContainer()
{
    if (alloc) {
        delete[] sec;
        delete[] rjd;
        delete[] mag;
        delete[] err;
    }
}

// Allocate space for data (no initialization)
void DataContainer::allocate(size_t _size)
{
    if (size > 0)
        std::cerr << "Error: DataContainer is not empty.\n";
    else {
        sec = new int[_size];
        rjd = new double[_size];
        mag = new double[_size];
        err = new double[_size];
        size = _size;
        alloc = true;
    }
}

// Calculate the fractional magnitude differences
void DataContainer::calculate_mag_frac()
{
    if (!mag_frac) {
        mag_frac.reset(new double[size]);
        err_frac.reset(new double[size]);
    }

    mag_avg = 0;
    double weights_sum = 0;

    for (size_t i = 0; i < size; i++) {
        mag_avg += mag[i] / SQ(err[i]);
        weights_sum += 1 / SQ(err[i]);
    }

    mag_avg /= weights_sum;

    for (size_t i = 0; i < size; i++) {
        mag_frac[i] = mag[i] / mag_avg - 1;
        err_frac[i] = err[i] / mag_avg;
    }
}

// Return a cleaned version of the data set (baseline and flares removed)
// If given, rotation period P_rot alters the running median window
// If given, mask of removed values will be stored in mask
// Iterative flare removal will be done N_flares times
std::unique_ptr<DataContainer> DataContainer::clean(double P_rot, bool *mask, int N_flares)
{
    double hw = (P_rot ? std::max(P_rot / 20, 2.5 / 24) : 6. / 24);
    return clean_hw(hw, mask, N_flares);
}

std::unique_ptr<DataContainer> DataContainer::clean_hw(double hw, bool *mask, int N_flares)
{
    std::unique_ptr<DataContainer> data;
    std::unique_ptr<bool[]> mask_;
    std::vector<double> med = running_median(hw);

    if (N_flares == 0) {
        data.reset(new DataContainer);
        data->store(sec, rjd, mag, err, size);
        mask_.reset(new bool[size] {false});
    }
    else {
        for (int i = 0; i < N_flares; i++) {
            mask_ = find_flares(med.data());
            data.reset(new DataContainer);
            data->store_unmasked(sec, rjd, mag, err, mask_.get(), size);
            med = data->running_median_eval(hw, rjd, size);
        }
        med = data->running_median(hw);
    }

    for (size_t i = 0; i < data->size; i++) data->mag[i] /= med[i];

    if (mask != nullptr)
        std::copy_n(mask_.get(), size, mask);
    return data;
}

// Return a boolean mask indicating any flare data
// Flares are defined as a series (4) of consequtive bright 2-sigma deviants (in magn.)
// double[size] mag0 is the baseline flux to be divided out before detection
// Note that (m-m0 < -a*dm) is equivalent to log(F/F0) > a*log(e)*dF/F
std::unique_ptr<bool[]> DataContainer::find_flares(const double *mag0)
{
    std::unique_ptr<bool[]> mask(new bool[size] {false});
    size_t j, k;

    for (size_t i = 0; i < size; i++) {
        j = 0;
        while ((i + j < size) &&
               (log10(mag[i + j] / mag0[i + j]) >= 2 * LOG10E * err[i + j] / mag[i + j]))
            j++;
        if (j >= 4) {
            k = i;
            while ((k > 0) && (mag[k - 1] > mag0[k - 1])) k--;
            while ((i + j < size) && (mag[i + j] > mag0[i + j])) j++;
            for (; k < i + j; k++) mask[k] = true;
            i += j;
        }
    }

    return mask;
}

std::unique_ptr<bool[]> DataContainer::find_flares()
{
    double mag0[size];
    for (size_t i = 0; i < size; i++) mag0[i] = 1;
    return find_flares(mag0);
}

// Get a sorted list of unique sectors in the data
std::set<int> DataContainer::get_sectors()
{
    if (sec == nullptr)
        std::cerr << "Data has no sector information.";

    std::set<int> sectors;
    for (size_t i = 0; i < size; i++) sectors.emplace(sec[i]);

    return sectors;
}

// Calculate the time range in days
double DataContainer::get_time_range() const
{
    if (!size)
        return 0;

    double t_min = rjd[0];
    double t_max = t_min;

    for (size_t i = 0; i < size; i++) {
        t_min = std::min(t_min, rjd[i]);
        t_max = std::max(t_max, rjd[i]);
    }

    return (t_max - t_min);
}

// Imprint a signal to existing data
// Assumes the signal mags are in the same order (and size) as mag
void DataContainer::imprint(double *signal, size_t size_)
{
    if (size_ != size) {
        std::cerr << "Failed to insert signal into LC: signal array has the wrong size";
        return;
    }

    for (size_t i = 0; i < size; i++) mag[i] *= signal[i];
}

// Return a phase-folded data container sorted by phase
// If save_indices is given, store indices in that array
// If t_extend is given, extend for t<0 and t>P_rot by duplicating periodic values
std::unique_ptr<DataContainer>
    DataContainer::phase_folded(double P_rot, double t_extend, int *save_indices)
{
    std::unique_ptr<DataContainer> data(new DataContainer());
    double phases[size];
    int indices[size];

    for (size_t i = 0; i < size; i++) {
        phases[i] = fmod(rjd[i], P_rot);
        indices[i] = i;
    }

    std::sort(indices, indices + size, [&phases](int j, int k) { return (phases[j] < phases[k]); });

    // Number of data points to extend
    size_t j = 0, k = 0;
    while ((j < size) && (phases[indices[j]] < t_extend)) j++;
    while ((k < size) && (phases[indices[size - k - 1]] > P_rot - t_extend)) k++;

    data->rjd = new double[size + j + k];
    data->mag = new double[size + j + k];
    data->err = new double[size + j + k];
    data->size = size + j + k;
    data->alloc = true;
    for (size_t i = 0; i < size; i++) {
        data->rjd[i + k] = phases[indices[i]];
        data->mag[i + k] = mag[indices[i]];
        data->err[i + k] = err[indices[i]];
    }
    if (j > 0) {
        for (size_t i = 0; i < j; i++) data->rjd[size + i + k] = data->rjd[i + k] + P_rot;
        std::copy_n(data->mag + k, j, data->mag + size + k);
        std::copy_n(data->err + k, j, data->err + size + k);
    }
    if (k > 0) {
        for (size_t i = 0; i < k; i++) data->rjd[i] = data->rjd[size + i] - P_rot;
        std::copy_n(data->mag + size, k, data->mag);
        std::copy_n(data->err + size, k, data->err);
    }

    if (save_indices)
        std::copy_n(indices, size, save_indices);

    return data;
}

void DataContainer::read_from_file(std::string filename, const std::vector<int> *sectors)
{
    if (size != 0) {
        std::cerr << "Error reading from file: data container must be empty.";
        return;
    }

    // Open input data file
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Error: Unable to open " << filename << "!";
        return;
    }

    int sector;
    std::string line;
    std::istringstream linestream;

    size_t N_lines = file_count_lines(filename, sectors);
    allocate(N_lines);

    // Read in data
    for (size_t i = 0; i < N_lines;) {
        std::getline(fin, line);

        if (fin.eof()) {
            std::cout
                << "Something went wrong: end of data file reached before all data read in.\n";
            break;
        }

        if ((line[0] == '#') || (line.empty()))
            continue;

        linestream.str(line);
        linestream.clear();
        linestream >> sector;
        if (sectors) {
            if (std::find(sectors->begin(), sectors->end(), sector) != sectors->end()) {
                sec[i] = sector;
                linestream >> rjd[i] >> mag[i] >> err[i];
                i++;
            }
        }
        else {
            sec[i] = sector;
            linestream >> rjd[i] >> mag[i] >> err[i];
            i++;
        }
    }

    fin.close();
    std::cout << "Read in " << N_lines << " data points.\n";
    return;
}

// Running median with a half-width hw, evaluated at data points
// Assumes that the data is time-sorted
std::vector<double> DataContainer::running_median(double hw)
{
    return running_median_eval(hw, rjd, size);
}

// Running median with a half-width hw, evaluated at times rjd_eval
// Assumes that both the data and rjd_eval are time-sorted
std::vector<double> DataContainer::running_median_eval(double hw, double *rjd_eval, size_t N_eval)
{
    std::vector<double> med(N_eval, 0);
    size_t i = 0, j = 0;

    for (size_t k = 0; k < N_eval; k++) {
        while ((i < size) && (rjd[i] < rjd_eval[k] - hw)) i++;
        while ((j < size) && (rjd[j] < rjd_eval[k] + hw)) j++;
        if (i < j)
            med[k] = median(std::vector<double>(mag + i, mag + j));
    }

    return med;
}

// Periodic running median (using phase-folded data)
std::vector<double> DataContainer::running_median_per(double hw, double P_rot)
{
    std::vector<double> med(size);
    int indices[size];
    int N_extend = 0;

    std::unique_ptr<DataContainer> data = phase_folded(P_rot, hw, indices);
    while (data->rjd[N_extend] < 0) N_extend++;
    std::vector<double> med_per = data->running_median(hw);
    for (size_t i = 0; i < size; i++) med[indices[i]] = med_per[i + N_extend];

    return med;
}

void DataContainer::set(double *rjd, double *mag, double *err, size_t size)
{
    set(nullptr, rjd, mag, err, size);
}

void DataContainer::set(int *sec, double *rjd, double *mag, double *err, size_t size)
{
    this->sec = sec;
    this->rjd = rjd;
    this->mag = mag;
    this->err = err;
    this->size = size;
}

// Periodic cubic spline fit
// Removed due to ALGLIB dependency
/*std::vector<double> DataContainer::splfit(double P_rot, int M) {

    std::vector<double> spl_(size);
    alglib::real_1d_array x, y;
    alglib::spline1dinterpolant spl;
    alglib::spline1dfitreport rep;

    x.setlength(size);
    y.setcontent(size, mag);

    for (size_t i = 0; i < size; i++)
        x[i] = fmod(rjd[i], P_rot);

    alglib::spline1dfitcubic(x, y, size, M, spl, rep);
    for (size_t i = 0; i < size; i++)
        spl_[i] = alglib::spline1dcalc(spl, x[i]);

    return spl_;

}*/

// Cubic spline fit with M nodes evaluated at rjd_eval
// Removed due to ALGLIB dependency
/*std::vector<double> DataContainer::splfit_eval(int M, double* rjd_eval, size_t N_eval) {

    std::vector<double> spl_(N_eval);
    alglib::real_1d_array x, y;
    alglib::spline1dinterpolant spl;
    alglib::spline1dfitreport rep;

    x.setcontent(size, rjd);
    y.setcontent(size, mag);
    alglib::spline1dfitcubic(x, y, size, M, spl, rep);

    for (size_t i = 0; i < N_eval; i++)
        spl_[i] = alglib::spline1dcalc(spl, rjd_eval[i]);

    return spl_;

}*/

// Split data by sector (returns a map of sector -> data)
std::unordered_map<int, std::unique_ptr<DataContainer>> DataContainer::split_by_sector()
{
    std::unique_ptr<bool[]> mask(new bool[size]);
    std::set<int> sectors(get_sectors());
    std::unordered_map<int, std::unique_ptr<DataContainer>> data_map;

    for (auto &sector : sectors) {
        DataContainer *pdata = data_map.emplace(sector, new DataContainer).first->second.get();
        for (size_t i = 0; i < size; i++) mask[i] = (sec[i] != sector);
        pdata->store_unmasked(sec, rjd, mag, err, mask.get(), size);
    }

    return data_map;
}

// Allocate memory and copy a dataset (sec filled with zeros)
void DataContainer::store(double *rjd, double *mag, double *err, size_t size)
{
    store(nullptr, rjd, mag, err, size);
}

// Allocate memory and copy a dataset
void DataContainer::store(int *sec, double *rjd, double *mag, double *err, size_t size)
{
    allocate(size);
    if (sec == nullptr)
        for (size_t i = 0; i < size; i++) this->sec[i] = 0;
    else
        std::copy_n(sec, size, this->sec);
    std::copy_n(rjd, size, this->rjd);
    std::copy_n(mag, size, this->mag);
    std::copy_n(err, size, this->err);
}

// Equivalent to store but only stores members where mask evaluates to false
void DataContainer::store_unmasked(double *rjd, double *mag, double *err, bool *mask, size_t size)
{
    store_unmasked(nullptr, rjd, mag, err, mask, size);
}

void DataContainer::store_unmasked(
    int *sec, double *rjd, double *mag, double *err, bool *mask, size_t size)
{
    size_t new_size = 0;
    for (size_t i = 0; i < size; i++) new_size += (!mask[i]);
    allocate(new_size);
    for (size_t i = 0, j = 0; i < size; i++) {
        if (!mask[i]) {
            this->sec[j] = (sec == nullptr ? 0 : sec[i]);
            this->rjd[j] = rjd[i];
            this->mag[j] = mag[i];
            this->err[j] = err[i];
            j++;
        }
    }
}

// Count lines in a file
size_t file_count_lines(std::string filename, const std::vector<int> *sectors)
{
    // Open input data file
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Error: Unable to open " << filename << "!";
        return 0;
    }

    size_t N_lines = 0;
    int sector;
    std::string line;
    std::istringstream linestream;

    while (!fin.eof()) {
        std::getline(fin, line);
        if ((line[0] == '#') || (line.empty()))
            continue;
        if (sectors) {
            linestream.str(line);
            linestream.clear();
            linestream >> sector;
            if (std::find(sectors->begin(), sectors->end(), sector) != sectors->end())
                N_lines++;
        }
        else
            N_lines++;
    }

    return N_lines;
}

// Estimate the surface gravity (in cm/s2)
double Target::logg()
{
    return log10(M) - 2 * log10(R) + 4.437;
}

// Estimate the effective temperature
double Target::Teff()
{
    return 5772 * pow(L, 0.25) / pow(R, 0.5);
}
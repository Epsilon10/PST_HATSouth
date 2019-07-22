#pragma once

#include <iostream>
#include <cmath>
#include <vector> 
#include <numeric> 

/**
 * Optical transfer function
 * Represents the fourier transform of the airy function
 * @param roh_prime - equal to roh / roh_c
 */ 
double get_otf(double roh_prime) {
    const double c1 = 2 / M_PI;
    if (roh_prime >= 1.0) return 0.0;
    const double inside = acos(roh_prime) - roh_prime*sqrt(1 - pow(roh_prime, 2));
    return c1 * inside;
}

/**
 * Returns the standard deviation
 * @param input - input data
 * @param mean - mean of the data
 */
double get_variance(const std::vector<double>& input, double mean) {
    const double N = input.size();
    double res = 0.0;
    for (double x : input) 
        res += (x - mean) * (x - mean);
    res *= 1.0 / N;
    return sqrt(res);
}

/**
 * Returns the value of a fourier transformed gaussian function for some input and value
 * @param variance
 * @param x - wavelength to evaluate the gaussian at
 */ 
double fourier_gaussian1d(double variance, double x) {
    const double K = 1.0 / sqrt(2 * M_PI * (variance * variance));
    const double a = 1.0 / (2 * (variance * variance));
    return K * sqrt(M_PI / a) * std::exp(-1 *(M_PI * M_PI) * (x * x) / a);
}

/**
 * Returns the value of the point spread function for some x and y
 * @param wavelength - wavelength being shined on the spots
 * @param f number - focal ration
 * @param k_x - x wavelength
 * @param k_y - y wavelength
 * @param s_x - standard deviation of the x coords of centroids
 * @param s_y - standard deviation of the y coords of centroids
 */ 
double fourier_psf(double wavelength, double f_number, double k_x, double k_y, 
       double s_x,double s_y) {
    const double roh = sqrt(k_x * k_x + k_y * k_y);
    const double roh_c = 1 / (f_number * wavelength);
    const double otf_val = get_otf (roh / roh_c);
    const double x_gauss = fourier_gaussian1d(s_x, k_x);
    const double y_gauss = fourier_gaussian1d(s_y, k_y);
    // std::cout << otf_val << std::endl;
    return otf_val * x_gauss * y_gauss;
}
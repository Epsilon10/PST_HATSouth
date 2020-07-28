#pragma once

#include "include/intensity.hh"
#include <utility>
#include <cstdint>
#include "include/alglib/src/interpolation.h"
#include <vector>

class PRF {
    public:
        explicit PRF(
            double wavelength, double focal_ratio, 
            std::pair<double,double> variance,
            std::pair<double,double> standard_deviation,
            std::pair<std::uint32_t, std::uint32_t> subpixel_res, 
            std::pair<std::uint32_t, std::uint32_t> map_res
        );
    
        void initialize(std::pair<double,double> max_wavenumber);
    private:
        Intensity intensity;
        alglib::spline2dinterpolant spline;
        std::vector<double> x_vals, y_vals;
};
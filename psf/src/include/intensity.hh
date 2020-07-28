#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include "psf.hh"
#include "tensor.hh"
#include <algorithm>
#include <utility>

constexpr std::complex<double> i(0.0,1.0);

/**
 * Class that 
 *
 */
class Intensity {
    public:
    Intensity(
        double wavelength, double focal_ratio, std::pair<double,double> standard_deviations, std::pair<uint32_t, uint32_t> subpixel_resolutions,
            std::pair<uint32_t, uint32_t> map_resolutions
        ) : wavelength(wavelength), focal_ratio(focal_ratio), standard_deviations(standard_deviations), subpixel_resolutions(subpixel_resolutions), 
            map_resolutions(map_resolutions), dx(1.0 / subpixel_resolutions.first), dy(1.0 / subpixel_resolutions.second) {

        I_vals.resize(map_resolutions.first * map_resolutions.second);
        psi_vals.resize(map_resolutions.first * map_resolutions.second);

        p = fftw_plan_r2r_2d(map_resolutions.first, map_resolutions.second, psi_vals.data(), 
                I_vals.data(), FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE
        );
        ky_vals.reserve(map_resolutions.second);
        kx_vals.reserve(map_resolutions.first);
        psi_amps.resize(psi_vals.size());
    }

    ~Intensity() {
        fftw_destroy_plan(p);
    }

    void execute(std::pair<double, double> wavenum_ranges) {
        this->wavenum_ranges = wavenum_ranges;
        fill_psi();
        fftw_execute(p);
    }

    /**
     * Returns the value of the integral at some i,j
     * @param i -> row
     * @param j -> column
     * TODO: switch to tensor implementation
     */
    double get_intensities(size_t i, size_t j) {
        return I_vals[j + i * map_resolutions.second];
    }

    /**
     * @return x wavenumbers
     */
    const std::vector<double>& get_x_wavenumbers_vals() const { return kx_vals; }

    /**
     * @return y wavenumbers
     */
    const std::vector<double>& get_y_wavenumbers_vals() const { return ky_vals; }

    /**
     * @return subpixel integrals
     */
    const std::vector<double>& get_all_intensities() const { return I_vals; }

    /**
     * @return psi amplitudes
     */
    const std::vector<double>& get_psi_amps() const { return psi_amps; }

    /**
     * Fill map x coordinates
     * @param in -> pointer to input array
     */
    void fill_x_wavenumbers(double* in) const {
        for (size_t k = 0; k < map_resolutions.first; k++) {
            double x = k / (2.0 * wavenum_ranges.first);
            in[k] = x;
        }
    }

    /**
     * Fill map y coordinates
     * @param in -> pointer to input array
     */
    void fill_y_wavenumbers(double* in) const {
        for (size_t l = 0; l < map_resolutions.second; l++) {
            double y =  l / (2.0 * wavenum_ranges.second);
            in[l] = y;
        }
    }

    /**
     * Get value of psi and a set of indices
     * @param i -> row index
     * @param -> column index
     */
    std::complex<double> get_psi_vals(size_t i, size_t j) const {
        return psi_vals[j + i * map_resolutions.second];
    }

    private:
    std::pair<double,double> standard_deviations;
    std::pair<uint32_t,uint32_t> subpixel_resolutions;
    std::pair<uint32_t, uint32_t> map_resolutions;
    std::pair<double, double> wavenum_ranges;

    double wavelength, focal_ratio;
    double dx;
    double dy;
    
    std::vector<double> I_vals;
    std::vector<double> psi_vals;
    std::vector<double> psi_amps;
    std::vector<double> kx_vals;
    std::vector<double> ky_vals;

    fftw_plan p;
    
    double psi(double m, double n, double k_x, double k_y) const {
        double sin_term = std::sin(M_PI * dx * k_x) * std::sin(M_PI * dx * k_y);
        double denom = M_PI * M_PI * (m + 0.5) * (n + 0.5);
        return fourier_psf(wavelength, focal_ratio, k_x, k_y, standard_deviations.first, standard_deviations.second) *
            sin_term / denom;
    }
    
    void fill_psi() {
        size_t m,n;
        for (m = 0; m < map_resolutions.first; m++) {
            double k_x = (m + 0.5) * (wavenum_ranges.first / map_resolutions.first);
            kx_vals.push_back(k_x);
            for (n = 0; n < map_resolutions.second; n++) {
                double k_y = (n + 0.5) * (wavenum_ranges.second / map_resolutions.second);
                if (ky_vals.size() < ky_vals.capacity()) 
                    ky_vals.push_back(k_y);

                double psi_res = psi(m, n, k_x, k_y);
                psi_vals[n + m * map_resolutions.second] = psi_res;
                psi_amps[n + m * map_resolutions.second] = std::abs(psi_res);
            }
        }
    }
     
};

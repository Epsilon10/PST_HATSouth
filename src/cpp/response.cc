#include "../include/response.hh"
#include <cstring>

extern "C" {
    void fill_I_vals(
        double wavelength, double N, double s_x, double s_y, uint32_t N_x, uint32_t N_y,
        uint32_t nu_x, uint32_t nu_y, double K_x, double K_y,  double* I_vals, double* kx_vals,
        double* ky_vals, double* psi_val_amps, double* x_vals, double* y_vals
    ) {
        Response response(wavelength, N,std::make_pair(s_x, s_y), std::make_pair(N_x, N_y), std::make_pair(nu_x, nu_y));
        response.get_response(std::make_pair(K_x, K_y));
        
        std::memcpy(I_vals, response.get_I_vals().data(), sizeof(double) * response.get_I_vals().size());
        std::memcpy(kx_vals, response.get_kx_vals().data(), sizeof(double) * response.get_kx_vals().size());
        std::memcpy(ky_vals, response.get_ky_vals().data(), sizeof(double) * response.get_ky_vals().size());
        std::memcpy(psi_val_amps, response.get_psi_vals_amp().data(), sizeof(double) * response.get_psi_vals_amp().size());
        response.fill_x_vals(x_vals);
        response.fill_y_vals(y_vals);
    }
}


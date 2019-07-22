#include <iostream>
#include<fftw3.h>
#include <vector>
#include <complex> 
#include <algorithm>

void fill_vals(std::vector<double>& vec);
int main() {
    std::vector<double> in_vals;
    in_vals.resize(100 * 100);
    std::vector<std::complex<double>> out_vals;
    out_vals.resize(100 * 51);

    fftw_plan p = fftw_plan_dft_r2c_2d(100, 100, in_vals.data(), reinterpret_cast<fftw_complex*>(out_vals.data()), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    fill_vals(in_vals);
    fftw_execute(p);
        std::vector<std::complex<double>> out_2(out_vals);

    fftw_plan p2 = fftw_plan_dft_c2r_2d(100, 100, reinterpret_cast<fftw_complex*>(out_vals.data()), in_vals.data(), FFTW_ESTIMATE);
    fftw_execute(p2);
    fftw_plan p3 = fftw_plan_dft_r2c_2d(100, 100, in_vals.data(), reinterpret_cast<fftw_complex*>(out_vals.data()), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    fftw_execute(p3);
    fftw_destroy_plan(p);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
    for (size_t i = 0; i < 10; i++) {
        for (size_t j = 0; j < 10; j++) {
            auto v1 = out_vals[j + i * 51];
            auto v2 = out_2[j + i * 51];
            std::cout << "first: " << v1 << " second: " << v2 << std::endl;
        }  
    }
}

void fill_vals(std::vector<double>& vec) {
    double dx = 0.1;
    double dy = 0.1;
    auto func = [](double x, double y) { return x * x + y * y; };
    for (size_t i = 0; i < 100; i++) {
        for (size_t j = 0; j < 100; j++) {
            vec[j + i * 100] = func(i*dx, j*dy);
        }
    }
}
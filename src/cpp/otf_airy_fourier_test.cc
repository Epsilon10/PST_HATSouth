#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <complex>
#include <ctime>
#include <boost/math/special_functions/bessel.hpp>
constexpr double PI = 3.141592653589793238463;

// aperture and setup constants
constexpr double a = 2.5;
constexpr double F = 10;
constexpr double N = F / (2 *a);
constexpr double wavelength = 0.6;
constexpr uint32_t max_intensity = 400;
constexpr double roh_c = 1 / (wavelength * N);


double u;
double v;
double roh;
int resolution;

std::complex<double> airy_fourier_numeric(double lower, double upper);
double otf(double roh_prime);
double bessel(double x);
int main(int argc, char** argv) {
    double lower =  std::atof(argv[1]);
    double upper =  std::atof(argv[2]);
    resolution = std::atof(argv[3]);
    u = std::atof(argv[4]);
    v = std::atof(argv[5]);
    roh = sqrt(u*u + v*v);
    
    std::clock_t begin = clock();
    std::complex<double> res = airy_fourier_numeric(lower, upper);

    std::cout << "u: " << u <<" v: " << v << std::endl;
    std::cout << "Lower bound: " << lower << " upper bound: " << upper << std::endl;
    std::cout << "wavelength: " << wavelength << " micrometers" << std::endl;
    std::cout << "roh_prime: " << roh / roh_c << std::endl;
    std::cout << "Numeric fourier result: " << res << std::endl;
    double otf_res = otf(roh / roh_c);
    std::cout << "otf value: " << otf_res<< std::endl;
    std::cout << "F(u,v) / otf: " << res.real() / otf_res << std::endl;
    std::clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Finished in: " << elapsed_secs << " seconds." << std::endl;
    
}

/**
 * Optical transfer function
 * @param roh_prime - equal to roh / roh_c
 */ 
double otf(double roh_prime) {
    double c1 = 2 / PI;
    double inside = acos(roh_prime) - roh_prime*sqrt(1 - pow(roh_prime, 2));
    return c1 * inside;
}

/**
 * Numerically calculates the bessel function of order 1 [J_1(x)]
 * @param x - function input
 */ 
double bessel(double x) {
    /*
    auto func = [=](const double val) { return cos(val - x * sin(val)); };
    double sum = 0.0;

    double step = PI / 1000;
    double last_val = func(0);
    for (uint32_t i = 1; i <= 1000; i++) {
        double new_val = func(i*step);
        sum += ((last_val + new_val) / 2) * step;
        last_val = new_val;
    }

    return (1/PI) * sum;
    */
   return boost::math::cyl_bessel_j(1, x);
}

/**
 * Calculates the intensity at a certain point given x and y
 * @param x - x coordinate on the aperture
 * @param y - y coordinate on the aperture
 */ 

double airy_func(double x, double y) {
    double theta = sqrt(pow(x,2) + pow(y,2)) / F;
    double b_in = ((2 * PI) / wavelength) *a*theta;
    if (b_in != 0) {
        double J1 = bessel(b_in);
        double base = ((2*J1)/b_in);
        return max_intensity * pow(base, 2);
    }
    // since (2J1(x)/x)^2 approaches 1 as x approaches 0
    return max_intensity;
}

/**
 * Represents the real portion of the integrand of the fourier transform of an airy function
 * @param x - x coordinate on the aperture
 * @param y - y coordinate on the aperture
 */ 

double real_airy(double x, double y) {
    const double airy = airy_func(x,y);
    return airy * cos(2 * PI * (u * x + v * y));
}

/**
 * Represents the complex portion of the integrand of the fourier transform of an airy function
 * @param x - x coordinate on the aperture
 * @param y - y coordinate on the aperture
 */ 
double complex_airy(double x, double y) {
    const double airy = airy_func(x,y);
    return -1 * airy * sin(2 * PI * (u * x + v * y));
}
  
/**
 * Numerically calculates a multiple integral over a square of two variables
 * @param function - two variable function of integration
 * @param lower - the lower x and y bound
 * @param upper - the upper x and y bound
 */ 
double double_integral(auto& function, double lower, double upper) { 
    double answer = 0.0;
    
    double step = (upper - lower) /resolution;
    std::vector<double> ax;
    ax.reserve(resolution + 1);

    for (int i = 0; i <= resolution; ++i) { 
        ax[i] = 0; 
        for (int j = 0; j <= resolution; ++j) { 
            if (j == 0 || j == resolution) 
                ax[i] += function(lower + i*step, lower + j*step); 
            else if (j % 2 == 0) 
                ax[i] += 2*function(lower + i*step, lower + j*step); 
            else
                ax[i] += 4 * function(lower + i*step, lower + j*step); 
        } 
        ax[i] *= (step / 3); 
    } 
  
    for (int i = 0; i <= resolution; ++i) { 
        if (i == 0 || i ==resolution) 
            answer += ax[i]; 
        else if (i % 2 == 0) 
            answer += 2 * ax[i]; 
        else
            answer += 4 * ax[i]; 
    } 
    answer *= (step / 3); 
  
    return answer; 
} 

/**
 * Returns the result of the numeric foureir transform of an airy function
 * @param lower - the lower bound of the fourier integral
 * @param upper - the upper bound of the fourier integral
 */ 
std::complex<double> airy_fourier_numeric(double lower, double upper) {
    double real_integral = double_integral(real_airy, lower, upper);
    double complex_integral = double_integral(complex_airy, lower, upper);
    std::complex<double> ans(real_integral, complex_integral);
    return ans;
}

double simpsons(auto& func,int lower, int upper, uint32_t steps) {
    double dx = (upper - lower) / steps;
    double sum = func(lower);
    for (uint32_t i = lower + 1; i < steps; i++) {
        if (i % 2 == 0) 
            sum += 2 * func(i * dx);
        else
            sum += 4 * func(i * dx);
    }
    sum += func(upper);

    return (dx / 3) * sum;
}
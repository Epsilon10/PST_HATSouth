#include "../include/psf.hh"
#include <iostream>

int main() {
    std::cout << fourier_psf(0.6, 2.0,.13,.14, 0.1, 0.1) << std::endl;
}
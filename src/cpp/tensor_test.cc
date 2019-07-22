#include <iostream>
#include "../include/tensor.hh"

int main() {
    std::vector<double> x;
    Tensor<2, double> t({3,3}, x);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            t[i][j] = i + j;
            std::cout << "i: " << i << " j: " << j << " res: " << &t[i][j] << std::endl;
         }
    }
    std::cout << "------------------------------------" << std::endl;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            std::cout << "i: " << i << " j: " << j << " res: " << t[i][j] << std::endl;
         }
    }
}

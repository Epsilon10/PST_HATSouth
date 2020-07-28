#pragma once

#include <vector>
#include <utility>
#include <cstdint>
#include <type_traits>

template<typename T>
std::pair<std::vector<T>,std::vector<T>> meshgrid(std::vector<T> const& x, std::vector<T> const& y) {
    static_assert(std::is_arithmetic_v<T>, "Must be of numeric types.");

    std::vector<T> x_vals(x.size() * y.size());
    std::vector<T> y_vals(x.size() * y.size());

    for (std::size_t i = 0; i < x.size(); i++) {
        for (std::size_t j = 0; j < y.size(); j++) {
            x_vals.push_back(x[i]);
            y_vals.push_back(y[j]);        
        }
    }

    return std::make_pair(x_vals, y_vals);
}
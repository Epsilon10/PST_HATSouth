#include <vector>
#include <initializer_list>
#include <numeric>
#include <array>
#include <iostream>

/**
 * Container to hold an n dimensional set of values in contiguous memory (useful for things like fftw)
 * @tparam D - dimension
 * @tparam T - data type
 * I need to fix a bug that this doesn't work with nonsquare dimensions (it overrites memory)
 */

template <size_t D, typename T>
class TensorProxy {
    public:
    TensorProxy(std::vector<T>& v, std::array<size_t, D-1> dims, size_t order) :
        v(v), order(order), dims(dims) {}
    
    TensorProxy<D-1, T> operator[](size_t i) {
        std::array<size_t, D-2> new_dims_;
        std::copy(dims.begin() + 1, dims.end(), new_dims_.begin());
        return TensorProxy<D-1, T>(v, new_dims_, order + i*dims[0]);
    }
    private:
    std::vector<T>& v;
    size_t order;
    std::array<size_t, D-1> dims;
};

template <typename T>
class TensorProxy<1, T> {
    public:
    TensorProxy(std::vector<T>& v, std::array<size_t, 0> dims, size_t order) :
        v(v), order(order) {}
     T& operator[](size_t i) {
        return v[order + i];
    }

    private:
    std::vector<T>& v;
    size_t order;
};

template <size_t D, typename T>
class Tensor {
    public:
    Tensor(std::initializer_list<size_t> a_dims, std::vector<T>& vec) : data(vec) {
        std::copy(a_dims.begin() + 1, a_dims.end(), dims.begin());
        std::cout << dims[0] << std::endl;
        data.resize(std::accumulate(a_dims.begin(), a_dims.end(), 1, std::multiplies<size_t>()));
    }
    TensorProxy<D-1, T> operator[](size_t i) {
        std::array<size_t, D-2> new_dims_;
        std::copy(dims.begin() + 1, dims.end(), new_dims_.begin());
        return TensorProxy<D-1, T>(data, new_dims_, i*dims[0]);
    }

    T* ptr() { return data.data(); }

    private:
    std::array<size_t, D-1> dims;
    std::vector<T> data;
};
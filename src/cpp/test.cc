#include <iostream>
#include <vector>
#include <utility>

int main() {
    std::vector<int> v(2);
    v.push_back(1);
    v.push_back(2);
    
    int x = 4;
    int const& y = x;
    int z = 3;
    std::pair<double,double> p(1.0, 3.0);
    auto* ptr = &p;
    std::cout << sizeof(ptr) << std::endl;
}


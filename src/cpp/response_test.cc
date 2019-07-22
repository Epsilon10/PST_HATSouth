#include "../include/response.hh"

int main() {
    Response response(0.6, 2.0, std::make_pair(0.1, 0.1), std::make_pair(300, 300), std::make_pair(5000, 5000));
    response.get_response(0.3,0.3);
    response.test_fourier();
}
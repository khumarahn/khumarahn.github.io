//#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>

#include <boost/math/tools/roots.hpp>
#include <tuple>

#include "cheb.h"

typedef long double real_t;
const int PREC = 64;
const int real_digits = std::numeric_limits<real_t>::digits;

const int CHEB_N = 128;

// find a root of a function given by Cheb with derivative using Newton-Raphson
// in interval [a,b]
real_t cheb_root(const cheb_ns::Cheb<real_t> &f, const cheb_ns::Cheb<real_t> &fp,
        const real_t &low, const real_t &high, const real_t &value) {
    auto FF = [&f, &fp, &value] (real_t y) {
        double fy = f.value(y) - value;
        double dfy = fp.value(y);
        return std::make_tuple(fy, dfy);
    };

    real_t guess = (low + high) / 2;
    real_t root = boost::math::tools::newton_raphson_iterate(
            FF, guess, low, high, real_digits);
    return root;
}

int main() {
    // construct a Chebyshev approximation and try to invert it!
    
    using std::cout;

    auto f = [] (real_t x) -> real_t {
        return x + x*x;
    };

    cheb_ns::Cheb<real_t> f_cheb(f, 0.0, 1.0, CHEB_N);
    auto fp_cheb = f_cheb.derivative();

    real_t val = 0.5;
    real_t root = cheb_root(f_cheb, fp_cheb, 0.0, 1.0, val);

    cout << "root: " << root
        << ", err: " << f(root) - val
        << "\n";

    for (real_t y = 0.25; y <= 0.75; y += 0.125) {
        cout << "gamma: " << f_cheb.value(y)
            << " vs " << f(y)
            << " err " << f_cheb.value(y) - f(y)
            << "\n";
    }

    return 0;
}

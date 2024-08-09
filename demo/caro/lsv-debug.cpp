#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>

#include "lsv.h"

// CHOOSE: boost multiprecision or (long) double;
// !! we need a patched version of Eigen for anything more precise than double
// !! because of https://gitlab.com/libeigen/eigen/-/issues/2841
#if 1
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/eigen.hpp>
namespace bmp = boost::multiprecision;
typedef bmp::number<bmp::mpfr_float_backend<64>> real_t; // in decimal digits
const int PREC = 180; // in bits
#else
typedef long double real_t;
const int PREC = 64;
#endif

#include "lsv-common.cpp"

int main() {
    LSV lsv;

    using std::cout;

    auto h = [&lsv] (real_t x) -> real_t {
        return lsv.h(x);
    };

    for (real_t gamma = 0.25; gamma <= 4.0; gamma += 0.25) {
        lsv.set_gamma(gamma);
        cout << "gamma: " << lsv.gamma()
            << "\n"
            << "h(1/2) / h(1) sanity check: " << h(0.5) / h(1) - (gamma + 2) / 2 << "\n"
            << "leading eigenvalue sanity check: " << abs(lsv.R_evalues(0) - complex_t(1)) << "\n"
            << "h(1/128) / h(1): " << (lsv.h_full(1./128).array() / lsv.h_full(1.).array()).transpose()  << "\n"
            << "\n";
    }

    return 0;
}

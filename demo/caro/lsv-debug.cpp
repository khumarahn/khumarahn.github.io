#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>

#include "lsv.h"

// CHOOSE: boost multiprecision or (long) double;
// !! we need a patched version of Eigen for anything more precise than double
// !! because of https://gitlab.com/libeigen/eigen/-/issues/2841
#if 0
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
    using std::abs;

    auto h = [&lsv] (real_t x) -> real_t {
        return lsv.h(x);
    };

    for (real_t gamma = 0.25; gamma <= 4.0; gamma += 0.25) {
        lsv.set_gamma(gamma);
        cout << "gamma: " << lsv.gamma()
            << "\n"
            << "h(1/2) / h(1) sanity check: " << h(0.5) / h(1) - (gamma + 2) / 2 << "\n"
            << "leading eigenvalue sanity check: " << abs(lsv.R_evalues(0) - complex_t(1)) << "\n"
            << "h(1/128) / h(1): " << (lsv.h_full(1./128).array() / lsv.h_full(1.).array()).transpose()  << "\n";

        // --- ABEL FUNCTION TESTS ---
        real_t x = 0.3; // Test point on the left branch

        // 1. Calculate f(x), f'(x), and f''(x)
        LSV::Vector2r fx_vec = lsv.left(x);
        real_t fx = fx_vec(0);
        real_t f_p = fx_vec(1);
        real_t f_pp = 2 * gamma * (gamma + 1) * pow(2 * x, gamma - 1);

        // 2. Evaluate A, A', A'' at x and f(x)
        real_t Ax = lsv.full_abel(x);
        real_t Ax_p = lsv.full_abel_p(x);
        real_t Ax_pp = lsv.full_abel_pp(x);

        real_t Afx = lsv.full_abel(fx);
        real_t Afx_p = lsv.full_abel_p(fx);
        real_t Afx_pp = lsv.full_abel_pp(fx);

        // 3. Functional Equation Checks
        cout << "  Abel Functional Equation Errors:\n"
             << "    A(f(x)) - (A(x) - 1): " << abs(Afx - (Ax - 1.0)) << "\n"
             << "    A'(f(x))*f'(x) - A'(x): " << abs(Afx_p * f_p - Ax_p) << "\n"
             << "    A''(f(x))*(f'(x))^2 + A'(f(x))*f''(x) - A''(x): "
             << abs(Afx_pp * f_p * f_p + Afx_p * f_pp - Ax_pp) << "\n";

        // 4. Finite Difference Checks (for A' and A'')
        // Note: eps = 1e-5 is reasonable for 64-bit float/long double finite differences.
        // Decrease it to ~1e-20 if you enable boost::multiprecision with 180 bits.
        real_t eps = 1e-5;
        real_t Ax_plus = lsv.full_abel(x + eps);
        real_t Ax_minus = lsv.full_abel(x - eps);

        real_t fd_p = (Ax_plus - Ax_minus) / (2 * eps);
        real_t fd_pp = (Ax_plus - 2 * Ax + Ax_minus) / (eps * eps);

        cout << "  Finite Difference Errors (eps=" << eps << "):\n"
             << "    A'(x) vs FD: " << abs(Ax_p - fd_p) << "\n"
             << "    A''(x) vs FD: " << abs(Ax_pp - fd_pp) << "\n\n";
    }

    return 0;
}

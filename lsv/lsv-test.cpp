//#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>

//#include "web.h"

#include "lsv.h"

const int PREC = 32;

#include "lsv-common.cpp"

using std::cout;

using interval_t = LSV::interval_t;

int main() {

    LSV lsv;

    auto UPPER = [&lsv](const interval_t &x) { return lsv.UPPER(x); };
    //auto LOWER = [&lsv](const interval_t &x) { return lsv.LOWER(x); };
    const interval_t ONE = lsv.ONE, HALF = lsv.HALF;

    // 0
    lsv.set_alpha(0.875);
    // 1
    lsv.compute_L();
    // 2
    lsv.compute_h_meta();
    // 3
    lsv.compute_h_cheb();
    // 4
    lsv.compute_F();

    // a sanity check for h
    {
        cout << "\nSanity check for h:\n";
        interval_t x = ONE / 1024,
                   hx = lsv.H(x);
        while (UPPER(x) < HALF) {
            auto Y = lsv.left(x);
            interval_t hy = lsv.H((Y(0) + 1) / 2) / 2 + hx / Y(1);
            x = Y(0);
            hx = hy;
        }
        interval_t Hx = lsv.H(x);
        interval_t err = bmp::width(interval_t(Hx - hx));
        cout << "  x: " << x << ", h(x) computed from infinite sum: " << hx
            << "  h(x) from Chebyshev: " << lsv.H(x) << "\n"
            << "  diff: " << err << "\n";
        verify(err < 1e-4);
    }

    // 5
    lsv.compute_derivative_signs_right();
    // 6
    lsv.compute_derivative_bounds();

    std::cout << lsv.oracle("min_hp_h_prime") << "\n"
        << lsv.oracle("max_hpp_h_prime") << "\n"
        << lsv.oracle("alpha") << "\n"
        << lsv.oracle("gamma") << "\n";

    return 0;
}

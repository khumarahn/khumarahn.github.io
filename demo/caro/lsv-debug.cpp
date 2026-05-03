//#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>

const int PREC = 48;   // in bits

// headers for multiprecision and interval arithmetics

#include <Eigen/Dense>

#include "lsv.h"

namespace bmp = boost::multiprecision;

int main() {
    using LSV = lsv_ns::LSV<PREC>;

    // real and interval types
    using real_t = LSV::real_t;
    using interval_t = LSV::interval_t;
    //using complex_interval_t = LSV::complex_interval_t;

    // Matrix and vector types in interval arithmetic
    using MatrixXi = LSV::MatrixXi;
    //using MatrixXix = LSV::MatrixXix;
    //using MatrixXr = LSV::MatrixXr;
    using VectorXi = LSV::VectorXi;
    //using Vector2ci = LSV::Vector2ci;

    const real_t real_eps = std::numeric_limits<real_t>::epsilon();

    using std::cout;

    // Chebyshev interpolation class in interval arithmetic
    using Cheb_i = cheb_ns::Cheb<interval_t>;

    auto uncertainty = [] (const MatrixXi &M) {
        real_t err = 0;
        for (const auto &x : M.reshaped())
            err += bmp::width(x);
        return err;
    };

    for (real_t gamma = real_t(8) / 7; gamma <= 5; gamma += 50) {
        //cout << "gamma: " << gamma << "\n";
        LSV lsv(gamma);

        // Look at the coefficients of the Abel function
        {
            real_t err = 0;

            VectorXi abel_coef = lsv.abel_coef();
            err = 0;
            for (int i=0; i<abel_coef.size(); i++) {
                err += bmp::width(abel_coef(i));
            }
            //cout << "L1 error in Abel coefficients: " << err << "\n";

            //cout << "Abel coeff: " << abel_coef.transpose() << "\n";
        }

        auto h_meta = lsv.h_meta();
        // Retrieve the transfer operator in interval form,
        // as an N x N matrix acting on Chebyshev polynomials on [0.5,1]
        // with the first coeff doubled as in Numerical Recipes
        MatrixXi L = h_meta.L;
        //cout << "Transfer operator matrix computed...\n";
        //cout << L.topLeftCorner(5,5) << "\n...\n";
        //const int N = L.cols();
        //cout << "N: " << N << "\n";

        auto L_uncertainty = uncertainty(L);
        //cout << "Transfer operator retrieved, L1 uncertainty: " << L_uncertainty << " ...\n";

        interval_t a = 0.5, b = 1.0;

        // invariant density h
        Cheb_i h = h_meta.h;

        VectorXi hv = h.coef();

        //cout << "Error in Lh - h: " << (L * hv - hv).norm() << "\n";
        //cout << "First coeff of h: \n" << hv.head(5).transpose()
        //    << "\n";
        //cout << "L1 uncertainty of coefficients: " << uncertainty(hv) << "\n";

        cout << "\nrho_A: " << h_meta.rho_A
            << "\nERROR: " << h_meta.err << "\n";
    }

    return 0;
}

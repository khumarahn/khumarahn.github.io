//#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>

const int PREC = 32;   // in bits

// headers for multiprecision and interval arithmetics

#include <Eigen/Dense>

#include "lsv.h"

namespace bmp = boost::multiprecision;



int main() {
    using LSV = lsv_ns::LSV<PREC>;
    LSV lsv;

    // real and interval types
    using real_t = LSV::real_t;
    using interval_t = LSV::interval_t;
    using complex_interval_t = LSV::complex_interval_t;

    // Matrix and vector types in interval arithmetic
    using MatrixXi = LSV::MatrixXi;
    using MatrixXix = LSV::MatrixXix;
    using MatrixXr = LSV::MatrixXr;
    using VectorXi = LSV::VectorXi;
    using Vector2ci = LSV::Vector2ci;

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
    auto uncertainty_x = [] (const MatrixXix &M) {
        real_t err = 0;
        for (const auto &x : M.reshaped())
            err += real_t(bmp::width(x));
        return err;
    };

    if (0) {
        auto f = [] (complex_interval_t x) {
            return Vector2ci(x * x - complex_interval_t(0,2), complex_interval_t(2) * x);
        };

        complex_interval_t guess(interval_t(0.75, 1.25), interval_t(0.75,1.25));

            cout << "xxx: " << lsv.complex_krawczyk(f, guess)
            << "\n";
            return 0;
    }

    for (real_t gamma = 0.75; gamma <= 0.75; gamma += 0.25) {
        cout << "gamma: " << gamma << "\n";
        lsv.set_gamma(gamma);

        cout << "abel(1): " << lsv.abel(interval_t(1)) << "\n";
        // Look at the coefficients of the Abel function
        {
            real_t err = 0;

            MatrixXix abel_mat = lsv.abel_matrix();
            cout << "L1 error in Abel matrix: " << uncertainty_x(abel_mat) << "\n";

            VectorXi abel_coef = lsv.abel_coef();
            err = 0;
            for (int i=0; i<abel_coef.size(); i++) {
                err += bmp::width(abel_coef(i));
            }
            cout << "L1 error in Abel coefficients: " << err << "\n";

            cout << "Abel coeff: " << abel_coef.transpose() << "\n";
        }

        // Retrieve the transfer operator in interval form,
        // as an RN x RN matrix acting on Chebyshev polynomials on [0.5,1]
        // with the first coeff doubled as in Numerical Recipes
        MatrixXi R = lsv.Lind();
        cout << "Transfer operator matrix computed...\n";
        const int RN = R.cols();
        cout << "RN: " << RN << "\n";

        // Make it a bit uncertain
        cout << "Transfer operator retrieved, L1 uncertainty: " << uncertainty(R) << " ...\n";

        interval_t a = 0.5, b = 1.0;

        Cheb_i hh;
        // computation of invariant density without eigen decomposition
        {
            // On the vector space W of first RN Chebyshev polynomials,
            // we compute iota : W \to R, w \mapsto \int_{0.5}^1 w
            // and u \in W so that iota(u) = 1
            VectorXi iota(RN);
            VectorXi u(RN);

            // inefficiently compute values of integrals on [1/2,1] of the basis Chebyshev polynomials
            for (int i=0; i<RN; i++) {
                // higher dimension because when taking an integral, the coefficients shift
                VectorXi co(RN+1);
                co(i) = interval_t(1);

                Cheb_i cheb(co, a, b);
                Cheb_i ii = cheb.integral();

                iota(i) = ii.value(b) - ii.value(a);
            }
            cout << "Integrals of first basis vectors: \n" << iota.head(5).transpose() << "\n";

            u(0) = 1 / iota(0);

            // Now the invariant density in Chebyshev basis
            // is hh = (I - R + u iota)^{-1} u
            MatrixXi S = MatrixXi::Identity(RN,RN) - R + u * iota.transpose();
            MatrixXr Ss(RN, RN);
            for (int i=0; i<RN; i++)
                for (int j=0; j<RN; j++)
                    Ss(i,j) = bmp::median(S(i,j));
            MatrixXr Ssi = Ss.inverse();
            MatrixXi Si(RN,RN);
            for (int i=0; i<RN; i++)
                for (int j=0; j<RN; j++)
                    Si(i,j) = Ssi(i,j);

            //VectorXi hv = (Si * S).householderQr().solve(Si * u);
            VectorXi hv = interval_root_ns::linear_krawczyk(S, u);

            cout << "Error in computation of hv: " << (S * hv - u).norm() << "\n";

            hh = Cheb_i(hv, a, b);

            cout << "First coeff of h: \n" << hv.head(5).transpose()
                << "\n";
            cout << "L1 uncertainty of coefficients: " << uncertainty(hv) << "\n";
        }

    }

    return 0;
}

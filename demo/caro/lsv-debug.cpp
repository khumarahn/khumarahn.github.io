//#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>

// precision
const int DIGITS = 16;                   // in decimal digits
const int PREC = (DIGITS * 332) / 100;   // in bits

// headers for multiprecision and interval arithmetics
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpfi.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/eigen.hpp>

namespace bmp = boost::multiprecision;

// real and interval types
using real_t = bmp::number<bmp::mpfr_float_backend<DIGITS>>;
using interval_t = bmp::number<bmp::mpfi_float_backend<DIGITS>>;

const real_t real_epsilon = std::numeric_limits<real_t>::epsilon();

#include "lsv.h"
#include "lsv-common.cpp"

int main() {
    LSV lsv;

    using std::cout;

    // Matrix and vector types in interval arithmetic
    using MatrixXri = Eigen::Matrix<interval_t, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXri = Eigen::Matrix<interval_t, Eigen::Dynamic, 1>;

    // Chebyshev interpolation class in interval arithmetic
    using Cheb_i = cheb_ns::Cheb<interval_t>;


    for (real_t gamma = 0.75; gamma <= 0.75; gamma += 0.25) {
        cout << "gamma: " << gamma << "\n";
        lsv.set_gamma(gamma);
        cout << "h(1/2) / h(1) sanity check: " << lsv.h(0.5) / lsv.h(1) - (gamma + 2) / 2 << "\n"
            << "leading eigenvalue sanity check: " << abs(lsv.R_evalues(0) - complex_t(1)) << "\n"
            << "\n";

        cout << "LSV basic things computed...\n";

        // Retrieve the transfer operator in interval form,
        // as an RN x RN matrix acting on Chebyshev polynomials on [0.5,1]
        // with the first coeff doubled as in Numerical Recipes
        //
        // Make it a bit uncertain
        const int RN = lsv.R_cols();
        cout << "RN: " << RN << "\n";
        MatrixXri R(RN,RN);
        for (int i=0; i<RN; i++) {
            for (int j=0; j<RN; j++) {
                R(i,j) = interval_t(lsv.R_coef(i,j));
                // add uncertainty
                R(i,j) *= interval_t(1 - 1024 * real_epsilon, 1 + 1024 * real_epsilon);
            }
        }
        cout << "Transfer operator retrieved...\n";

        interval_t a = 0.5, b = 1.0;

        Cheb_i hh;
        // computation of invariant density without eigen decomposition
        {
            // On the vector space W of first RN Chebyshev polynomials,
            // we compute iota : W \to R, w \mapsto \int_{0.5}^1 w
            // and u \in W so that iota(u) = 1
            VectorXri iota(RN);
            VectorXri u(RN);

            // inefficiently compute values of integrals on [1/2,1] of the basis Chebyshev polynomials
            for (int i=0; i<RN; i++) {
                // higher dimension because when taking an integral, the coefficients shift
                VectorXri co(RN+1);
                co(i) = interval_t(1);

                Cheb_i cheb(co, a, b);
                Cheb_i ii = cheb.integral();

                iota(i) = ii.value(b) - ii.value(a);
            }
            cout << "Integrals of first basis vectors: \n" << iota.head(5).transpose() << "\n";

            u(0) = 1 / iota(0);

            // Now the invariant density in Chebyshev basis
            // is hh = (I - R + u iota)^{-1} u
            MatrixXri S = MatrixXri::Identity(RN,RN) - R + u * iota.transpose();
            VectorXri hv = S.householderQr().solve(u);

            cout << "Error in computation of hv: " << (S * hv - u).norm() << "\n";

            hh = Cheb_i(hv, a, b);

            cout << "First coeff of h: \n" << hv.head(5).transpose()
                << "\n";
            real_t err = 0, errr;
            for (int i=0; i<RN; i++) {
                interval_t rr = abs(hv(i) - lsv.h_coef(i));
                err += bmp::upper(rr);
                errr += bmp::upper(hv(i)) - bmp::lower(hv(i));
            }
            cout << "L1 error compated to computation using eigen decomposition: " << err << "\n"
                << "L1 uncertainty of coefficients: " << errr << "\n";
        }

    }

    return 0;
}

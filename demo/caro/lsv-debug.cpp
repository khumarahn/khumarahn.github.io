//#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>

const int PREC = 96;   // in bits

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

    const real_t real_epsilon = std::numeric_limits<real_t>::epsilon();

    using std::cout;

    // Matrix and vector types in interval arithmetic
    using MatrixXri = Eigen::Matrix<interval_t, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXri = Eigen::Matrix<interval_t, Eigen::Dynamic, 1>;

    // Chebyshev interpolation class in interval arithmetic
    using Cheb_i = cheb_ns::Cheb<interval_t>;


    for (real_t gamma = 0.75; gamma <= 0.75; gamma += 0.25) {
        cout << "gamma: " << gamma << "\n";
        lsv.set_gamma(gamma);

        // Retrieve the transfer operator in interval form,
        // as an RN x RN matrix acting on Chebyshev polynomials on [0.5,1]
        // with the first coeff doubled as in Numerical Recipes
        auto Lind = lsv.Lind();
        cout << "Transfer operator matrix computed...\n";
        const int RN = Lind.cols();
        cout << "RN: " << RN << "\n";

        // Make it a bit uncertain
        MatrixXri R(RN,RN);
        real_t errrr = 0;
        for (int i=0; i<RN; i++) {
            for (int j=0; j<RN; j++) {
                R(i,j) = Lind(i,j);
                errrr += bmp::upper(R(i,j)) - bmp::lower(R(i,j));
                // add uncertainty
            }
        }
        cout << "Transfer operator retrieved, L1 uncertainty: " << errrr << " ...\n";

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
            //VectorXri hv = S.householderQr().solve(u);
            VectorXri hv = S.partialPivLu().solve(u);

            cout << "Error in computation of hv: " << (S * hv - u).norm() << "\n";

            hh = Cheb_i(hv, a, b);

            cout << "First coeff of h: \n" << hv.head(5).transpose()
                << "\n";
            real_t err = 0;
            for (int i=0; i<RN; i++) {
                err += bmp::upper(hv(i)) - bmp::lower(hv(i));
            }
            cout << "L1 uncertainty of coefficients: " << err << "\n";
        }

    }

    return 0;
}

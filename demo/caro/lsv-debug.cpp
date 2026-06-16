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
    using cheb_t = LSV::interval_cheb_t;

    auto uncertainty = [] (const MatrixXi &M) {
        real_t err = 0;
        for (const auto &x : M.reshaped())
            err += bmp::width(x);
        return err;
    };

    for (interval_t gamma = interval_t(8) / 7; gamma <= 5; gamma += 50) {
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
        cheb_t h = h_meta.h;

        VectorXi hv = h.coef();

        //cout << "Error in Lh - h: " << (L * hv - hv).norm() << "\n";
        //cout << "First coeff of h: \n" << hv.head(5).transpose()
        //    << "\n";
        //cout << "L1 uncertainty of coefficients: " << uncertainty(hv) << "\n";

        auto H = [&lsv, &h_meta] (const interval_t &x) {
            interval_t r;
            // sum without an error
            interval_t S = lsv.cheb_sum(x).transpose() * h_meta.h.coef();
            S /= 2;
            // an error
            interval_t E = lsv.eps_sum(x, h_meta.err);

            r = S + E;
            return r;
        };

        // try computing h
        for (interval_t x = 0.75; x > 1.0 / 1024; x /= 16) {
            //interval_t x = interval_t(1.0 / 1024, 1.0 / 1024 + 1e-6);
            interval_t hx = h.value(x) + interval_t(-h_meta.err, h_meta.err),
                       Hx = H(x),
                       Hhx = Hx - hx;

            cout << "\n"
                << "h(" << x << "):" << hx << ", width: " << bmp::width(hx)
                << "\n"
                << "H(" << x << "):" << Hx << ", width: " << bmp::width(Hx)
                << "\n"
                << "diff: " << Hhx << ", width: " << bmp::width(Hhx)
                << "\n";

        }

        cout << "\nrho_A: " << h_meta.rho_A
            << "\nERROR: " << h_meta.err << "\n";

        {
            cout << "\nComputing F_1, F_2:\n";
            const auto abel = lsv.abel_meta();
            const auto &h = h_meta.h;

            interval_t h_norm_A = h.ellipse_norm(h_meta.rho_A) + h_meta.err;
            h_norm_A = bmp::upper(h_norm_A);

            // h(1/2)
            interval_t h_half = h.value(interval_t(1)/2)
                + interval_t(-h_meta.err, h_meta.err);
            cout << "  h(1/2): " << h_half
                << ", width: " << bmp::width(h_half) << "\n";

            interval_t delta_star = abel.dist_P1_to_A;
            cout << "  delta_star: " << delta_star << "\n";

            interval_t h_prime_hash;
            {
                // before bounding the derivative,
                // subtract the constant part from h
                auto h0_coef = h.coef();
                h0_coef(0) = 0;
                cheb_t h0(h0_coef, h.a(), h.b());

                h_prime_hash = (h0.ellipse_norm(h_meta.rho_A) + h_meta.err)
                    / delta_star;
                h_prime_hash = bmp::upper(h_prime_hash);
            }
            cout << "  bound on sup norm of h in A: " << h_norm_A
                << ", and h' delta_star-inside A: " << h_prime_hash << "\n";

            interval_t abel_am1 = abel.coef(0);
            interval_t varkappa1 = abel.varkappa0
                * 2 * (abel_am1 * abel.r1 + abel.nu) * abel.nu
                / (pow(abel.nu, 2) + pow(abel_am1 * abel.r1 + abel.nu, 2));
            varkappa1 = bmp::lower(varkappa1);
            cout << "  varkappa1: " << varkappa1 << "\n";

            interval_t q = abel.C1 / abel_am1;
            interval_t R = 1 - pow(q + varkappa1 * (1 + q), 2);
            assert(R > 0);
            R = sqrt(R);

            interval_t r_star_1 = abel.r1 + abel.nu / abel_am1,
                       r_star_2 = varkappa1 * abel.nu
                           / ((abel_am1 - abel.C1) * (1 - R)),
                       r_star = max(r_star_1, r_star_2);
            interval_t x_star = pow( interval_t(2), - 1 - gamma / 2)
                * pow(r_star, - 1 / gamma);

            cout << "  R: " << R
                << ", r_star: " << r_star
                << ", x_star: " << x_star << "\n";

            interval_t C_psi = h_norm_A * (abel_am1 + abel.C1)
                / (abel_am1 - abel.C1)
                * pow(R, - (1 + 1 / gamma));
            interval_t R_L = 6 * C_psi
                / (pow(2 * lsv.pi_ * varkappa1, 3) * pow(abel.nu, 2));
            R_L = bmp::upper(R_L);

            interval_t C_Delta =
                abs(gamma - 1) / (gamma * r_star * (abel_am1 - abel.C1))
                + (2 * abel_am1 + abel.C2) / (r_star * pow(abel_am1 - abel.C1, 2));
            interval_t Delta_star =
                (gamma / 2 * abel.r1 * abel.C1 + interval_t(1) / 4
                 + C_Delta / 24) * h_norm_A
                + pow(r_star, - 1 - 1 / gamma) / (48 * gamma
                        * (abel_am1 - abel.C1)) * h_prime_hash
                + R_L / 2;
            Delta_star = bmp::upper(Delta_star);

            cout << "  R_L: " << R_L << ", Delta_star: " << Delta_star << "\n";
            interval_t F1, F2;
            {
                interval_t theta_star = 1 / sin(min(
                            lsv.pi_ / 2, lsv.pi_ / (4 * gamma)
                            ));
                interval_t D0 = gamma * abel_am1 * h_prime_hash / 8,
                           D1 = D0 * (gamma + 3),
                           D2 = D0 * (pow(gamma, 2) + 7 * gamma + 8);
                interval_t zh = Delta_star * pow(x_star, gamma);
                interval_t delta0 = D0 * x_star + zh,
                           delta1 = D1 * x_star + 2 * zh * theta_star,
                           delta2 = D2 * x_star + 8 * zh * pow(theta_star, 2);
                interval_t M0 = gamma * abel_am1 * h_half / 2;

                F1 = (gamma * M0 - delta1) / (M0 + delta0);
                F1 = bmp::lower(F1);
                F2 = (gamma * (gamma + 1) * M0 - delta2) / (M0 + delta0);
                F2 = bmp::lower(F2);
            }

            cout << "  F1: " << F1 << ", F2: " << F2 << "\n";
            interval_t F1_x = F1 / x_star,
                       F2_x = F2 / pow(x_star, 2);
            F1_x = bmp::lower(F1_x);
            F2_x = bmp::lower(F2_x);

            cout << "\nAlmost Theorem. Let gamma = " << gamma
                << ". For all 0 < x <= " << bmp::lower(x_star) << ",\n"
                << "  h'(x) / h(x) <= - " << F1 << " / x <= - " << F1_x << ", and\n"
                << "  h''(x) / h(x) >= " << F2 << " / x^2 >= " << F2_x << ".\n"
                << "(would be a theorem if the numbers were rounded correctly)\n";

            // sanity check
            {
                interval_t hd = x_star * 0.001,
                           Hm = H(x_star - hd),
                           H0 = H(x_star),
                           Hp = H(x_star + hd);
                interval_t hp = (Hp - Hm) / (2 * hd),
                           hpp = (Hm - 2 * H0 + Hp) / (hd * hd);
                cout << "\n\nSanity check: approximately, at x = x_star,\n"
                    << "  h'(x) / h(x) * x : "
                    << hp / H0 * x_star << ",\n"
                    << "  h''(x) / h(x) * x^2 : "
                    << hpp / H0 * pow(x_star, 2) << "\n\n";
            }
        } // computing F1, F2
    }

    return 0;
}

//#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>

const int PREC = 48;   // in bits

// headers for multiprecision and interval arithmetics

#include <Eigen/Dense>

#include "lsv.h"

namespace bmp = boost::multiprecision;

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

// Bounds f'(x) on [x1, x2]
// m2 is the bound on |f''| over the interval
interval_t derivative_range(const interval_t& x1, const interval_t& x2,
        const interval_t& y1, const interval_t& y2,
        const interval_t& m2);

// Bounds f''(x) on [x1, x3] (assuming x1 < x2 < x3)
// m3 is the bound on |f'''| over the interval
interval_t second_derivative_range(const interval_t& x1, const interval_t& x2, const interval_t& x3,
        const interval_t& y1, const interval_t& y2, const interval_t& y3,
        const interval_t& m3);

int main() {
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

        const auto h_meta = lsv.h_meta();
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

        const interval_t a = 0.5, b = 1.0,
              rho_A = h_meta.rho_A;

        cout << "\nrho_A: " << h_meta.rho_A
            << "\nERROR: " << h_meta.err << "\n";

        // h(x) on (0,1] with error bounds
        auto H = [&lsv, &h_meta] (const interval_t &x) -> interval_t {
            if (bmp::lower(x) >= 0.5) {
                return h_meta.h.value(x) + interval_t(-1, 1) * h_meta.err;
            } else {
                interval_t r;
                // sum without an error
                interval_t S = lsv.cheb_sum(x).transpose() * h_meta.h.coef();
                // the error
                interval_t E = lsv.eps_sum(x, h_meta.err);

                r = (S + E) / 2;
                return r;
            }
        };
        // derivatives of h on [1/2,1]
        cheb_t hp_cheb = h_meta.h.derivative(),
               hpp_cheb = hp_cheb.derivative();
        interval_t dist_A = (rho_A + 1 / rho_A - 2) / 8,
                   hp_cheb_err = h_meta.err / dist_A,
                   hpp_cheb_err = 2 * h_meta.err / pow(dist_A, 2);
        hp_cheb_err  = bmp::upper(hp_cheb_err);
        hpp_cheb_err = bmp::upper(hpp_cheb_err);
        auto Hp = [&hp_cheb, &hp_cheb_err] (const interval_t &x) {
            assert(bmp::lower(x) >= 0.5 && bmp::upper(x) <= 1);
            return hp_cheb.value(x) + interval_t(-1,1) * hp_cheb_err;
        };
        auto Hpp = [&hpp_cheb, &hpp_cheb_err] (const interval_t &x) {
            assert(bmp::lower(x) >= 0.5 && bmp::upper(x) <= 1);
            return hpp_cheb.value(x) + interval_t(-1,1) * hpp_cheb_err;
        };

        // a sanity check for h
        {
            cout << "\nSanity check for h:\n";
            interval_t x = 0.001,
                       hx = H(x);
            while (bmp::upper(x) < 0.5) {
                auto Y = lsv.left(x);
                interval_t hy = H((Y(0) + 1) / 2) / 2 + hx / Y(1);
                x = Y(0);
                hx = hy;
            }
            interval_t Hx = H(x);
            cout << "  x: " << x << ", h(x) computed from infinite sum: " << hx
                << "  h(x) from Chebyshev: " << H(x)
                << ", diff: " << bmp::width(interval_t(Hx - hx))
                << "\n";
        }

        interval_t x_star,
                   h_norm_A;

        {
            cout << "\n\nComputing F_1, F_2:\n";
            const auto abel = lsv.abel_rough_meta();
            const auto &h = h_meta.h;

            h_norm_A = h.ellipse_norm(h_meta.rho_A) + h_meta.err;
            h_norm_A = bmp::upper(h_norm_A);

            // h(1/2)
            interval_t h_half = h.value(interval_t(1)/2)
                + interval_t(-h_meta.err, h_meta.err);
            cout << "  h(1/2): " << h_half
                << ", width: " << bmp::width(h_half) << "\n";

            interval_t delta_hash = abel.dist_P1_to_A;
            cout << "  delta_hash: " << delta_hash << "\n";

            interval_t h_prime_hash;
            {
                // before bounding the derivative,
                // subtract the constant part from h
                auto h0_coef = h.coef();
                h0_coef(0) = 0;
                cheb_t h0(h0_coef, h.a(), h.b());

                h_prime_hash = (h0.ellipse_norm(h_meta.rho_A) + h_meta.err)
                    / delta_hash;
                h_prime_hash = bmp::upper(h_prime_hash);
            }
            cout << "  bound on sup norm of h in A: " << h_norm_A
                << ", and h' delta_hash-inside A: " << h_prime_hash << "\n";

            // choose nu arbitrarily
            interval_t nu = 4;
            cout << "nu: " << nu << "\n";
            interval_t abel_am1 = abel.coef(0);
            interval_t varkappa1 = abel.varkappa0
                * 2 * (abel_am1 * abel.r1 + nu) * nu
                / (pow(nu, 2) + pow(abel_am1 * abel.r1 + nu, 2));
            varkappa1 = bmp::lower(varkappa1);
            cout << "  varkappa1: " << varkappa1 << "\n";

            interval_t q = abel.C1 / abel_am1;
            interval_t R = 1 - pow(q + varkappa1 * (1 + q), 2);
            assert(R > 0);
            R = sqrt(R);

            interval_t r_star_1 = abel.r1 + nu / abel_am1,
                       r_star_2 = varkappa1 * nu
                           / ((abel_am1 - abel.C1) * (1 - R)),
                       r_star = max(r_star_1, r_star_2);
            x_star = pow(interval_t(2), - 1 - 1 / (2 * gamma))
                * pow(r_star, - 1 / gamma);

            cout << "  R: " << R
                << ", r_star: " << r_star
                << ", x_star: " << x_star << "\n";

            interval_t C_psi = h_norm_A * (abel_am1 + abel.C1)
                / (abel_am1 - abel.C1)
                * pow(R, - (1 + 1 / gamma));
            interval_t R_L = 6 * C_psi
                / (pow(2 * lsv.pi_ * varkappa1, 3) * pow(nu, 2));
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
            assert(F1 > 0 && F2 > 0);

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
                           H_m = H(x_star - hd),
                           H0 = H(x_star),
                           H_p = H(x_star + hd);
                interval_t hp = (H_p - H_m) / (2 * hd),
                           hpp = (H_m - 2 * H0 + H_p) / (hd * hd);
                cout << "\nTheorem sanity check: approximately, at x = " << x_star << ",\n"
                    << "  h'(x) / h(x) * x: "
                    << hp / H0 * x_star << ",\n"
                    << "  h''(x) / h(x) * x^2: "
                    << hpp / H0 * pow(x_star, 2) << "\n\n";
            }
        } // computing F1, F2

        {   // verifying conditions away from zero

            // bound on |h^{(k)}| on [x1, x2]
            auto hp_max = [&gamma, &h_meta, &h_norm_A]
                (interval_t x1, interval_t x2, int k) -> interval_t {
                assert(0 < x1 && x1 <= x2 && x2 <= 1);

                interval_t rho_A = h_meta.rho_A,
                           theta_C = h_meta.theta_C;
                // distance from [1/2,1] to \partial A
                interval_t di = (rho_A + 1 / rho_A - 2) / 8;
                // radius around [x1, x2] where we control |h|
                interval_t dx = min(di, bmp::lower(x1) * sin(theta_C));

                interval_t C = pow(2, -gamma),
                           M = (gamma + 1) * pow(2, gamma),
                           sum_Jk = 1 + C * pow(x1 - dx, - gamma)
                               * exp(M * pow(x2 + dx, gamma));

                interval_t r = sum_Jk * h_norm_A;
                for (int j = 1; j <= k; j++)
                    r *= j / dx;

                return bmp::upper(r);
            };

            interval_t x2 = bmp::lower(x_star),
                       h2 = H(x2),
                       dx = 0.1 * h2 * pow(x2, 2) / hp_max(x2, x2, 2);
            assert(dx < x2 / 2);
            interval_t x1 = x2 - dx,
                       x3 = x2 + dx;
            x1 = bmp::lower(x1);
            x3 = bmp::upper(x3);
            interval_t h1 = H(x1),
                       h3 = H(x3);
            interval_t x13(x1, x3),
                       h13(bmp::lower(h3), bmp::upper(h1));
            interval_t K2 = hp_max(x1, x3, 2),
                       K3 = hp_max(x1, x3, 3);
            cout << "K2: " << K2 << ", K3: " << K3 << "\n";

            // estimate h'(t), h''(t) on [x1, x3]
            interval_t hp = derivative_range(x1, x3, h1, h3, K2);
            interval_t hpp = second_derivative_range(x1, x2, x3,
                    h1, h2, h3, K3);

            cout << "Derivatives at x: " << x2
                << ", x3 - x1: " << x3 - x1 << "\n"
                << " x h'(x1, x3) / h(x): " << x13 * hp / h13 << "\n"
                << " x^2 h''(x1, x3) / h(x): " << pow(x13, 2) * hpp / h13
                << "\n";
        }
    }

    return 0;
}

interval_t derivative_range(const interval_t& x1, const interval_t& x2,
        const interval_t& y1, const interval_t& y2,
        const interval_t& m2) {

    interval_t dx = x2 - x1;
    interval_t mean_slope = (y2 - y1) / dx;

    // The maximum deviation of f'(x) from f'(c) on [x1, x2] is m2 * dx
    interval_t max_dev = m2 * dx;
    max_dev = bmp::upper(max_dev);
    interval_t error_term(-max_dev, max_dev);

    return mean_slope + error_term;
}

interval_t second_derivative_range(const interval_t& x1, const interval_t& x2, const interval_t& x3,
        const interval_t& y1, const interval_t& y2, const interval_t& y3,
        const interval_t& m3) {

    interval_t dx12 = x2 - x1;
    interval_t dx23 = x3 - x2;
    interval_t dx13 = x3 - x1;

    interval_t slope12 = (y2 - y1) / dx12;
    interval_t slope23 = (y3 - y2) / dx23;

    // Using divided differences, f''(c) = 2 * f[x1, x2, x3]
    interval_t mean_second_deriv = interval_t(2) * (slope23 - slope12) / dx13;

    interval_t max_dev = m3 * dx13;
    max_dev = bmp::upper(max_dev);
    interval_t error_term(-max_dev, max_dev);

    return mean_second_deriv + error_term;
}

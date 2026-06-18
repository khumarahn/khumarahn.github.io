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
                << ", h(x) from Chebyshev: " << H(x)
                << ", diff: " << bmp::width(interval_t(Hx - hx))
                << "\n";
        }

        interval_t x_star,
                   h_norm_A;

        {
            cout << "\n\nComputing F_1, F_2:\n";
            const auto abel = lsv.abel_meta();
            const auto &h = h_meta.h;

            h_norm_A = h.ellipse_norm(h_meta.rho_A) + h_meta.err;
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
            x_star = pow( interval_t(2), - 1 - gamma / 2)
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
                           H_m = H(x_star - hd),
                           H0 = H(x_star),
                           H_p = H(x_star + hd);
                interval_t hp = (H_p - H_m) / (2 * hd),
                           hpp = (H_m - 2 * H0 + H_p) / (hd * hd);
                cout << "\nSanity check: approximately, at x = " << x_star << ",\n"
                    << "  h'(x) / h(x) * x: "
                    << hp / H0 * x_star << ",\n"
                    << "  h''(x) / h(x) * x^2: "
                    << hpp / H0 * pow(x_star, 2) << "\n\n";
            }
        } // computing F1, F2

        {   // verifying conditions away from zero

            // bounds on derivatives of h:
            // |h^{(k)}(x)| \leq C_k x^{-\gamma - k}
            // where C_k = \frac{k! C_0}{\eta^k}
            interval_t C0, C1, C2, C3, C4;
            {
                interval_t rho_A = h_meta.rho_A,
                           theta_C = h_meta.theta_C;
                interval_t eta = min(sin(theta_C), (rho_A + 1 / rho_A - 2) / 4);
                assert(0 < eta && eta < 1);
                C0 = h_norm_A / 2 * (1 + 2 * pow(2 * (1 - eta), - gamma)
                        * exp( (gamma + 1) * pow(2 * (1 + eta), gamma)));
                C1 = C0 / eta;
                C2 = C0 * 2 / pow(eta, 2);
                C3 = C0 * 6 / pow(eta, 3);
                C4 = C0 * 24 / pow(eta, 4);
            }
            interval_t sigma = 0.2;
            interval_t x = bmp::lower(x_star),
                       x_delta = pow(x, gamma + 1) * sigma,
                       x_m = x - x_delta,
                       x_p = x + x_delta;
            x_m = bmp::lower(x_m);
            x_p = bmp::upper(x_p);
            interval_t h_max_p1 = C1 * pow(x - x_delta, - gamma - 1),
                       h_max_p2 = C2 * pow(x - x_delta, - gamma - 2),
                       h_max_p3 = C3 * pow(x - x_delta, - gamma - 3),
                       h_max_p4 = C4 * pow(x - x_delta, - gamma - 4);

            // estimate h'(x) on [x_m, x_p]
            // using g(t) = h(t^{-1/\gamma}) / t
            // @gemini: create a code to put here
            // use the function H(x) to compute h(x)

            interval_t H_m = H(x_m);
            interval_t H_p = H(x_p);

            interval_t t_m = pow(x_m, -gamma);
            interval_t t_p = pow(x_p, -gamma);
            interval_t delta = t_m - t_p;

            interval_t g_m = pow(x_m, gamma) * H_m;
            interval_t g_p = pow(x_p, gamma) * H_p;
            interval_t S = (g_m - g_p) / delta;

            interval_t M2 = pow(x_p, 2 * gamma) * (2 * C0 + (3 * gamma + 1) / pow(gamma, 2) * C1 * x_p + C2 / pow(gamma, 2) * pow(x_p, 2));

            interval_t g_prime = S + interval_t(-1, 1) * delta * M2;

            cout << "  t_m: " << t_m << ", delta: " << delta
                << ", S: " << S << ", M2: " << M2
                << ", g_prime: " << g_prime << "\n";

            interval_t X(bmp::lower(x_m), bmp::upper(x_p));
            interval_t T(bmp::lower(t_p), bmp::upper(t_m));

            // h is monotone decreasing
            interval_t H_X(bmp::lower(H_p), bmp::upper(H_m));
            interval_t g_T = pow(X, gamma) * H_X;

            interval_t x_h_prime_over_h = -gamma * (interval_t(1) + T * g_prime / g_T);

            cout << "Derivative of h(x) the interval x in " << interval_t(x_m, x_p)
                << " (width " << bmp::width(interval_t(x_m, x_p)) << ")"
                << " is " << x_h_prime_over_h << "\n"
                ;


        }
    }

    return 0;
}

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

// Bounds f(x) on [x1, x2]
// m1 is the bound on |f'| over the interval
interval_t function_range(const interval_t& x1, const interval_t& x2,
        const interval_t& y1, const interval_t& y2,
        const interval_t& m1);
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

interval_t LOWER(const interval_t &x) { return bmp::lower(x); }
interval_t UPPER(const interval_t &x) { return bmp::upper(x); }

void TO_LOWER(interval_t &x) { x = bmp::lower(x); }
void TO_UPPER(interval_t &x) { x = bmp::upper(x); }

// Chebyshev interpolation class in interval arithmetic
using cheb_t = LSV::interval_cheb_t;

// get a rigorous but rough estimate for min and max of a cheb approximation
interval_t cheb_range(const cheb_t &p);

int main() {
    using std::cout;

    for (interval_t gamma = 1 + interval_t(-1, 1) * 1e-80; gamma <= 5; gamma += 50) {
        //cout << "gamma: " << gamma << "\n";
        LSV lsv(gamma);

        // derivatives of f(x) = x (1 + 2^gamma x^gamma)
        auto fp = [&gamma](int n, const interval_t &x) -> interval_t {
            const interval_t c = pow(2, gamma) * (gamma + 1);
            if (n == 0) {
                return x * (1 + pow(2 * x, gamma));
            } else if (n == 1) {
                return 1 + c * pow(x, gamma);
            } else if (n == 2) {
                return c * gamma * pow(x, gamma - 1);
            } else if (n == 3) {
                return c * gamma * (gamma - 1) * pow(x, gamma - 2);
            } else if (n==4) {
                return c * gamma * (gamma - 1) * (gamma - 2) * pow(x, gamma - 3);
            } else {
                assert(n >= 0 && n <= 4);
                return 0;
            }
        };

        const auto h_meta = lsv.h_meta();
        const cheb_t h_cheb = h_meta.h;

        const interval_t HALF = interval_t(1) / 2,
              h_cheb_err_A = h_meta.err,
              rho_A = h_meta.rho_A,
              // distance from 1/2 to the boundary of A
              dist_half_A = (rho_A + 1 / rho_A - 2) / 8;
        // sup norm of h on A
        interval_t h_sup_A = UPPER(h_cheb.ellipse_norm(rho_A) + h_cheb_err_A);

        cout << "\nrho_A: " << rho_A
            << "\nERROR: " << h_cheb_err_A << "\n";

        // h(x) on (0,1] with error bounds
        auto H = [&lsv, &h_cheb, &h_cheb_err_A, &HALF] (const interval_t &x) -> interval_t {
            if (LOWER(x) >= HALF) {
                return h_cheb.value(x) + interval_t(-1, 1) * h_cheb_err_A;
            } else {
                interval_t r;
                // sum without an error
                interval_t S = lsv.cheb_sum(x).transpose() * h_cheb.coef();
                // the error
                interval_t E = lsv.eps_sum(x, h_cheb_err_A);

                r = (S + E) / 2;
                return r;
            }
        };

        // a sanity check for h
        {
            cout << "\nSanity check for h:\n";
            interval_t x = 0.001,
                       hx = H(x);
            while (UPPER(x) < HALF) {
                auto Y = lsv.left(x);
                interval_t hy = H((Y(0) + 1) / 2) / 2 + hx / Y(1);
                x = Y(0);
                hx = hy;
            }
            interval_t Hx = H(x);
            cout << "  x: " << x << ", h(x) computed from infinite sum: " << hx
                << "  h(x) from Chebyshev: " << H(x) << "\n"
                << "  diff: " << bmp::width(interval_t(Hx - hx))
                << "\n";
        }

        // derivatives of h on [1/2,1]
        cheb_t hp1_cheb = h_cheb.derivative(),
               hp2_cheb = hp1_cheb.derivative(),
               hp3_cheb = hp2_cheb.derivative();
        interval_t hp1_cheb_err = UPPER(1 * h_cheb_err_A / pow(dist_half_A, 1)),
                   hp2_cheb_err = UPPER(2 * h_cheb_err_A / pow(dist_half_A, 2)),
                   hp3_cheb_err = UPPER(6 * h_cheb_err_A / pow(dist_half_A, 3));
        auto Hp1 = [&hp1_cheb, &hp1_cheb_err, &HALF] (const interval_t &x) {
            assert(LOWER(x) >= HALF && UPPER(x) <= 1);
            return hp1_cheb.value(x) + interval_t(-1,1) * hp1_cheb_err;
        };
        auto Hp2 = [&hp2_cheb, &hp2_cheb_err, &HALF] (const interval_t &x) {
            assert(LOWER(x) >= HALF && UPPER(x) <= 1);
            return hp2_cheb.value(x) + interval_t(-1,1) * hp2_cheb_err;
        };
        auto Hp3 = [&hp3_cheb, &hp3_cheb_err, &HALF] (const interval_t &x) {
            assert(LOWER(x) >= HALF && UPPER(x) <= 1);
            return hp3_cheb.value(x) + interval_t(-1,1) * hp3_cheb_err;
        };

        // and sup norms of derivatives on [1/2,1]
        interval_t
            hp1_sup = UPPER(hp1_cheb.ellipse_norm(1) + hp1_cheb_err),
            hp2_sup = UPPER(hp2_cheb.ellipse_norm(1) + hp2_cheb_err),
            hp3_sup = UPPER(hp3_cheb.ellipse_norm(1) + hp3_cheb_err);

        cout << "\n\ngamma is: " << gamma << " (interval width: " << bmp::width(gamma) << ")"
            << "\nTarget precision: " << std::floor(PREC * std::log(2.0) / std::log(10.0)) << " decimal digits, "
            << "working precision: " << LSV::DIGITS << " decimal digits."
            << "\nAll results below hold up to rounding errors when printing numbers.\n";

        interval_t x_star,
                   F1_minus, F2_minus, F3_minus,
                   F1_plus, F2_plus, F3_plus;
        {
            //cout << "\n\nComputing F_1, F_2, F_3:\n";
            const auto abel = lsv.abel_rough_meta();

            interval_t delta_hash = abel.dist_P1_to_A;
            //cout << "  delta_hash: " << delta_hash << "\n";

            // sup norm of h', h'' delta_hash-inside A (could be done better)
            interval_t
                hp1_hash = hp1_cheb.ellipse_norm(rho_A) + 1 * h_cheb_err_A / pow(delta_hash, 1),
                hp2_hash = hp2_cheb.ellipse_norm(rho_A) + 2 * h_cheb_err_A / pow(delta_hash, 2);
            //cout << "  bound on sup norm of h in A: " << h_cheb_err_A
            //    << ", and h' delta_hash-inside A: " <<  hp1_hash
            //    << ", and h'' delta_hash-inside A: " << hp2_hash
            //    << "\n";

            // choose nu arbitrarily
            const interval_t nu = 4,
                  abel_am1 = abel.coef(0),
                  varkappa1 = LOWER(abel.varkappa0 * 2 * (abel_am1 * abel.r1 + nu) * nu
                          / (pow(nu, 2) + pow(abel_am1 * abel.r1 + nu, 2))
                          );
            //cout << "  nu: " << nu
            //    << ", varkappa1: " << varkappa1 << "\n";

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

            //cout << "  R: " << R
            //    << ", r_star: " << r_star
            //    << ", x_star: " << x_star << "\n";

            interval_t C_psi = h_cheb_err_A * (abel_am1 + abel.C1)
                / (abel_am1 - abel.C1)
                * pow(R, - (1 + 1 / gamma));
            interval_t R_L = 6 * C_psi
                / (pow(2 * lsv.pi_ * varkappa1, 3) * pow(nu, 2));
            R_L = UPPER(R_L);

            interval_t C_Delta =
                abs(gamma - 1) / (gamma * r_star * (abel_am1 - abel.C1))
                + (2 * abel_am1 + abel.C2) / (r_star * pow(abel_am1 - abel.C1, 2));
            interval_t Delta_star = UPPER(
                    (gamma / 2 * abel.r1 * abel.C1 + interval_t(1) / 4
                     + C_Delta / 24) * h_cheb_err_A
                    + pow(r_star, - 1 - 1 / gamma) / (48 * gamma * (abel_am1 - abel.C1)) * hp1_hash
                    + R_L / 2
                    );

            //cout << "  R_L: " << R_L << ", Delta_star: " << Delta_star << "\n";
            {
                interval_t theta_star = 1 / sin(min(
                            lsv.pi_ / 2, lsv.pi_ / (4 * gamma)
                            ));
                interval_t D0 = gamma * abel_am1 * hp2_hash / 48,
                           D1 = D0 * (gamma + 4),
                           D2 = D0 * (pow(gamma, 2) + 9 * gamma + 14),
                           D3 = D0 * (pow(gamma, 3) + 15 * pow(gamma, 2) + 56 * gamma + 48);
                interval_t zh = Delta_star * pow(x_star, gamma);
                interval_t M0 = gamma * abel_am1 * H(HALF) / 2,
                           M1 = gamma * abel_am1 * Hp1(HALF) / 8;
                interval_t delta0 = abs(M1) * x_star + D0 * pow(x_star, 2) + zh,
                           delta1 = abs(1 - gamma) * abs(M1) * x_star
                               + D1 * pow(x_star, 2) +  zh * theta_star,
                           delta2 = gamma * abs(1 - gamma) * abs(M1) * x_star
                               + D2 * pow(x_star, 2) +  2 * zh * pow(theta_star, 2),
                           delta3 = gamma * (gamma + 1) * abs(1 - gamma) * abs(M1) * x_star
                               + D3 * pow(x_star, 2) +  6 * zh * pow(theta_star, 3);

                F1_minus = LOWER((gamma * M0 - delta1) / (M0 + delta0));
                F2_minus = LOWER((gamma * (gamma + 1) * M0 - delta2) / (M0 + delta0));
                F3_minus = LOWER((gamma * (gamma + 1) * (gamma + 2) * M0 - delta3) / (M0 + delta0));

                assert(LOWER(M0 - delta0) > 0);

                F1_plus = UPPER((gamma * M0 + delta1) / (M0 - delta0));
                F2_plus = UPPER((gamma * (gamma + 1) * M0 + delta2) / (M0 - delta0));
                F3_plus = UPPER((gamma * (gamma + 1) * (gamma + 2) * M0 + delta3) / (M0 - delta0));
            }

            // Verify the lower bounds are strictly positive
            assert(F1_minus > 0 && F2_minus > 0 && F3_minus > 0);

            interval_t F1_minus_x = LOWER(F1_minus / pow(x_star, 1)),
                       F2_minus_x = LOWER(F2_minus / pow(x_star, 2)),
                       F3_minus_x = LOWER(F3_minus / pow(x_star, 3));
            interval_t F1_plus_x = UPPER(F1_plus / pow(x_star, 1)),
                       F2_plus_x = UPPER(F2_plus / pow(x_star, 2)),
                       F3_plus_x = UPPER(F3_plus / pow(x_star, 3));

            cout << "\nLemma. For all 0 < x <= " << LOWER(x_star) << ",\n"
                << "  - " << F1_plus << " / x <= h'(x) / h(x) <= - " << F1_minus << " / x\n"
                << "  "   << F2_minus << " / x^2 <= h''(x) / h(x) <= "  << F2_plus << " / x^2\n"
                << "  - " << F3_plus << " / x^3 <= h'''(x) / h(x) <= - " << F3_minus << " / x^3\n";

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
        } // computing F_k

        {
            cout << "Lemma. h'(x) < 0, h''(x) > 0, h'''(x) < 0 on [1/2,1].\n"
                << "Proving...\n";
            interval_t
                hp1_max = UPPER(cheb_range(hp1_cheb) + hp1_cheb_err),
                hp2_min = LOWER(cheb_range(hp2_cheb) - hp2_cheb_err),
                hp3_max = UPPER(cheb_range(hp3_cheb) + hp3_cheb_err);
            assert(hp1_max < 0 && hp2_min > 0 && hp3_max < 0);
            cout << "  ... done. h'(x) <= " << hp1_max << ", h''(x) >= " << hp2_min
                << ", h'''(x) <= " << hp3_max << " ∎\n";
            cout << "\nSanity check: h'(1) is " << Hp1(1)
                << ", h''(1) is " << Hp2(1)
                << ", h'''(1) is " << Hp3(1)
                << "\n\n";
        }

        {
            cout << "Lemma. h'(x) < 0, h''(x) > 0, h'''(x) < 0 on (0, 1].\n"
                << "Proving...\n";

            // get a decreasing sequence P where first N points
            // are 1, 1 - 0.5/N, ..., 0.5 + 0.5/N, and P(k+N) is the preimage of P(k)
            const int N = 16;
            std::vector<interval_t> P;

            // first N points
            for (int k = 0; k < N; k++) {
                P.push_back(1 - interval_t(k) / (2 * N));
            }
            // point N+1 is exactly 1/2
            P.push_back(HALF);

            // the rest, with at least N + 3 points in (0, x_star]
            // P(P_bound) and onward are in (0, x_star]
            int P_bound = -1;
            while (true) {
                int k = P.size();
                P.push_back(lsv.left_inv(P[k - N])(0));

                if (P_bound == -1 && UPPER(P.back()) <= LOWER(x_star))
                    P_bound = k;
                if (P_bound != -1 && P.size() >= P_bound + N + 3)
                    break;
            }
            const int M = P.size();

            // compute values of h iteratively from the left
            VectorXi H_val(M);

#pragma omp parallel for
            for (int k = M - 1; k >= M - N; k--) {
                H_val(k) = H(P[k]);
            }

            for (int k = M - N - 1; k >= 0; k--) {
                interval_t f_prime = lsv.left(P[k + N])(1);
                H_val(k) = HALF * H((1 + P[k]) / 2) + H_val(k + N) / f_prime;
            }

            // array of bounds (C_0, C_1, C_2, C_3),
            // so that (-1)^k h^{(k)} \geq C_k
            // C(k) are the bounds that hold on the interval [P(k+1), P(k)]
            MatrixXi C(M - 1, 4);

            // init in (0, x_star] using the F1--F3 bounds
            for (int k = M - 2; k >= P_bound; k--) {
                interval_t q = P[k];
                interval_t hq = LOWER(H_val(k));
                C(k, 0) = hq;
                C(k, 1) = LOWER(F1_minus / q * hq);
                C(k, 2) = LOWER(F2_minus / pow(q, 2) * hq);
                C(k, 3) = LOWER(F3_minus / pow(q, 3) * hq);
            }

            auto max0 = [](const interval_t &x) -> interval_t {
                return LOWER(x) > 0 ? x : 0;
            };

            // apply Lemma lem:hmon from left to right
            for (int k = P_bound - 1; k >= 0; k--) {
                const interval_t
                    &o = P[k + N + 2],
                    &p = P[k + N + 1],
                    &q = P[k + N];

                interval_t pq(LOWER(p), UPPER(q));

                const interval_t
                    &ho = H_val(k + N + 2),
                    &hp = H_val(k + N + 1);
                    //&hq = H_val(k + N);

                interval_t hp_u = UPPER(hp),
                           h_fq = (1 + P[k]) / 2;

                interval_t d1 = fp(1, pq),
                           d2 = fp(2, pq),
                           d3 = fp(3, pq),
                           d4 = fp(4, pq);

                interval_t f1 = LOWER(d1),
                           F1 = UPPER(d1),
                           f2 = LOWER(d2),
                           F2 = UPPER(d2),
                           f3m = max0(-UPPER(d3)),
                           F3p = max0(UPPER(d3)),
                           f4p = max0(LOWER(d4)),
                           F4m = max0(-LOWER(d4));

                const interval_t
                    &c0 = C(k + N, 0),
                    &c1 = C(k + N, 1),
                    &c2 = C(k + N, 2),
                    &c3 = C(k + N, 3);

                // upper bound on -h'
                interval_t C_prime = UPPER(-(hp - ho) / (p - o));

                C(k, 0) = LOWER(H_val(k));

                C(k, 1) = LOWER( c1/pow(F1, 2) + c0*f2/pow(F1, 3) - interval_t(1)/4 * Hp1(h_fq) );

                C(k, 2) = LOWER( c2/pow(F1, 3) + 3*c1*f2/pow(F1, 4) + 3*c0*pow(f2, 2)/pow(F1, 5)
                        + c0*f3m/pow(F1, 4) - hp_u*F3p/pow(f1, 4) + interval_t(1)/8 * Hp2(h_fq) );

                C(k, 3) = LOWER( c3/pow(F1, 4) + 6*c2*f2/pow(F1, 5) + 15*c1*pow(f2, 2)/pow(F1, 6)
                        + 15*c0*pow(f2, 3)/pow(F1, 7) + (4*c1*f3m + c0*f4p)/pow(F1, 5)
                        + 10*c0*f2*f3m/pow(F1, 6) - (4*C_prime*F3p + hp_u*F4m)/pow(f1, 5)
                        - 10*hp_u*F2*F3p/pow(f1, 6) - interval_t(1)/16 * Hp3(h_fq) );

                assert(C(k, 0) > 0 && C(k, 1) > 0 && C(k, 2) > 0 && C(k, 3) > 0);
            }

            // find the maximum of h'''
            interval_t min_hp3 = -C(0, 3);
            for (int k = 1; k < M - 1; k++) {
                min_hp3 = max(-C(k, 3), min_hp3);
            }

            cout << "  ... done. h'''(x) <= " << min_hp3 << " ∎\n";

            // lower bound (h'/h)' and upper bound (h''/h)'
            interval_t min_hp_h_prime,
                       max_hpp_h_prime;

            {
                // start with a bound on (0, x_star]
                interval_t B1 = F2_minus - pow(F1_plus, 2);
                min_hp_h_prime = LOWER(B1 / pow(x_star, 2));

                interval_t B2 = F3_minus - F2_plus * F1_plus;
                max_hpp_h_prime = UPPER(-B2 / pow(x_star, 3));

                // do [x_star, 1], iterating over [P[k+1], P[k]]
                for (int k = 0; k <= P_bound - 1; k++) {
                    const interval_t
                        &n = P[k + 3],
                        &o = P[k + 2],
                        &p = P[k + 1];
                        //&q = P[k];
                    const interval_t
                        &hn = H_val(k + 3),
                        &ho = H_val(k + 2),
                        &hp = H_val(k + 1);
                        //&hq = H_val(k);
                    // Local lower bounds for h, h'', -h'''
                    const interval_t
                        &c0 = C(k, 0),
                        &c2 = C(k, 2),
                        &c3 = C(k, 3);

                    // max of -h' as a finite difference on the left
                    interval_t M1 = UPPER(- derivative_range(o, p, ho, hp, 0));

                    // max of h'' as a divided difference on the left
                    interval_t M2 = UPPER(second_derivative_range(n, o, p, hn, ho, hp, 0));

                    // bounds for the derivative ratios on [p, q]
                    interval_t pq_hp_h_prime = LOWER((c2 * c0 - pow(M1, 2)) / pow(hp, 2)),
                               pq_hpp_h_prime = UPPER((-c3 * c0 + M2 * M1) / pow(hp, 2));

                    min_hp_h_prime = min(min_hp_h_prime, pq_hp_h_prime);
                    max_hpp_h_prime = max(max_hpp_h_prime, pq_hpp_h_prime);
                }
            }

            cout << "\nTheorem. On (0,1]:\n"
                 << "  (h'/h)' >= " << min_hp_h_prime << "\n"
                 << "  (h''/h)' <= " << max_hpp_h_prime << "\n";

            cout << "\nSanity check:\n"
                << "  (h'/h)'(1): " << Hp2(1) / H(1) - pow(Hp1(1) / H(1), 2) << "\n"
                << "  (h''/h)'(1): " << Hp3(1) / H(1) - Hp1(1) * Hp2(1) / pow(H(1), 2) << "\n";

            assert(min_hp_h_prime > 0 && max_hpp_h_prime < 0);

        }

        cout << "\n";
    }

    return 0;
}

interval_t function_range(const interval_t& x1, const interval_t& x2,
        const interval_t& y1, const interval_t& y2,
        const interval_t& m1) {
    assert(x1 <= x2);
    interval_t mean_y = (y1 + y2) / 2,
               max_dev = m1 * (x2 - x1) / 2;

    return interval_t(LOWER(mean_y - max_dev), UPPER(mean_y + max_dev));
}

interval_t derivative_range(const interval_t& x1, const interval_t& x2,
        const interval_t& y1, const interval_t& y2,
        const interval_t& m2) {
    assert(x1 <= x2);

    interval_t dx = x2 - x1;
    interval_t mean_slope = (y2 - y1) / dx;

    // The maximum deviation of f'(x) from f'(c) on [x1, x2] is m2 * dx
    interval_t max_dev = m2 * dx;
    max_dev = UPPER(max_dev);
    interval_t error_term(-max_dev, max_dev);

    return mean_slope + error_term;
}

interval_t second_derivative_range(const interval_t& x1, const interval_t& x2, const interval_t& x3,
        const interval_t& y1, const interval_t& y2, const interval_t& y3,
        const interval_t& m3) {
    assert(x1 <= x2 && x2 <= x3);

    interval_t dx12 = x2 - x1;
    interval_t dx23 = x3 - x2;
    interval_t dx13 = x3 - x1;

    interval_t slope12 = (y2 - y1) / dx12;
    interval_t slope23 = (y3 - y2) / dx23;

    // Using divided differences, f''(c) = 2 * f[x1, x2, x3]
    interval_t mean_second_deriv = interval_t(2) * (slope23 - slope12) / dx13;

    interval_t max_dev = m3 * dx13;
    max_dev = UPPER(max_dev);
    interval_t error_term(-max_dev, max_dev);

    return mean_second_deriv + error_term;
}

interval_t cheb_range(const cheb_t &p) {
    cheb_t dp = p.derivative();
    interval_t d_sup = dp.ellipse_norm(1);

    VectorXi nodes = p.nodes();
    int N = nodes.size();

    VectorXi X(N + 2),
             Y(N + 2);

    X(0) = interval_t(p.a());
    Y(0) = p.value(X(0));

    for(int i = 0; i < N; i++) {
        // reverse because nodes are descending
        X(i + 1) = nodes(N - 1 - i);
        Y(i + 1) = p.value(X(i + 1));
    }

    X(N + 1) = interval_t(p.b());
    Y(N + 1) = p.value(X(N + 1));

    interval_t min_v = LOWER(Y(0));
    interval_t max_v = UPPER(Y(0));

    for(int i = 0; i < N + 1; i++) {
        interval_t range = function_range(X(i), X(i+1), Y(i), Y(i+1), d_sup);
        min_v = min(min_v, LOWER(range));
        max_v = max(max_v, UPPER(range));
    }

    return interval_t(min_v, max_v);
}

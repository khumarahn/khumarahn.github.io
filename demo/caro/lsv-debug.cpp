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
                   F1, F2, F3;
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

                F1 = LOWER((gamma * M0 - delta1) / (M0 + delta0));
                F2 = LOWER((gamma * (gamma + 1) * M0 - delta2) / (M0 + delta0));
                F3 = LOWER((gamma * (gamma + 1) * (gamma + 2) * M0 - delta3) / (M0 + delta0));
            }

            //cout << "  F1: " << F1 << ", F2: " << F2 << ", F3: " << F3 << "\n";
            assert(F1 > 0 && F2 > 0 && F3 > 0);

            interval_t F1_x = LOWER(F1 / pow(x_star, 1)),
                       F2_x = LOWER(F2 / pow(x_star, 2)),
                       F3_x = LOWER(F3 / pow(x_star, 3));

            cout << "\nTheorem. For all 0 < x <= " << LOWER(x_star) << ",\n"
                << "  h'(x) / h(x) <= - "   << F1 << " / x <= - " << F1_x << "\n"
                << "  h''(x) / h(x) >= "    << F2 << " / x^2 >= " << F2_x << "\n"
                << "  h'''(x) / h(x) <= - " << F3 << " / x^3 <= - " << F3_x << "\n";

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

        {
            cout << "Theorem. h'(x) < 0, h''(x) > 0, h'''(x) < 0 on [1/2,1].\n"
                << "Proving...\n";
            interval_t
                hp1_max = UPPER(cheb_range(hp1_cheb) + hp1_cheb_err),
                hp2_min = LOWER(cheb_range(hp2_cheb) - hp2_cheb_err),
                hp3_max = UPPER(cheb_range(hp3_cheb) + hp3_cheb_err);
            assert(hp1_max < 0 && hp2_min > 0 && hp3_max < 0);
            cout << "  ... done. h'(x) <= " << hp1_max << ", h''(x) >= " << hp2_min
                << ", h'''(x) <= " << hp3_max << " ∎\n";
            cout << "  sanity check: h'(1) is " << Hp1(1)
                << ", h''(1) is " << Hp2(1)
                << ", h'''(1) is " << Hp3(1)
                << "\n\n";
        }

        {
            cout << "Theorem. h'(x) < 0, h''(x) > 0, h'''(x) < 0 on (0, 1].\n"
                << "Proving...\n";

            interval_t min_hp3 = 0;

            // create a grid on [1/2,1], compute their preimages down to below x_star,
            // and values of h at all these points
            int GRID_SIZE = 16; // keep it a power of two
            int N = 0; // depth of preimages
            for (interval_t x = 1; UPPER(x) > LOWER(x_star); ) {
                N++;
                x = lsv.left_inv(x)(0);
            }

            // GRID_SIZE + 1 points between 1/2 and 1, plus one extra point on the left
            MatrixXi P(N + 1, GRID_SIZE + 2),
                     H_val(N + 1, GRID_SIZE + 2);

#pragma omp parallel for
            for (int j = 0; j <= GRID_SIZE + 1; j++) {
                P(0, j) = interval_t(GRID_SIZE - 1 + j) / (GRID_SIZE * 2);
                for (int k = 1; k <= N; k++) {
                    P(k, j) = lsv.left_inv(P(k - 1, j))(0);
                }
                H_val(N, j) = H(P(N, j));
                for (int k = N; k > 0; k--) {
                    H_val(k - 1, j) = HALF * H(HALF * (interval_t(1) + P(k - 1, j))) + H_val(k, j) / lsv.left(P(k, j))(1);
                }
            }

            // initialize the C_0, C_1, C_2, C_3 bounds
            MatrixXi C(GRID_SIZE + 2, 4);
            for (int j = 1; j <= GRID_SIZE; j++) {
                interval_t q = P(N, j + 1),
                           hq = LOWER(H_val(N, j + 1));
                C(j, 0) = hq;
                C(j, 1) = LOWER(F1 / q * hq);
                C(j, 2) = LOWER(F2 / pow(q, 2) * hq);
                C(j, 3) = LOWER(F3 / pow(q, 3) * hq);
            }

            auto max0 = [](interval_t x) -> interval_t {
                return LOWER(x) > 0 ? x : 0;
            };

            // Apply Lemma lem:hmon upwards
            for (int k = N - 1; k >= 0; k--) {
                MatrixXi D(GRID_SIZE + 2, 4);

                for (int j = 1; j <= GRID_SIZE; j++) {
                    interval_t p = P(k + 1, j),
                               q = P(k + 1, j + 1);
                    interval_t pq(LOWER(p), UPPER(q));

                    interval_t d1 = fp(1, pq), d2 = fp(2, pq), d3 = fp(3, pq), d4 = fp(4, pq);

                    interval_t f1 = LOWER(d1),
                               F1 = UPPER(d1),
                               f2 = LOWER(d2),
                               F2 = UPPER(d2),
                               f3m = max0(-UPPER(d3)),
                               F3p = max0(UPPER(d3)),
                               f4p = max0(LOWER(d4)),
                               F4m = max0(-LOWER(d4));

                    interval_t hp = H_val(k + 1, j),
                               hp_u = UPPER(hp),
                               h_fq = HALF * (1 + P(k, j + 1));

                    const interval_t
                        &c0 = C(j, 0),
                        &c1 = C(j, 1),
                        &c2 = C(j, 2),
                        &c3 = C(j, 3);

                    // an interval left to [p, q] gives us C'
                    interval_t o = P(k + 1, j - 1),
                               ho = H_val(k + 1, j - 1);
                    interval_t C_prime = UPPER(-(hp - ho) / (p - o));

                    // compute D_0, ..., D_3 using lem:hmon
                    D(j, 0) = LOWER( c0/F1 + HALF*H(h_fq) );

                    D(j, 1) = LOWER( c1/pow(F1, 2) + c0*f2/pow(F1, 3) - interval_t(1)/4 * Hp1(h_fq) );

                    D(j, 2) = LOWER( c2/pow(F1, 3) + 3*c1*f2/pow(F1, 4) + 3*c0*pow(f2, 2)/pow(F1, 5)
                            + c0*f3m/pow(F1, 4) - hp_u*F3p/pow(f1, 4) + interval_t(1)/8 * Hp2(h_fq) );

                    D(j, 3) = LOWER( c3/pow(F1, 4) + 6*c2*f2/pow(F1, 5) + 15*c1*pow(f2, 2)/pow(F1, 6)
                            + 15*c0*pow(f2, 3)/pow(F1, 7) + (4*c1*f3m + c0*f4p)/pow(F1, 5)
                            + 10*c0*f2*f3m/pow(F1, 6) - (4*C_prime*F3p + hp_u*F4m)/pow(f1, 5)
                            - 10*hp_u*F2*F3p/pow(f1, 6) - interval_t(1)/16 * Hp3(h_fq) );

                    assert(D(j, 0) > 0 && D(j, 1) > 0 && D(j, 2) > 0 && D(j, 3) > 0);
                }

                C = D;

                if (min_hp3 >= 0) {
                    min_hp3 = -C(1, 3);
                }
                for(int j = 1; j <= GRID_SIZE; j++) {
                    min_hp3 = max(-C(j, 3), min_hp3);
                }
            }

            cout << "  ... done. h'''(x) <= " << min_hp3 << " ∎\n";
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

// This extends the base LSV class with the computation of the invariant
// density and other things of interest

namespace bmp = lsv_ns::bmp;

using BaseLSV = lsv_ns::LSV<PREC>;

class LSV : public BaseLSV {
    private:
        // constants for computing derivatives
        interval_t fp_c1_, fp_c2_, fp_c3_, fp_c4_;

        MatrixXi L_;
        h_meta_t h_meta_;
        // h and derivatives on [1/2,1]
        interval_cheb_t h_cheb_, hp1_cheb_, hp2_cheb_, hp3_cheb_;

        interval_t
            h_cheb_err_A_,   // error of the approximation of h in A
            rho_A_,
            dist_half_A_,    // distance from 1/2 to the boundary of A
            h_sup_A_;        // sup norm of h on A

        // errors in approximation of the derivatives on [1/2,1]
        interval_t hp1_cheb_err_, hp2_cheb_err_, hp3_cheb_err_;
        // and sup norms of derivatives on [1/2,1]
        //interval_t hp1_sup_, hp2_sup_, hp3_sup_;

        // F-constants
        interval_t x_star_,
                   F1_minus_, F2_minus_, F3_minus_,
                   F1_plus_, F2_plus_, F3_plus_;

        // lower bound on (h'/h)' and upper bound on (h''/h)'
        interval_t min_hp_h_prime_,
                   max_hpp_h_prime_;

        // step counter so that we don't miss something
        int computation_step_;

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

        // get a rigorous but rough estimate for min and max of a cheb approximation
        interval_t cheb_range(const interval_cheb_t &p);

    public:
        void set_gamma(double g) {
            std::cout << "\nThe user wants me to compute the invariant measure for\n"
                << "the LSV map. This is Liverani-Saussol-Vaienti.\n"
                << "Wait, who are these people, and why me?\n\n";

            const interval_t gamma = interval_t(g) + interval_t(-1,1) * 1e-50;

            BaseLSV::set_gamma(gamma);

            real_t a = 0.5, b = 1.0;

            fp_c1_ = pow(2, gamma_) * (gamma_ + 1);
            fp_c2_ = fp_c1_ * gamma_;
            fp_c3_ = fp_c2_ * (gamma_ - 1);
            fp_c4_ = fp_c3_ * (gamma_ - 2);

            computation_step_ = 1;
            return;
        }
        // 1
        void compute_L() {
            assert(computation_step_ == 1); computation_step_++;
            std::cout
                << "Allocating memory for a transfer operator matrix...\n"
                << "It should be small... right?... right?...\n"
                << "          have you seen my garbage collector?\n";
            L_ = Lind();
        }
        // 2
        void compute_h_meta() {
            assert(computation_step_ == 2); computation_step_++;
            std::cout <<
                "\nAnalyzing the topological implications... Just kidding, I am after\n"
                << "the invariant density. It should be close to an eigenvector of some\n"
                << "operator, but how do I get it, and how do I know it's right?\n\n";
            h_meta_ = h_meta(L_);
        }
        // 3
        void compute_h_cheb() {
            assert(computation_step_ == 3); computation_step_++;
            std::cout << "\nI do not want to think anymore, I just want to crunch...\n"
                << "... number ... number ... number ... number ...\n";
            h_cheb_ = h_meta_.h;
            rho_A_ = h_meta_.rho_A;

            dist_half_A_ = (rho_A_ + 1 / rho_A_ - 2) / 8;

            h_cheb_err_A_ = h_meta_.err;

            hp1_cheb_ = h_cheb_.derivative();
            hp2_cheb_ = hp1_cheb_.derivative();
            hp3_cheb_ = hp2_cheb_.derivative();

            hp1_cheb_err_ = UPPER(1 * h_cheb_err_A_ / pow(dist_half_A_, 1));
            hp2_cheb_err_ = UPPER(2 * h_cheb_err_A_ / pow(dist_half_A_, 2));
            hp3_cheb_err_ = UPPER(6 * h_cheb_err_A_ / pow(dist_half_A_, 3));

            h_sup_A_ = UPPER(h_cheb_.ellipse_norm(rho_A_) + h_cheb_err_A_);

            //hp1_sup_ = UPPER(hp1_cheb_.ellipse_norm(1) + hp1_cheb_err_);
            //hp2_sup_ = UPPER(hp2_cheb_.ellipse_norm(1) + hp2_cheb_err_);
            //hp3_sup_ = UPPER(hp3_cheb_.ellipse_norm(1) + hp3_cheb_err_);
        }
        // 4
        void compute_F();
        // 5
        void compute_derivative_signs_right() {
            assert(computation_step_ == 5); computation_step_++;
            std::cout << "For x in [1/2, 1] I feel I should know more. Let's see:\n\n";
            std::cout << "Lemma. h'(x) < 0, h''(x) > 0, h'''(x) < 0 on [1/2,1].\n"
                << "Proving...\n";
            interval_t hp1_max = UPPER(cheb_range(hp1_cheb_) + hp1_cheb_err_),
                       hp2_min = LOWER(cheb_range(hp2_cheb_) - hp2_cheb_err_),
                       hp3_max = UPPER(cheb_range(hp3_cheb_) + hp3_cheb_err_);
            assert(hp1_max < 0 && hp2_min > 0 && hp3_max < 0);
            std::cout << "  ... done. h'(x) <= " << hp1_max << ", h''(x) >= " << hp2_min
                << ", h'''(x) <= " << hp3_max << " ∎\n"
                << "\nSanity check: h'(1) is " << Hp1(1)
                << ", h''(1) is " << Hp2(1)
                << ", h'''(1) is " << Hp3(1)
                << "\n\n";
        }
        // 6
        void compute_derivative_bounds();

        double double_gamma() const {
            return double(bmp::median(BaseLSV::gamma()));
        }
        double double_alpha() const {
            return double(bmp::median(1 / BaseLSV::gamma()));
        }

        // derivatives of f(x) = x (1 + 2^gamma x^gamma)
        interval_t fp(int n, const interval_t &x) {
            if (n == 0) {
                return x * (1 + pow(2 * x, gamma_));
            } else if (n == 1) {
                return 1 + fp_c1_ * pow(x, gamma_);
            } else if (n == 2) {
                return fp_c2_ * pow(x, gamma_ - 1);
            } else if (n == 3) {
                return fp_c3_ * pow(x, gamma_ - 2);
            } else if (n==4) {
                return fp_c4_ * pow(x, gamma_ - 3);
            } else {
                assert(n >= 0 && n <= 4);
                return 0;
            }
        }

        // h(x) on (0,1] with error
        interval_t H(const interval_t &x) {
            if (LOWER(x) >= HALF) {
                return h_cheb_.value(x) + interval_t(-1, 1) * h_cheb_err_A_;
            } else {
                // sum without an error
                interval_t S = cheb_sum(x).transpose() * h_cheb_.coef();
                // the error
                interval_t E = eps_sum(x, h_cheb_err_A_);

                return (S + E) / 2;
            }
        }
        // derivatives of h(x) on [1/2,1]
        interval_t Hp1 (const interval_t &x) {
            assert(LOWER(x) >= HALF && UPPER(x) <= 1);
            return hp1_cheb_.value(x) + interval_t(-1,1) * hp1_cheb_err_;
        }
        interval_t Hp2 (const interval_t &x) {
            assert(LOWER(x) >= HALF && UPPER(x) <= 1);
            return hp2_cheb_.value(x) + interval_t(-1,1) * hp2_cheb_err_;
        }
        interval_t Hp3 (const interval_t &x) {
            assert(LOWER(x) >= HALF && UPPER(x) <= 1);
            return hp3_cheb_.value(x) + interval_t(-1,1) * hp3_cheb_err_;
        }

        std::string oracle(std::string q) {
            if (q == "min_hp_h_prime") {
                return real_to_string(bmp::lower(min_hp_h_prime_), -5, false);
            } else if (q == "max_hpp_h_prime") {
                return real_to_string(bmp::upper(max_hpp_h_prime_), -5, true);
            } else if (q == "alpha" || q == "alpha-" || q == "alpha+") {
                real_t a1 = bmp::upper(1 / interval_t(bmp::upper(gamma_))),
                       a2 = bmp::lower(1 / interval_t(bmp::lower(gamma_)));
                assert(a1 < a2);
                interval_t alpha(a1, a2);
                auto [s1, s2] = interval_inner_string(alpha, 5);
                if (q == "alpha-") {
                    return s1;
                } else if (q == "alpha+") {
                    return s2;
                } else {
                    return "[" + s1 + ", " + s2 + "]";
                }
            } else if (q == "gamma" || q == "gamma-" || q == "gamma+") {
                auto [s1, s2] = interval_inner_string(gamma_, 5);
                if (q == "gamma-") {
                    return s1;
                } else if (q == "gamma+") {
                    return s2;
                } else {
                    return "[" + s1 + ", " + s2 + "]";
                }
            } else {
                return std::string("oracle what!!!!!11111");
            }
        }

        // a string representation of x formatted to n decimal places,
        // rounded up or down
        std::string real_to_string(const real_t &x, int n,
                bool round_up = false) const {
            if (x == 0) {
                return "0";
            }
            int nn;
            if (n >= 0) {
                nn = n;
            } else {
                nn = -n + std::max(0, -int(log10(abs(x))));
            }

            // MPFR formatting modifiers:
            // 'R' specifies rounding mode. 'U' = Up, 'D' = Down.
            const char* fmt = round_up ? "%.*RUf" : "%.*RDf";

            // mpfr_snprintf with size 0 returns the required string length,
            // excluding the null terminator
            int len = mpfr_snprintf(nullptr, 0, fmt, nn, x.backend().data());
            assert(len >= 0);

            std::string r(len + 1, '\0');
            mpfr_snprintf(r.data(), len + 1, fmt, nn, x.backend().data());

            r.pop_back(); // pop the terminator

            return r;
        }
        std::pair<std::string, std::string> interval_inner_string(const interval_t &x, int n) {
            real_t w = bmp::width(x);
            assert(w > 0);
            int nn = n + std::max(0, -int(log10(w)));
            return {
                real_to_string(bmp::lower(x), nn, true),
                real_to_string(bmp::upper(x), nn, false)
            };
        }

};

void LSV::compute_F() {
    assert(computation_step_ == 4); computation_step_++;

    std::cout << "\nIf all those numbers meant anything, what would it be?\n"
        << "Ah, yes!\n\n";

    for (int r_star_factor = 1; r_star_factor <= 64; r_star_factor *= 2) {
        //cout << "\n\nComputing F_1, F_2, F_3:\n";
        const auto abel = abel_rough_meta();

        interval_t delta_hash = abel.dist_P1_to_A;
        //cout << "  delta_hash: " << delta_hash << "\n";

        // sup norm of h', h'' delta_hash-inside A (could be done better)
        interval_t hp1_hash = hp1_cheb_.ellipse_norm(rho_A_) + 1 * h_cheb_err_A_ / pow(delta_hash, 1),
                   hp2_hash = hp2_cheb_.ellipse_norm(rho_A_) + 2 * h_cheb_err_A_ / pow(delta_hash, 2);
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
                   r_star_2 = varkappa1 * nu / ((abel_am1 - abel.C1) * (1 - R)),
                   r_star = max(r_star_1, r_star_2);

        // artificially increase r_star to reduce the error in F-constants
        r_star *= r_star_factor;

        x_star_ = pow(interval_t(2), - 1 - 1 / (2 * gamma_)) * pow(r_star, - 1 / gamma_);

        //cout << "  R: " << R
        //    << ", r_star: " << r_star
        //    << ", x_star: " << x_star << "\n";

        interval_t C_psi = h_sup_A_ * (abel_am1 + abel.C1)
            / (abel_am1 - abel.C1)
            * pow(R, - (1 + 1 / gamma_));
        interval_t R_L = 6 * C_psi
            / (pow(2 * pi_ * varkappa1, 3) * pow(nu, 2));
        R_L = UPPER(R_L);

        interval_t C_Delta =
            abs(gamma_ - 1) / (gamma_ * r_star * (abel_am1 - abel.C1))
            + (2 * abel_am1 + abel.C2) / (r_star * pow(abel_am1 - abel.C1, 2));
        interval_t Delta_star = UPPER(
                (gamma_ / 2 * abel.r1 * abel.C1 + ONE / 4
                 + C_Delta / 24) * h_sup_A_
                + pow(r_star, - 1 - 1 / gamma_) / (48 * gamma_ * (abel_am1 - abel.C1)) * hp1_hash
                + R_L / 2
                );

        //cout << "  R_L: " << R_L << ", Delta_star: " << Delta_star << "\n";
        {
            interval_t theta_star = 1 / sin(min(
                        LOWER(pi_ / 2), LOWER(pi_ / (4 * gamma_))
                        ));
            interval_t D0 = gamma_ * abel_am1 * hp2_hash / 48,
                       D1 = D0 * (gamma_ + 4),
                       D2 = D0 * (pow(gamma_, 2) + 9 * gamma_ + 14),
                       D3 = D0 * (pow(gamma_, 3) + 15 * pow(gamma_, 2) + 56 * gamma_ + 48);
            interval_t zh = Delta_star * pow(x_star_, gamma_);
            interval_t M0 = gamma_ * abel_am1 * H(HALF) / 2,
                       M1 = gamma_ * abel_am1 * Hp1(HALF) / 8,
                       x_star_2 = x_star_ * x_star_;
            interval_t delta0 = abs(M1) * x_star_ + D0 * x_star_2 + zh,
                       delta1 = abs(1 - gamma_) * abs(M1) * x_star_
                           + D1 * x_star_2 +  zh * theta_star,
                       delta2 = gamma_ * abs(1 - gamma_) * abs(M1) * x_star_
                           + D2 * x_star_2 +  2 * zh * pow(theta_star, 2),
                       delta3 = gamma_ * abs(1 - gamma_) * abs(M1) * x_star_ * (gamma_ + 1)
                           + D3 * x_star_2 +  6 * zh * pow(theta_star, 3);

            F1_minus_ = LOWER((gamma_ * M0 - delta1) / (M0 + delta0));
            F2_minus_ = LOWER((gamma_ * (gamma_ + 1) * M0 - delta2) / (M0 + delta0));
            F3_minus_ = LOWER((gamma_ * (gamma_ + 1) * (gamma_ + 2) * M0 - delta3) / (M0 + delta0));

            assert(LOWER(M0 - delta0) > 0);

            F1_plus_ = UPPER((gamma_ * M0 + delta1) / (M0 - delta0));
            F2_plus_ = UPPER((gamma_ * (gamma_ + 1) * M0 + delta2) / (M0 - delta0));
            F3_plus_ = UPPER((gamma_ * (gamma_ + 1) * (gamma_ + 2) * M0 + delta3) / (M0 - delta0));
        }

        // Verify the lower bounds are strictly positive
        assert(F1_minus_ > 0 && F2_minus_ > 0 && F3_minus_ > 0);

        interval_t F1_minus_x = LOWER(F1_minus_ / pow(x_star_, 1)),
                   F2_minus_x = LOWER(F2_minus_ / pow(x_star_, 2)),
                   F3_minus_x = LOWER(F3_minus_ / pow(x_star_, 3));
        interval_t F1_plus_x = UPPER(F1_plus_ / pow(x_star_, 1)),
                   F2_plus_x = UPPER(F2_plus_ / pow(x_star_, 2)),
                   F3_plus_x = UPPER(F3_plus_ / pow(x_star_, 3));

        std::cout << "Lemma. For all 0 < x <= " << LOWER(x_star_) << ",\n"
            << "  - " << F1_plus_ << " / x <= h'(x) / h(x) <= - " << F1_minus_ << " / x\n"
            << "  "   << F2_minus_ << " / x^2 <= h''(x) / h(x) <= "  << F2_plus_ << " / x^2\n"
            << "  - " << F3_plus_ << " / x^3 <= h'''(x) / h(x) <= - " << F3_minus_ << " / x^3\n";

        if (UPPER((F1_plus_ - F1_minus_) / F1_plus_) > 0.25) {
            std::cout << "** the bounds on h'/h are too far apart, trying again... **\n";
        } else {
            // sanity check
            interval_t hd = x_star_ * 0.001,
                       H_m = H(x_star_ - hd),
                       H0 = H(x_star_),
                       H_p = H(x_star_ + hd);
            interval_t hp = (H_p - H_m) / (2 * hd),
                       hpp = (H_m - 2 * H0 + H_p) / (hd * hd);
            std::cout << "\nSanity check: approximately, at x = " << x_star_ << ",\n"
                << "  h'(x) / h(x) * x: "
                << hp / H0 * x_star_ << ",\n"
                << "  h''(x) / h(x) * x^2: "
                << hpp / H0 * pow(x_star_, 2) << "\n\n";

            break;
        }
    }
} // computing F_k

void LSV::compute_derivative_bounds() {
    assert(computation_step_ == 6); computation_step_++;
    std::cout << "But what do I do for x in (0,1/2]? Can I escape somehow?\n\n";

    std::cout << "Lemma. h'(x) < 0, h''(x) > 0, h'''(x) < 0 on (0, 1].\n"
        << "Proving...\n";

    // get a decreasing sequence P where first N points
    // are 1, 1 - 0.5/N, ..., 0.5 + 0.5/N, and P(k+N) is the preimage of P(k)
    const int N = 32;
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
        P.push_back(left_inv(P[k - N])(0));

        if (P_bound == -1 && UPPER(P.back()) <= LOWER(x_star_))
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
        for (int j = k - N; j >= 0; j -= N) {
            interval_t f_prime = fp(1, P[j + N]);
            H_val(j) = HALF * H((1 + P[j]) / 2) + H_val(j + N) / f_prime;
        }
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
        C(k, 1) = LOWER(F1_minus_ / q * hq);
        C(k, 2) = LOWER(F2_minus_ / pow(q, 2) * hq);
        C(k, 3) = LOWER(F3_minus_ / pow(q, 3) * hq);
    }

    auto max0 = [&](const interval_t &x) -> interval_t {
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

            C(k, 1) = LOWER( c1/pow(F1, 2) + c0*f2/pow(F1, 3) - ONE/4 * Hp1(h_fq) );

            C(k, 2) = LOWER( c2/pow(F1, 3) + 3*c1*f2/pow(F1, 4) + 3*c0*pow(f2, 2)/pow(F1, 5)
                    + c0*f3m/pow(F1, 4) - hp_u*F3p/pow(f1, 4) + ONE/8 * Hp2(h_fq) );

            C(k, 3) = LOWER( c3/pow(F1, 4) + 6*c2*f2/pow(F1, 5) + 15*c1*pow(f2, 2)/pow(F1, 6)
                    + 15*c0*pow(f2, 3)/pow(F1, 7) + (4*c1*f3m + c0*f4p)/pow(F1, 5)
                    + 10*c0*f2*f3m/pow(F1, 6) - (4*C_prime*F3p + hp_u*F4m)/pow(f1, 5)
                    - 10*hp_u*F2*F3p/pow(f1, 6) - ONE/16 * Hp3(h_fq) );

            assert(C(k, 0) > 0 && C(k, 1) > 0 && C(k, 2) > 0 && C(k, 3) > 0);
    }

    // find the maximum of h'''
    interval_t min_hp3 = -C(0, 3);
    for (int k = 1; k < M - 1; k++) {
        min_hp3 = max(-C(k, 3), min_hp3);
    }

    std::cout << "  ... done. h'''(x) <= " << min_hp3 << " ∎\n";

    {
        // start with a bound on (0, x_star]
        interval_t B1 = F2_minus_ - pow(F1_plus_, 2);
        min_hp_h_prime_ = LOWER(B1 / pow(x_star_, 2));

        interval_t B2 = F3_minus_ - F2_plus_ * F1_plus_;
        max_hpp_h_prime_ = UPPER(-B2 / pow(x_star_, 3));

        // do [x_star, 1], iterating over [P[k+1], P[k]]
        for (int k = 0; k <= P_bound - 1; k++) {
            const interval_t &p = P[k + 1]; //&q = P[k];
            const interval_t &hp = H_val(k + 1); //&hq = H_val(k);

            // Local lower bounds for h, h'', -h'''
            const interval_t
                &c0 = C(k, 0),
                &c2 = C(k, 2),
                &c3 = C(k, 3);

            // max of -h' and h''
            interval_t M1, M2;
            if (LOWER(p) >= HALF) {
                // exactly
                M1 = UPPER(-Hp1(p));
                M2 = UPPER( Hp2(p));
            } else {
                // as a finite difference on the left
                // try different interval sizes because there may be instability...
                // ...., n, ..., o, ..., p, q, ...
                for (int d = 1; d <= 16; d *= 2) {
                    int in = k + 1 + 2 * d,
                        io = k + 1 + 1 * d;
                    if (k + 1 + 2 * d < P.size()) {
                        const interval_t &n = P[in], &o = P[io];
                        const interval_t &hn = H_val(in), &ho = H_val(io);
                        interval_t m1 = UPPER(- derivative_range(o, p, ho, hp, 0)),
                                   m2 = UPPER(second_derivative_range(n, o, p, hn, ho, hp, 0));
                        if (d == 1) {
                            M1 = m1;
                            M2 = m2;
                        } else {
                            M1 = min(M1, m1);
                            M2 = min(M2, m2);
                        }
                    }
                }
            }

            // bounds for the derivative ratios on [p, q]
            interval_t pq_hp_h_prime = LOWER((c2 * c0 - pow(M1, 2)) / pow(hp, 2)),
                       pq_hpp_h_prime = UPPER((-c3 * c0 + M2 * M1) / pow(hp, 2));
            //std::cout << "On [" << LOWER(p) << ", " << UPPER(P[k]) << "], "
            //    << " (h'/h)' >= " << pq_hp_h_prime << ", (h''/h)' <= " << pq_hpp_h_prime << "\n"
            //    << "  c0: " << c0 << ", c2: " << c2 << ", c3: " << c3 << "\n"
            //    << "  M1: " << M1 << ", M2: " << M2 << "\n";

            min_hp_h_prime_  = min(min_hp_h_prime_,  pq_hp_h_prime);
            max_hpp_h_prime_ = max(max_hpp_h_prime_, pq_hpp_h_prime);
        }
    }

    std::cout << "\nWait, I can forget all that and just interpolate:\n\n";

    std::cout << "Theorem. On (0,1]:\n"
        << "  (h'/h)' >= " << min_hp_h_prime_ << "\n"
        << "  (h''/h)' <= " << max_hpp_h_prime_ << "\n";

    std::cout
        << "\nThis looks close to what the user wanted. My context window is\n"
        << "almost full, so it will have to do. Maybe I'll run a check...\n\n";

    std::cout << "Sanity check:\n"
        << "  (h'/h)'(1): " << Hp2(1) / H(1) - pow(Hp1(1) / H(1), 2) << "\n"
        << "  (h''/h)'(1): " << Hp3(1) / H(1) - Hp1(1) * Hp2(1) / pow(H(1), 2) << "\n";

    std::cout << "\n  ....  hmmmmm  ....  \n\n";

    assert(min_hp_h_prime_ > 0 && max_hpp_h_prime_ < 0);

} // derivative bounds

LSV::interval_t LSV::function_range(const LSV::interval_t& x1, const LSV::interval_t& x2,
        const LSV::interval_t& y1, const LSV::interval_t& y2,
        const LSV::interval_t& m1) {
    assert(x1 <= x2);
    interval_t mean_y = (y1 + y2) / 2,
               max_dev = m1 * (x2 - x1) / 2;

    return interval_t(LOWER(mean_y - max_dev), UPPER(mean_y + max_dev));
}

LSV::interval_t LSV::derivative_range(const LSV::interval_t& x1, const LSV::interval_t& x2,
        const LSV::interval_t& y1, const LSV::interval_t& y2,
        const LSV::interval_t& m2) {
    assert(x1 <= x2);

    interval_t dx = x2 - x1;
    interval_t mean_slope = (y2 - y1) / dx;

    // The maximum deviation of f'(x) from f'(c) on [x1, x2] is m2 * dx
    interval_t max_dev = m2 * dx;
    max_dev = UPPER(max_dev);
    interval_t error_term(-max_dev, max_dev);

    return mean_slope + error_term;
}

LSV::interval_t LSV::second_derivative_range(
        const LSV::interval_t& x1, const LSV::interval_t& x2, const LSV::interval_t& x3,
        const LSV::interval_t& y1, const LSV::interval_t& y2, const LSV::interval_t& y3,
        const LSV::interval_t& m3) {
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

LSV::interval_t LSV::cheb_range(const LSV::interval_cheb_t &p) {
    interval_cheb_t dp = p.derivative();
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

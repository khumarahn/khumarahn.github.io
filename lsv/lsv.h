#pragma once

#include <concepts>
#include <algorithm>
#include <limits>
#include <vector>
#include <complex>
#include <cmath> // for pow, abs, log, ... for standard types
#include <type_traits>

// We need a version of eigen where the constant PI in kissFFT has
// been fixed and has full precision, like:
// git clone --depth 1 --branch 5.0.1 https://gitlab.com/libeigen/eigen.git eigen
#include <Eigen/Core>

// arbitrary precision and intervals
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpfi.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/eigen.hpp>

// arbitrary precision integers (for Bernoulli)
#include <boost/multiprecision/gmp.hpp>

//#include "clenshawcurtis.h"
#include "cheb.h"
#include "interval-root.h"

// Range checking is disabled if NDEBUG or EIGEN_NO_DEBUG is defined

namespace lsv_ns {

namespace bmp = boost::multiprecision;

template <typename T, typename... TArgs>
concept type_one_of = (std::same_as<T, TArgs> || ...);

using std::max;
using std::pow;
using std::log;
using std::ceil;
using std::floor;
using std::exp;

using cheb_ns::Cheb;

template <int PREC = 64>
class LSV {
    public:
        // PREC is binary target precision;
        // DIGITS is decimal working precision,
        // we need it quite a bit higher than the target precision
        static constexpr int DIGITS = max(56, PREC / 2); // FIXME
        using real_t     = bmp::number<bmp::mpfr_float_backend<DIGITS>>;
        using interval_t = bmp::number<bmp::mpfi_float_backend<DIGITS>>;

        static const interval_t e_;
        static const interval_t pi_;
        static const real_t real_eps_;
        static const real_t real_eps_sqrt_;

        template <typename var_t>
            using Vector2 = Eigen::Matrix<var_t,2,1>;
        template <typename var_t>
            using VectorX = Eigen::Matrix<var_t,Eigen::Dynamic,1>;
        template <typename var_t>
            using MatrixX = Eigen::Matrix<var_t,Eigen::Dynamic,Eigen::Dynamic>;

        using VectorXr  = VectorX<real_t>;
        using Vector2r  = Vector2<real_t>;
        using Vector2i  = Vector2<interval_t>;
        using VectorXi  = VectorX<interval_t>;

        using MatrixXr  = MatrixX<real_t>;
        using MatrixXi  = MatrixX<interval_t>;

        using complex_t = std::complex<real_t>;
        using complex_interval_t = std::complex<interval_t>;

        using Vector2c  = Vector2<complex_t>;
        using Vector2ci = Vector2<complex_interval_t>;
        using VectorXci = VectorX<complex_interval_t>;

        using MatrixXci = MatrixX<complex_interval_t>;

        using interval_cheb_t = Cheb<interval_t>;

        // various constants associated with the Abel function
        // \tA = \tilde{A} and \hA = \hat{A},
        // \tA_n(t) = a_{-1} + a_\ell + a_0 + a_1 t + ... + a_n t^n
        struct abel_meta_t {
            // n
            int n;
            // array of coefficients
            VectorXi coef;
            VectorXr coef_ni;
            // |\tA(t) - \tA_n(t)| \leq C0 |t|^{-n}  when  Re(t) \geq r
            interval_t r, C0;
            // \tA(t) is computed with desired accuracy when Re(t) \geq r_good
            interval_t r_good;
            // |\tA'(t) - a_{-1}| \leq C1 when Re(t) \geq r1,
            // with C1 significantly smaller than a_{-1}
            interval_t r1, C1, am1_minus_C1;
            // |z^3 \hA''(z) - 2 a_{-1}| \leq C2
            interval_t C2;

            // CONSTANTS FOR Euler-Maclaurin
            // nu, a parameter for radius of circles
            interval_t nu;
            // \leq (a_{-1} - C_1) / (a_{-1} + C_1)
            interval_t varkappa0;
            // max of \kappa(w) / \kappa(w + s \zeta),
            // as used in the bound on C_\psi
            interval_t max_kappa_ratio;
            // lower bound on the distance from the preimage of the petal (\cP_1)^{1 / \gamma} under
            // the second branch to the boundary of ellipse A
            interval_t dist_P1_to_A;

            int Nstar;
            // Nstar-preimage of 1
            interval_t x_Nstar;

            int L, M, halfM;
        };

        // error and metadata for the invariant density
        struct h_meta_t {
            interval_cheb_t h;
            interval_t err, rho_A;
        };

    protected:
        interval_t gamma_;
    private:
        interval_t gamma_inv_;

        MatrixXi abel_matrix(int n) const;
        void compute_ellipses();
        abel_meta_t compute_abel_stuff(int n, bool rough = false) const;

        static constexpr int PREC_ = PREC;
        int NEED_DIGITS_;
        int N_;

        abel_meta_t abel_, abel_rough_;

        // angle of sector S_C
        interval_t theta_C_;
        // parameters \rho of the ellipses
        interval_t rho_A_, rho_B_, rho_C_, rho_B_plus_, rho_C_plus_;
        // points, sum at which bounds the transfer operator norm
        interval_t norm_point_A_, norm_point_C_, norm_point_C_plus_;

        VectorXci derivatives_s_,
                  derivatives_c_;

        // Constant part of error in evaluation of sums for
        // first N_ basis Chebyshev polynomials
        VectorXi cheb_sum_small_const_error_;
        // Constant part of error in evaluation of sums of small functions
        interval_t eps_sum_small_const_error_;
        void compute_sum_small_const_error();

        void compute_derivatives_cs() {
            // check that required constants are set
            assert(abel_.halfM > 0 && abel_.M > 0
                    && abel_.Nstar > 0 && abel_.L > 0);

            VectorXci &s = derivatives_s_;
            VectorXci &c = derivatives_c_;

            s.resize(abel_.halfM);
            c.resize(abel_.halfM);

            interval_t tau = abel_.varkappa0 * abel_.nu / e_;

            const VectorXi B2 = bernoulli2k(abel_.L);    // B2(ell) = B_{2 ell}

#pragma omp parallel for schedule(dynamic)
            for (int m = 1; m <= abel_.halfM; m++) {
                interval_t q = 2 * pi_ * (2 * m - 1) / interval_t(2 * abel_.M);
                complex_interval_t qq (0, q);
                s(m - 1) = tau * exp(qq);
                c(m - 1) = 0;
                for (int ell = 1; ell <= abel_.L; ell++) {
                    interval_t w = - q * (2 * ell - 1);
                    complex_interval_t ww (0, w);
                    interval_t www = B2(ell) * pow(tau, -(2*ell - 1)) / (2 * ell);
                    c(m - 1) += www * exp(ww);
                }
            }
            return;
        }

    public:
        LSV() {
            // an empty constructor, nothing initialized,
            // run set_gamma first...
        }
        LSV(const interval_t &gamma) {
            set_gamma(gamma);
        }

        const interval_t ONE = 1, HALF = ONE / 2;
        interval_t LOWER(const interval_t &x) const { return bmp::lower(x); }
        interval_t UPPER(const interval_t &x) const { return bmp::upper(x); }

        void set_gamma(const interval_t &gamma) {
            std::cout << "Setting gamma_, computing Abel function coeffs and constants..."
                << std::flush;

            gamma_ = gamma;

            gamma_inv_ = 1 / gamma_;

            compute_ellipses();

            interval_t mlogeps = log(interval_t(2)) * PREC_;

            N_ = int(ceil(bmp::upper(interval_t(
                        mlogeps / log(rho_C_ / rho_A_)
                        ))));

            {   // NEED_DIGITS_ (approximately)
                // when computing the norm of the transfer operator matrix, we want for j < N:
                //   min(\rho_{C_+}^{-j} \rho_C^j, 10^{-DIGITS} \rho_C^j) < 0.001
                interval_t J = bmp::upper(interval_t(
                            13 / log (rho_C_plus_ / rho_C_)
                            ));
                if (J > N_) J = N_;
                NEED_DIGITS_ = int(ceil(bmp::upper(interval_t(
                            (13 + J * log(rho_C_)) / log(real_t(10))
                            ))));
            }
            if(DIGITS < NEED_DIGITS_) {
                std::cout << "Need more DIGITS: " << NEED_DIGITS_ << "\n";
                assert(!(DIGITS < NEED_DIGITS_));
            }

            int n = int(ceil(bmp::upper(mlogeps))); // TODO what is the best choice?
            abel_ = compute_abel_stuff(n);
            abel_rough_ = compute_abel_stuff(8, true);  // 八八八八八八八八

            compute_derivatives_cs();
            compute_sum_small_const_error();

            interval_t CB_ratio = rho_C_ / rho_B_,
                       CA_ratio = rho_C_ / rho_A_;
            std::cout
                << " done\n"
                << "  gamma_: " << gamma_
                << ", width: " << bmp::width(gamma_) << "\n"
                << "  PREC_: " << PREC_
                << ", DIGITS: " << DIGITS
                << ",  NEED_DIGITS_: " << NEED_DIGITS_ << "\n"
                << "  N_: " << N_ << "\n"
                << "  theta_C_: " << theta_C_
                << ", rho_A_: " << rho_A_ << ", rho_B_: " << rho_B_
                << ", rho_C_: " << rho_C_ << "\n"
                << "  rho_B_plus_: " << rho_B_plus_
                << ", rho_C_plus_: " << rho_C_plus_ << "\n"
                << "  rho_C_ / rho_B_: " << bmp::lower(CB_ratio)
                << ", rho_C_ / rho_A_: " << bmp::lower(CA_ratio) << "\n"
                << "  abel_:\n"
                << "    .coef_ni: " << abel_.coef_ni(0)
                << ", " << abel_.coef_ni(1) << ", " << abel_.coef_ni(2)
                << ", " << abel_.coef_ni(3) << ", " << abel_.coef_ni(4)
                << ", ...\n"
                << "    .n: " << abel_.n << "\n"
                << "    .r: " << abel_.r
                << ", .C0: " << abel_.C0
                << ", .r_good: " << abel_.r_good << "\n"
                << "    .r1: " << abel_.r1
                << ", .C1: " << abel_.C1
                << ", .C2: " << abel_.C2 << "\n"
                << "    .varkappa0: " << abel_.varkappa0 << "\n"
                << "    .am1_minus_C1: " << abel_.am1_minus_C1
                << ", .dist_P1_to_A: " << abel_.dist_P1_to_A << "\n"
                << "    .L: " << abel_.L << ", .M: " << abel_.M << "\n"
                << "    .nu: " << abel_.nu
                << ", .Nstar: " << abel_.Nstar
                << ", .x_Nstar: " << abel_.x_Nstar << "\n"
                << "  abel_rough_:\n"
                << "    .coef_ni: " << abel_rough_.coef_ni(0)
                << ", " << abel_rough_.coef_ni(1) << ", " << abel_rough_.coef_ni(2)
                << ", " << abel_rough_.coef_ni(3)
                << ", ...\n"
                << "    .n: " << abel_rough_.n << "\n"
                << "    .r: " << abel_rough_.r
                << ", .C0: " << abel_rough_.C0
                << ", .r_good: " << abel_rough_.r_good << "\n"
                << "    .r1: " << abel_rough_.r1
                << ", .C1: " << abel_rough_.C1
                << ", .C2: " << abel_rough_.C2 << "\n"
                << "    .varkappa0: " << abel_rough_.varkappa0 << "\n"
                << "    .am1_minus_C1: " << abel_rough_.am1_minus_C1
                << ", .dist_P1_to_A: " << abel_rough_.dist_P1_to_A << "\n\n";

            return;
        }

        int NCheb() const { return N_; };
        int prec() const { return PREC_; };
        int KAbel() const { return abel_.n + 1; };
        int abel_n() const { return abel_.n; };
        interval_t alpha() const { return gamma_inv_; };
        interval_t gamma() const { return gamma_; };
        VectorXi abel_coef() const { return abel_.coef; };

        abel_meta_t abel_meta() { return abel_; };
        abel_meta_t abel_rough_meta() { return abel_rough_; };

        // branches
        Vector2i left(const interval_t &x) const {
            interval_t t = pow(2 * x, gamma_);
            return Vector2i(x * (1 + t), 1 + (gamma_ + 1) * t);
        }
        Vector2i right(const interval_t &x) const {
            return Vector2i(2 * x - 1, 2);
        }

        // inverses
        Vector2i left_inv(const interval_t &x) const {
            assert(bmp::lower(x) >= 0 && bmp::upper(x) <= 1);
            // max derivative is
            interval_t md = left(HALF)(1);
            // an interval enclosing the root
            interval_t guess(
                    bmp::lower(interval_t(x / md)),
                    bmp::upper(x) < HALF ? bmp::upper(x) : HALF
                    );

            auto f = [this, &x] (const interval_t &z) {
                Vector2<interval_t> r = left(z);
                r(0) -= x;
                return r;
            };
            return interval_root_ns::interval_newton(f, guess);
        }

        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> right_inv(const var_t &x) const {
            return Vector2<var_t>(
                    (x + var_t(1)) / var_t(2),
                    var_t(1) / var_t(2)
                    );
        }

        // sum S(x) = \sum_{k \geq 0} \varphi((x_k + 1) / 2) / J_k(x_k)
        // where \varphi is the vector of first N_ Chebyshev basis polynomials
        VectorXi cheb_sum(const interval_t &x) const {
            VectorXi r = VectorXi::Zero(N_);

            // for integrals only
            const interval_cheb_t cheb(interval_t(1) / 2, 1, 1);

            // values of first N Chebyshev polynomials at the
            // preimage on [1/2,1]
            auto phi = [this, &cheb]<typename var_t>(const var_t &x) {
                var_t y = this->right_inv(x)(0);
                // slow trigonometric evaluation is more precise
                // in interval arithmetic
                VectorX<var_t> r = cheb.basis_values_trig(y, N_);
                return r;
            };

            // for a **small** real x, return a rigorous approximation of S(z) = S(x^\gamma)
            // for the vector-valued observable of first N Chebyshev polynomials
            auto S_small = [this, &cheb, &phi](const interval_t &x) -> VectorXi {
                //interval_t z = pow(x, gamma_);
                VectorXi r = VectorXi::Zero(N_);
                Vector2i Ax = this->abel(x);

                {   // integral
                    interval_t y = this->right_inv(x)(0);
                    VectorXi bi = cheb.beta_integral(0, y, N_);
                    r -= 2 * Ax(1) * bi;
                }

                // \varphi(z) / 2
                r += phi(x) / 2;

                {   // derivatives
                    const VectorXci &s = derivatives_s_;
                    const VectorXci &c = derivatives_c_;

                    VectorXi der = VectorXi::Zero(N_);
                    for (int m = 1; m <= abel_.halfM; m++) {
//std::cout << __LINE__
//    << " x: " << x << ", width: " << bmp::width(x)
//    << " Ax(0): " << Ax(0) << ", width: " << bmp::width(Ax(0))
//    << "\n";
                        Vector2ci xm_inv = abel_inv(Ax(0) + s(m-1));
                        der += ((c(m-1) * xm_inv(1)) * phi(xm_inv(0))).real();
                    }

                    der *= - 2 * Ax(1) / abel_.M;
                    r += der;
                }

                r += cheb_sum_small_const_error_;

                return r;
            };

            interval_t xk = x,
                       inv_Jk = 1;

            // add branches naively up to Nstar
            while (bmp::upper(xk) > abel_.x_Nstar) {
                r += phi(xk) * inv_Jk;
                Vector2i y = left_inv(xk);
                xk = y(0);
                inv_Jk *= y(1);
            }

            // add the rest when xk is small
            r += S_small(xk) * inv_Jk;
            return r;
        }

        // a bound on the sum S(x) = \sum_{k \geq 0} \varphi((x_k + 1) / 2) / J_k(x_k)
        // where |\varphi| \leq \eps, and \varphi is analytic in \cP_1
        interval_t eps_sum(const interval_t &x, const interval_t &eps = 1) const {
            interval_t r = 0;

            interval_t xk = x,
                       inv_Jk = 1;

            // add branches naively up to Nstar
            while (bmp::upper(xk) > abel_.x_Nstar) {
                r += inv_Jk;
                Vector2i y = left_inv(xk);
                xk = y(0);
                inv_Jk *= y(1);
            }

            Vector2i Axk = abel(xk);

            // bound S(xk)
            interval_t w_minus_nu = Axk(0) - abel_.nu;
            interval_t tt = abel_t_inv(w_minus_nu)(0);
            assert(tt > abel_.r1);

            //r += abs(Axk(1) * xk) + interval_t(1) / 2;
            //r += eps_sum_small_const_error_;
            interval_t tail_bound = abs(Axk(1) * xk) + interval_t(1) / 2
                + eps_sum_small_const_error_;
            r += tail_bound * inv_Jk;

            r *= eps;
            return interval_t(-1, 1) * bmp::upper(r);
        }

        // Transfer operator for the induced map, the clever one.
        // This returns an approximation of the induced transfer operator by an
        // N_ by N_ matrix acting on the Chebyshev coefficients.
        MatrixXi Lind() const {
            const interval_cheb_t cheb(interval_t(1) / 2, 1, N_);
            const VectorXi x_nodes = cheb.nodes();

            // Compute (L T_n)(x) for n=0,...,N-1 and x in the node points
            MatrixXi L_values(N_, N_);
#pragma omp parallel for schedule(dynamic)
            for (int ix = 0; ix < x_nodes.size(); ix++) {
                interval_t x = x_nodes[ix];
                L_values.col(ix) = cheb_sum(x) / interval_t(2);
            }

            // Approximate the result with Chebyshev polynomials
            MatrixXi R(N_, N_);
            {
                interval_cheb_t c(interval_t(1) / 2, 1, N_);
                for (int k = 0; k < N_; k++) {
                    c.set_from_values(L_values.row(k).transpose());
                    VectorXi cc = c.coef();
                    R.col(k) = cc;
                }
            }
            {   // intersect the elements of R with their theoretical bounds
                interval_t norm_plus = cheb_sum(norm_point_C_plus_)(0),
                           C2N = pow(rho_C_plus_, - 2 * N_),
                           norm_plus_C2N = norm_plus / (1 - C2N);
                norm_plus_C2N = bmp::upper(norm_plus_C2N);
                VectorXi Cj(N_), Bbk(N_);
                for (int k = 0; k < N_; k++) {
                    Cj(k) = pow(rho_C_plus_, k);
                    Bbk(k) = (pow(rho_B_plus_, k) + pow(rho_B_plus_, -k))
                            * norm_plus_C2N;
                    Bbk(k) = bmp::upper(Bbk(k));
                }

                for (int k = 0; k < N_; k++) {
                    for (int j = 0; j < N_; j++) {
                        interval_t w = Bbk(k) * (1 / Cj(j) + Cj(j) * C2N),
                                   ww = interval_t(-1, 1) * bmp::upper(w);
                        R(j, k) = bmp::intersect(R(j, k), ww);
                        assert(bmp::width(R(j,k)) >= 0);
                    }
                }
            }

            return R;
        }

        h_meta_t h_meta() const {
            return h_meta(Lind());
        }
        h_meta_t h_meta(const MatrixXi &L) const {
            h_meta_t meta;

            interval_t a(interval_t(1) / 2),
                       b(1);
            int N = L.cols();

            // Cheb class for misc calculations
            interval_cheb_t cheb(a, b, 1);

            // iota is a vector of values of integrals on [1/2,1] of the basis
            // Chebyshev polynomials
            VectorXi iota = cheb.beta_integral(0, b, N);

            VectorXi u = VectorXi::Zero(N);
            u(0) = 1 / iota(0);

            // Now the invariant density in Chebyshev basis
            // is h = (I - L + u iota)^{-1} u
            MatrixXi S = MatrixXi::Identity(N,N) - L + u * iota.transpose();
            VectorXi hv = interval_root_ns::linear_krawczyk(S, u);

            // h as a function
            meta.h = interval_cheb_t(hv, a, b);

            // *** Now h and L are computed

            // norm of an operator given by a matrix
            // corresponding to Lemma \ref{lem:mnorm}
            auto norm_Q = [] (const MatrixXi& Q, interval_t rho_A, interval_t rho_C) {
                assert(Q.cols() == Q.rows());
                int n = Q.cols();

                VectorXi d_C(n), inv_d_A(n);
                d_C(0) = interval_t(1) / 2;
                inv_d_A(0) = 2;
                for (int k = 1; k < n; k++) {
                    d_C(k) = sqrt(pow(rho_C, 2*k) + pow(rho_C, -2*k)) / 2;
                    inv_d_A(k) = 2 / sqrt(pow(rho_A, 2*k) + pow(rho_A, -2*k));
                }

                interval_t K_AC_sq = 1;
                for (int k = 1; k < n; k++) {
                    interval_t num = pow(pow(rho_A, k) + pow(rho_A, -k), 2);
                    interval_t den = pow(rho_C, 2*k) + pow(rho_C, -2*k);
                    K_AC_sq += num / den;
                }
                interval_t K_AC = sqrt(K_AC_sq);

                MatrixXi DQD(n, n);
                real_t max_DQD_width = 0;
                for (int j = 0; j < n; j++) {
                    for (int k = 0; k < n; k++) {
                        DQD(j, k) = d_C(j) * Q(j, k) * inv_d_A(k);
                        max_DQD_width = max(max_DQD_width, bmp::width(DQD(j, k)));
                    }
                }
                //std::cout << "Max width of the intervals in DQD: " << max_DQD_width
                //    << "  [bad if large]\n";
                assert(max_DQD_width < 0.001);

                interval_t DQD_norm = interval_root_ns::matrix_L2_norm(DQD);

                interval_t r = K_AC * DQD_norm;
                r = bmp::upper(r);
                return r;
            };

            auto n_choose_k = [] (int n, int k) -> interval_t {
                if (k < 0 || k > n) return 0;
                if (k > n / 2)
                    k = n - k;

                interval_t r(1);
                for (int j = 1; j <= k; j++) {
                    r = (r * (n - j + 1)) / j;
                }
                return r;
            };

            // norm of I - \bpi between ellipses with
            // parameters rho_1 and rho_2
            auto norm_I_pi = [N] (interval_t rho_1, interval_t rho_2) {
                assert(rho_1 > rho_2);
                interval_t norm = 8 / (pow(rho_1 / rho_2, N - 1) * log(rho_1 / rho_2));
                norm = bmp::upper(norm);
                return norm;
            };

            // more norms
            interval_t norm_I_pi_C_A = norm_I_pi(rho_C_, rho_A_),
                       norm_I_pi_A_B = norm_I_pi(rho_A_, rho_B_),
                       norm_L_B_A = bmp::upper(cheb_sum(norm_point_A_)(0)),
                       norm_L_B_C = bmp::upper(cheb_sum(norm_point_C_)(0)),
                       norm_L_A_C = norm_L_B_C,
                       norm_u_iota_I_pi_A = norm_I_pi(rho_A_, 1);

            interval_t eps = norm_I_pi_C_A * norm_L_A_C
                + norm_I_pi_C_A * norm_L_B_C * norm_I_pi_A_B
                + norm_L_B_A * norm_I_pi_A_B
                + norm_u_iota_I_pi_A;
            eps = bmp::upper(eps);

            interval_t eps_prime = norm_I_pi_C_A * norm_L_A_C;
            eps_prime = bmp::upper(eps_prime);

            // print all norms above:
            std::cout << "Computing norms and eps:\n"
                << "  norm_I_pi_C_A: " << norm_I_pi_C_A
                << ", norm_I_pi_A_B: " << norm_I_pi_A_B << "\n"
                << "  norm_L_B_A: " << norm_L_B_A
                << ", norm_L_B_C: " << norm_L_B_C
                << ", norm_L_A_C: " << norm_L_A_C << "\n"
                << "  norm_u_iota_I_pi_A: " << norm_u_iota_I_pi_A << "\n"
                << "  eps: " << eps
                << ", eps_prime: " << eps_prime << "\n";

            assert(eps < 0.1);

            const MatrixXi bDelta = L - u * iota.transpose();


            std::vector<interval_t> delta;
            std::vector<interval_t> norm_bDelta;
            norm_bDelta.push_back(1);

            int n = 0;
            for (MatrixXi bDelta_n = bDelta; ; n++) {
                interval_t d = 0;
                for (int k = 0; k <= n; k++)
                    d += n_choose_k(n, k) * norm_bDelta[k] * pow(eps, n - k);
                delta.push_back(d);

                if (bmp::upper(delta[n]) < 0.5 || n > 4)
                    break;

                norm_bDelta.push_back(norm_Q(bDelta_n, rho_A_, rho_C_));
                bDelta_n *= bDelta;
            }

            assert(delta[n] < 0.5);
            std::cout << "  delta[" << n << "] = " << delta[n] << " < 0.5\n";

            meta.rho_A = rho_A_;

            interval_t norm_h_A = meta.h.ellipse_norm(rho_A_);
            norm_h_A = bmp::upper(norm_h_A);
            std::cout << "Bound on norm of h in A: " << norm_h_A << "\n";

            {   // error
                interval_t err = 0;
                for (int k = 0; k < n; k++)
                    err += delta[k];

                err *= 2 * eps_prime * norm_h_A;
                meta.err = bmp::upper(err);
            }

            //std::cout << "meta.err: " << meta.err << "\n";

            return meta;
        }

        // *** ABEL FUNCTIONS ***
        // * sum for \tA(t), or abel_t, with a provided array of coefficients, without the error term
        template <typename var_t, typename coef_t> requires
            (type_one_of<var_t, real_t, complex_t> && type_one_of<coef_t, VectorXr>) ||
            (type_one_of<var_t, interval_t, complex_interval_t> && type_one_of<coef_t, VectorXi>)
        Vector2<var_t> abel_t_sum(const var_t &t, const coef_t &coef) const {
            assert(coef.size() > 3);
            int n = coef.size() - 3;
            var_t A = coef(0) * t + coef(1) * log(t) + coef(2),
                  dA = coef(0) + coef(1) / t;

            const var_t ti = var_t(1) / t;

            var_t A_tail(0), dA_tail(0);

            // Horner's method
            for (int j = n; j >= 1; j--) {
                A_tail = coef(2 + j) + ti * A_tail;
                dA_tail = var_t(j * coef(2 + j)) + ti * dA_tail;
            }
            A += ti * A_tail;
            dA -= ti * ti * dA_tail;

            return Vector2<var_t>(A, dA);
        }
        // * non-interval version
        template <typename var_t> requires type_one_of<var_t, real_t, complex_t>
        Vector2<var_t> abel_t(const var_t &t) const {
            return abel_t_sum(t, abel_.coef_ni);
        }
        // * interval version
        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> abel_t(const var_t &t) const {
            return abel_t_sum(t, abel_.coef) + abel_t_error(t);
        }
        // * interval A(x)
        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> abel(const var_t &x) const {
            var_t t = pow(x, -gamma_);
            Vector2<var_t> A = abel_t(t);
            return Vector2<var_t>(A(0), -gamma_ * (t / x) * A(1));
        }

        // inverse
        // * non-interval version in t, Newton method
        template <typename var_t> requires type_one_of<var_t, real_t, complex_t>
        Vector2<var_t> abel_t_inv(const var_t &a) const {
            var_t t = a / abel_.coef_ni(0);
            Vector2<var_t> A;
            for (int i=0; i < 1000; i++) {
                A = abel_t(t);
                if (abs(A(0) - a) < 1000 * real_eps_) {
                    assert(i < 256);
                    break;
                }
                t -= (A(0) - a) / A(1);
            }
            return Vector2<var_t>(t, var_t(1) / A(1));
        }
        // * real interval version in t
        Vector2i abel_t_inv(const interval_t &a) const {
            Vector2r g = abel_t_inv<real_t>(bmp::median(a));
            // an interval enclosing the root
            // FIXME: look at the constants 16 and 0.0001
            interval_t guess = g(0) + 16 * g(1) * (bmp::width(a) + 0.0001) * interval_t(-1, 1);

            auto f = [this, &a] (const interval_t &t) {
                Vector2i r = abel_t(t);
                r(0) -= a;
                return r;
            };
            return interval_root_ns::interval_newton(f, guess);
        }
        // * complex interval version in t
        Vector2ci abel_t_inv(const complex_interval_t &a) const {
            auto m = [] (complex_interval_t x) {
                using bmp::median;
                return complex_t(median(x.real()), median(x.imag()));
            };
            Vector2c g = abel_t_inv<complex_t>(m(a));
            complex_interval_t g0 = g(0),
                               g1 = g(1);
            // FIXME: look at the constants 16 and 0.0001
            complex_interval_t guess = g0 + interval_t(16) * g1 * (
                    bmp::width(a.real()) + bmp::width(a.imag()) + interval_t(0.001)) *
                complex_interval_t(interval_t(-1,1), interval_t(-1,1));
            auto f = [this, &a] (const complex_interval_t &t) {
                Vector2ci r = abel_t(t);
                r(0) -= a;
                return r;
            };
            return interval_root_ns::complex_krawczyk(f, guess);
        }
        // * interval version in x
        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> abel_inv(const var_t &a) const {
            Vector2<var_t> Ai = abel_t_inv(a);
            const var_t &t = Ai(0), &dAt = Ai(1);
            var_t x = pow(t, -var_t(gamma_inv_));
            return Vector2<var_t>(x, -(x / t) * gamma_inv_ * dAt);
        }

        // error term for abel_t
        template <typename i_t> requires type_one_of<i_t, interval_t, complex_interval_t>
        Vector2<i_t> abel_t_error(const i_t &t) const {

            // check if we are in the right region where we can compute the
            // Abel function accurately
            if constexpr (type_one_of<i_t, complex_interval_t>)
                assert( bmp::lower(t.real()) >= abel_.r_good );
            else
                assert( bmp::lower(t) >= abel_.r_good );

            // ir_t is the underlying real interval type
            using ir_t = std::conditional_t<
                type_one_of<i_t,complex_interval_t>,  interval_t,  i_t
                >;

            ir_t factor0 = abel_.C0 * pow(abs(t), -abel_.n),
                 window0 (-bmp::upper(factor0), bmp::upper(factor0)),
                 factor1 = abel_.C0 * pow(abs(t) - 1, -abel_.n),
                 window1 (-bmp::upper(factor1), bmp::upper(factor1));

            Vector2<i_t> r;
            if constexpr (type_one_of<i_t, complex_interval_t>) {
                r(0) = complex_interval_t(window0, window0);
                r(1) = complex_interval_t(window1, window1);
            } else {
                r(0) = window0;
                r(1) = window1;
            }

            return r;
        }


        // Akiyama–Tanigawa algorithm for even Bernoulli numbers B_{2 k}, k = 0, 1, ..., p
        static VectorXi bernoulli2k(int p) {
            // we work in arbitrary precision integers
            using int_t = bmp::mpz_int;
            using VectorXint = Eigen::Matrix<int_t, Eigen::Dynamic, 1>;

            VectorXi B2(p+1);
            VectorXint A(2*p+1);

            int_t fact = 1;
            for (int m = 0; m <= 2*p; m++) {
                A(m) = fact;
                fact *= m + 1;
                for (int j = m; j >= 1; j--)
                    A(j - 1) = j * ((m + 1) * A(j - 1) - A(j));
                if (m % 2 == 0)
                    B2(m / 2) = interval_t(A(0)) / interval_t(fact);
            }
            return B2;
        }

};  // class LSV

template <int PREC>
void LSV<PREC>::compute_ellipses() {
    using i_t = interval_t;

    // sector
    theta_C_ = min(
            bmp::lower(i_t(pi_ / 2)),
            bmp::lower(i_t(pi_ / (2 * gamma_)))
            );

    i_t sin_theta_C = sin(theta_C_);

    // a rho when the ellipse would touch the sector boundary
    i_t rho_max = sqrt(1 + 8 * pow(sin_theta_C, 2))
        + sqrt(i_t(8)) * sin_theta_C;

    // C: decrease it a bit, so it is away from the sector boundary
    rho_C_ = pow(rho_max, interval_t(7) / 8);
    rho_C_ = bmp::lower(rho_C_);

    // B
    i_t s_B = 1 + (rho_C_ + 1 / rho_C_) / 2;
    rho_B_ = (s_B + sqrt(s_B * s_B - 4)) / 2;
    rho_B_ = bmp::upper(rho_B_);

    // C+
    rho_C_plus_ = exp((4 * log(rho_max) + 1 * log(rho_C_)) / 5);
    rho_C_plus_ = bmp::lower(rho_C_plus_);

    // B+
    i_t s_B_plus = 1 + (rho_C_plus_ + 1 / rho_C_plus_) / 2;
    rho_B_plus_ = (s_B_plus + sqrt(s_B_plus * s_B_plus - 4)) / 2;
    rho_B_plus_ = bmp::upper(rho_B_plus_);

    // A
    // a weighted geometric mean, perhaps maximizing
    // the ratio rho_C_ / rho_A_
    i_t mlogeps = log(i_t(2)) * PREC_;
    i_t wC = 12 / (mlogeps + 12);
    i_t wB = 1 - wC;
    rho_A_ = exp(wB * log(rho_B_) + wC * log(rho_C_));
    rho_A_ = bmp::median(rho_A_);

    auto norm_point = [g = gamma_, g_inv = gamma_inv_, pi = pi_]
        (const interval_t& rho) -> interval_t {
            interval_t a = (rho + 1 / rho) / 8,
                       b = (rho - 1 / rho) / 8,
                       q = 1 / ((g + 1) * pow(2, g)),
                       q2 = q * q;

            int steps = 512;
            interval_t dt = pi / (2 * steps);

            interval_t min_L, min_U;
            for (int i = 0; i < steps; ++i) {
                interval_t t = interval_t(i, i + 1) * dt;

                interval_t x = interval_t(3) / 4 - a * cos(t),
                           y = b * sin(t);

                interval_t r = sqrt(x * x + y * y),
                           phi = atan2(y, x),
                           rg = pow(r, g);

                interval_t dist = sqrt(
                        rg * rg + 2 * q * rg * cos(g * phi) + q2
                        );

                if (i == 0 || bmp::lower(dist) < min_L) min_L = bmp::lower(dist);
                if (i == 0 || bmp::upper(dist) < min_U) min_U = bmp::upper(dist);
            }
            assert(min_L > q);
            return pow(interval_t(min_L, min_U) - q, g_inv);
        };
    norm_point_C_plus_ = norm_point(rho_C_plus_);
    norm_point_C_ = norm_point(rho_C_);
    norm_point_A_ = norm_point(rho_A_);

    assert(0 < rho_B_ && rho_B_ < rho_A_ && rho_A_ < rho_C_ && rho_C_ < rho_C_plus_);
    assert(0 < bmp::lower(norm_point_C_plus_)
            && bmp::lower(norm_point_C_plus_) < bmp::lower(norm_point_C_)
            && bmp::lower(norm_point_C_) < bmp::lower(norm_point_A_));
}

template <int PREC>
LSV<PREC>::abel_meta_t LSV<PREC>::compute_abel_stuff(int n, bool rough) const {
    abel_meta_t abel;

    using r_t = real_t;
    using i_t = interval_t;

    assert(n > 0);
    abel.n = n;

    abel.coef.resize(n + 3);
    abel.coef_ni.resize(n + 3);

    assert(0 < rho_B_ && rho_B_ < rho_A_ && rho_A_ < rho_C_);

    {   // first compute non-constant coefficients (am1, al, a1, a2, ...) of the Abel function
        VectorXi b = VectorXi::Zero(abel.n + 2);
        b(0) = -1;

        MatrixXi A = abel_matrix(n); // lower triangular

        // brute force solve lower triangular matrix
        // instead of using Eigen that can fail in interval arithmetic:
        //VectorXi x = abel_matrix(n).template triangularView<Eigen::Lower>().solve(b);
        //VectorXi x = interval_root_ns::linear_krawczyk(A, b);

        VectorXi xx(abel.n + 2);
        for (int k = 0; k < n + 2; k++) {
            i_t s = 0;
            for (int j = 0; j < k; j++) {
                s += A(k,j) * xx(j);
            }
            xx(k) = (b(k) - s) / A(k,k);
        }

        abel.coef(0) = xx(0);
        abel.coef(1) = xx(1);
        abel.coef(2) = 0;
        for (int j = 2; j <= n + 1; j++) {
            abel.coef(j + 1) = xx(j);
        }
    }

    {   // compute r, C0  // TODO: RECHECK
        auto binom = [](const i_t& x, int k) -> i_t {
            i_t r = 1;
            for (int i = 0; i < k; i++) {
                r *= (x - i) / (i + 1);
            }
            return r;
        };

        i_t b = pow(2, gamma_);
        abel.r = b * (gamma_ + 1);

        i_t z = 1 / abel.r,
            bz = 1 / (gamma_ + 1); // b * z;

        // G(z) evaluated at z = r^{-1}
        i_t G_z = abs(abel.coef(0)) / z * (pow(1 - bz, -gamma_) - 1)
                - abs(abel.coef(1)) * gamma_ * log(1 - bz);
        for (int k = 1; k <= n; k++)
            G_z += abs(abel.coef(2 + k)) * pow(z, k)
                * (pow(1 - bz, -k * gamma_) - 1);

        // coefficients c(j) = c_j
        auto c_val = [&coef = abel.coef, n, &gamma = gamma_, &binom, &b]
            (int j) -> i_t {
            if (j == 0) return 0; // c_0 = 0 by design

            i_t r = coef(0) * pow(b, j + 1) * binom(-gamma, j + 1);

            // (-1)^{j-1}
            i_t sign = (j % 2 == 0) ? -1 : 1;
            r += coef(1) * gamma * sign * pow(b, j) / j;

            int max_k = std::min(n, j - 1);
            for (int k = 1; k <= max_k; k++) {
                r += coef(2 + k) * pow(b, j - k) * binom(k * gamma, j - k);
            }
            return r;
        };

        // coefficients g(j) = g_j
        auto g = [&coef = abel.coef, n, &gamma = gamma_, &binom, &b]
            (int j) -> i_t {
            if (j == 0)
                return abs(coef(0)) * b * gamma;

            i_t r = abs(coef(0)) * pow(b, j + 1)
                * binom(j + gamma, j + 1);
            r += abs(coef(1)) * gamma * pow(b, j) / j;

            int max_k = std::min(n, j - 1);
            for (int k = 1; k <= max_k; k++) {
                r += abs(coef(2 + k)) * pow(b, j - k)
                    * binom(j - k + k * gamma - 1, j - k);
            }
            return r;
        };

        // number of extra terms
        int K = rough ? 128 : 2 * n;

        // sum_{j=0}^{n+K} g_j z^j
        i_t sum_g = 0;
        for (int j = 0; j <= n + K; j++) {
            sum_g += g(j) * pow(z, j);
        }

        // sum_{j=n+1}^{n+K} |c_j| z^j
        i_t sum_c_abs = 0;
        for (int j = n + 1; j <= n + K; j++) {
            sum_c_abs += abs(c_val(j)) * pow(z, j);
        }

        i_t M = pow(abel.r, n + 1) * (sum_c_abs + G_z - sum_g);

        i_t I = sqrt(pi_) / 2
            * bmp::tgamma(i_t(n) / 2) / bmp::tgamma(i_t(n + 1) / 2);

        abel.C0 = 2 * M * I / (gamma_ * b);
        abel.C0 = bmp::upper(abel.C0);
    }

    // find a "good" value of r, for which  C0 / (r-1)^n is small
    if (rough) {
        abel.r_good = abel.r;
    } else {
        i_t accuracy = pow(i_t(10), - NEED_DIGITS_ - 2);
        abel.r_good = 1 + bmp::upper(pow(abel.C0 / accuracy, i_t(1) / i_t(n)));

        abel.r_good = max(abel.r_good, abel.r + r_t(1));
    }

    // now we can call abel_t_sum(*, abel.coef),
    // although the constant term is still zero

    const i_t &am1 = abel.coef(0);
    assert(am1 > 0);

    {   // r1, C1
        // increase r1 until C1 is significantly smaller than a_{-1}
        assert(am1 > 0);
        i_t max_C1 = am1 / (rough ? 4 : 16);
        max_C1 = bmp::lower(max_C1);
        for (abel.r1 = abel.r_good + 1; ; abel.r1 += 1) {
            abel.r1 = bmp::upper(abel.r1); // to be safe
            i_t iC1 = abs(abel.coef(1)) / abel.r1 + abel.C0 / pow(abel.r1 - 1, n);
            for (int k = 1; k <= n; k++)
                iC1 += k * abs(abel.coef(2 + k)) / pow(abel.r1, k + 1);

            if (bmp::upper(iC1) < max_C1) {
                abel.C1 = bmp::upper(iC1);
                break;
            }
        }

        i_t C1am1 = abel.C1 / am1;
        abel.varkappa0 = (am1 - abel.C1) / (am1 + abel.C1);
        abel.varkappa0 = bmp::lower(abel.varkappa0);

        abel.am1_minus_C1 = bmp::lower(i_t(
                    am1 - abel.C1
                    ));
    }

    {   // C2
        abel.C2 = 2 * abel.C1 + abs(abel.coef(1)) / abel.r1;
        for (int k = 1; k <= n; k++) {
            abel.C2 += k * (k + 1) * abs(abel.coef(2 + k))
                / pow(abel.r1, k + 1);
        }
        abel.C2 += 2 * abel.C0 * abel.r1 / pow(abel.r1 - 1, n);

        abel.C2 = bmp::upper(abel.C2);
    }

    // L, nu, M
    if (!rough) {
        assert(abel.r1 > 0 && rho_C_ > 0 && abel.varkappa0 > 0 && N_ > 0);

        // exponential factor which should be outweighed by
        // L / (pi e nu varkappa_1)
        i_t fatty = pow(i_t(10), NEED_DIGITS_ + 2);

        i_t L_nu_factor = 24 + PREC_ / 16;
        i_t L =  log(fatty) / (2 * log(L_nu_factor));
        i_t nu = L_nu_factor * L / (e_ * pi_ * abel.varkappa0);

        abel.L = int(ceil(bmp::upper(L)));
        abel.nu = bmp::upper(nu);

        // M is the number of points on the circle
        abel.halfM = int(ceil(bmp::upper(i_t(
                        log(fatty) / 2
                        ))));
        if (abel.halfM < abel.L)
            abel.halfM = abel.L;
        abel.M = 2 * abel.halfM;
    }

    // Nstar
    if (!rough) {
        // (instead of Caroline's Nstar = int(ceil(nu + mlogeps_)),
        // we slowly increase Nstar until the constraints are satisfied)

        // We need t \geq r_1 + \nu / (a_{-1} + C_1) for the EM
        // approximation to apply
        i_t t_good_for_EM = abel.r1 + abel.nu / (am1 + abel.C1);
        t_good_for_EM = bmp::upper(t_good_for_EM);

        // max error factor in C_\psi, we just choose
        abel.max_kappa_ratio = 8;

        // compute the minimal t for which this factor works,
        i_t t_good_for_kappa;
        {
            i_t C = abel.varkappa0 * abel.nu / abel.am1_minus_C1,
                D = 1 + 2 * abel.C1 / abel.am1_minus_C1,
                E = pow(abel.max_kappa_ratio / D, gamma_ / (gamma_ + 1));
            assert(E > 1);

            // want t / (t - C) <= E
            t_good_for_kappa = C * E / (E - 1);
            t_good_for_kappa = bmp::upper(t_good_for_kappa);
        }

        i_t t_good = max(max(abel.r_good, t_good_for_kappa), t_good_for_EM);

        i_t Ar1 = abel_t_sum(abel.r1, abel.coef)(0);
        i_t x = 1, t = 1;

        for (abel.Nstar = 1; ; abel.Nstar++) {
            x = left_inv(x)(0);
            t = pow(x, -gamma_);

            if (bmp::lower(t) > t_good) {
                r_t min_At = bmp::upper(i_t(
                        Ar1 + abel.nu + 1
                        ));
                r_t At = bmp::lower(abel_t_sum(t, abel.coef)(0));
                if (At > min_At)
                    break;
            }
        }
        abel.x_Nstar = bmp::lower(x);

        // Set the constant term so that A(1) is approximately 0,
        // not trying to control the accuracy...
        // Since we should have A(t) = Nstar:
        abel.coef(2) = 0;
        abel.coef(2) = abel.Nstar - bmp::median(abel_t_sum(t, abel.coef)(0));
    }

    // populate the non-interval coefficients
    for (int k = 0; k <= abel.n + 2; k++)
        abel.coef_ni(k) = bmp::median(abel.coef(k));

    // bound the distance from the preimage of the petal (\cP_1)^{1 / \gamma} under
    // the second branch to the boundary of ellipse A
    abel.dist_P1_to_A = (rho_A_ + 1 / rho_A_) / 8 - i_t(1) / 4 - pow(abel.r1, - gamma_inv_) / 2;
    abel.dist_P1_to_A = bmp::lower(abel.dist_P1_to_A);
    assert(abel.dist_P1_to_A > 0);

    // sanity checks
    assert( abel.r > 0 );
    assert( abel.C0 > 0 );
    assert( abel.r_good > 0 );
    assert( abel.r_good >= abel.r );
    assert( abel.r1 > abel.r_good );
    assert( abel.C1 > 0 );
    assert( abel.C2 > 0 );
    assert( abel.varkappa0 > 0 );
    assert( abel.am1_minus_C1 > 0 );
    if (!rough) {
        assert( abel.nu > 0 );
        assert( abel.x_Nstar > 0 );
        assert( abel.L > 0 && abel.halfM >= abel.L
                && abel.M == 2 * abel.halfM);
        assert(2 * pi_ * e_ * abel.nu * abel.varkappa0 > 2 * abel.L - 1);
    }

    return abel;
}

// precompute the constant part of the error in cheb_sum
template <int PREC>
void LSV<PREC>::compute_sum_small_const_error() {
    const interval_t
        &nu = abel_.nu,
        &vk = abel_.varkappa0;

    interval_t Lfac = 1;
    for (int ell = 1; ell <= 2 * abel_.L + 1; ell++)
        Lfac *= ell;

    // R_L term without C_\psi
    interval_t R_L = Lfac * nu
        / (abel_.L * pow(2 * pi_ * nu * vk, 2 * abel_.L + 1));

    // small observables
    {
        interval_t &c_eps = eps_sum_small_const_error_;
        c_eps = R_L;

        const VectorXi B2 = bernoulli2k(abel_.L);    // B2(ell) = B_{2 ell}

        for (int j = 1; j <= abel_.L; j++)
            c_eps += abs(B2(j)) / (2 * j * pow(nu * vk, 2 * j - 1));

        // kappa factor of C_\psi
        c_eps *= abel_.max_kappa_ratio;
    }

    // Chebyshev polynomials
    interval_t E = R_L;

    // error in D_L without C_\psi
    E += pi_ * pi_ * e_ * nu * vk /
        ( 3 *
          (pow((2 * pi_ * e_ * nu * vk) / (2 * abel_.L - 1), 2) - 1) *
          (exp(interval_t(abel_.M)) - 1)
        );

    // kappa factor of C_\psi
    E *= abel_.max_kappa_ratio;

    // maxima of Chebyshev polynomials in the petal Re(t) > r1
    VectorXi cheb_max_r1(N_);
    {
        interval_t rr = 2 * pow(abel_.r1, - gamma_inv_),
                   tt = 1 + rr + sqrt(pow(rr + 1, 2) - 1);
        for (int i = 0; i < N_; i++)
            cheb_max_r1(i) = (pow(tt, i) + pow(tt, -i)) / 2;
    }

    VectorXi EE = E * cheb_max_r1;

    for (int i = 0; i < N_; i++) {
        assert(EE(i) > 0);
        real_t ei = bmp::upper(EE(i));
        EE(i) = interval_t(-ei, ei);
    }

    cheb_sum_small_const_error_ = EE;

    return;
}

template <int PREC>
LSV<PREC>::MatrixXi LSV<PREC>::abel_matrix(int n) const {
    // X is a KAbel+1 by KAbel+1 matrix where the columns contain coefficients at 1, 1/t, 1/t^2, ... for:
    // * first column: coefficients of t(left(z))
    // * second column: coefficients of log(t(left(z)))
    // * 2+n column: coefficients of t(left(z))^{-k} - t^{-k}
    MatrixXi X = MatrixXi::Zero(n + 2, n + 2);
    const interval_t
        &g = gamma_,
        b = pow(2, g);
    interval_t c;

    c = - g * b;
    for (int j = 0; j <= n + 1; j++) {
        X(j, 0) = c;
        c *= b * (-g - (j + 1)) / (j + 2);
    }
    c = - g * b;
    for (int j = 1; j <= n + 1; j++) {
        X(j, 1) = c / j;
        c *= -b;
    }
    for (int k = 1; k <= n; k++) {
        c = 1;
        for (int j = k + 1; j <= n + 1; j++) {
            c *= (k * g - (j - k - 1)) * b / (j - k);
            X(j, k + 1) = c;
        }
    }
    return X;
}

template <int PREC>
const LSV<PREC>::interval_t LSV<PREC>::e_ = exp(interval_t(1));

template <int PREC>
const LSV<PREC>::interval_t LSV<PREC>::pi_ = 4 * atan(interval_t(1));

template <int PREC>
const LSV<PREC>::real_t LSV<PREC>::real_eps_ = std::numeric_limits<real_t>::epsilon();

template <int PREC>
const LSV<PREC>::real_t LSV<PREC>::real_eps_sqrt_ = sqrt(std::numeric_limits<real_t>::epsilon());

}  // namespace lsv_ns

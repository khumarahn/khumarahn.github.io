#pragma once

#include <limits>
#include <vector>
#include <complex>
#include <cmath> // for pow, abs, log, ... for standard types
#include <type_traits>
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
        // DIGITS is decimal working precision, it should be finer than the target precision
        static constexpr int DIGITS = PREC / 3; // PREC / 3;
        using real_t     = bmp::number<bmp::mpfr_float_backend<DIGITS>>;
        using interval_t = bmp::number<bmp::mpfi_float_backend<DIGITS>>;

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

    private:
        static const interval_t e_;
        static const interval_t pi_;
        static const real_t real_eps_;
        static const real_t real_eps_sqrt_;

        interval_t gamma_;
        VectorXi abel_coef_;
        VectorXr  abel_coef_ni_;

        MatrixXi abel_matrix() const;
        void compute_abel_stuff();

        static constexpr int PREC_ = PREC;
        interval_t mlogeps_;
        int N_, Nstar_, L_, halfM_, M_, Mhat_;

        // Constants for the Abel function, with \tA = \tilde{A} and \hA = \hat{A},
        // \tA_n(t) = a_{-1} + a_\ell + a_0 + a_1 t + ... + a_n t^n
        // * n:
        int abel_n_;
        // * |\tA(z) - \tA_n(z)| \leq C0 |t|^{-n}  when  Re(t) \geq r
        real_t abel_r_, abel_C0_;
        // * \tA(t) is computed with desired accuracy when Re(t) \geq r_good
        real_t abel_r_good_;
        // * |\tA'(t) - a_{-1}| \leq C1 when Re(t) \geq r1,
        //   with C1 significantly smaller than a_{-1}
        real_t abel_r1_, abel_C1_, abel_am1_minus_C1_;
        // * min distance from the Nstar-th preimage of 1 in the t-plane to r1
        real_t abel_nu_;
        // * misc
        real_t abel_delta1_,        // = C1 / sqrt(a_{-1}^2 - C1^2)
               abel_varkappa1_;     // = (1 + delta1^2)^{-1/2}  or  (1 - C_1^2 / a_{-1}^2)^{1/2}

        VectorXci derivatives_s_,
                  derivatives_c_;

        void compute_derivatives_cs() {
            // check that required constants are set
            assert(halfM_ > 0 && M_ > 0 && Nstar_ > 0 && L_ > 0);

            VectorXci &s = derivatives_s_;
            VectorXci &c = derivatives_c_;

            s.resize(halfM_);
            c.resize(halfM_);

            interval_t tau = abel_varkappa1_ * abel_nu_ / e_;

            const VectorXi B2 = bernoulli2k(L_);    // B2(ell) = B_{2 ell}

#pragma omp parallel for schedule(dynamic)
            for (int m = 1; m <= halfM_; m++) {
                interval_t q = 2 * pi_ * (2 * m - 1) / interval_t(2 * M_);
                complex_interval_t qq (0, q);
                s(m - 1) = tau * exp(qq);
                c(m - 1) = 0;
                for (int ell = 1; ell <= L_; ell++) {
                    interval_t w = - q * (2 * ell - 1);
                    complex_interval_t ww (0, w);
                    c(m - 1) += complex_interval_t(B2(ell))
                        * complex_interval_t(pow(tau, -(2*ell - 1)))
                        * exp(ww)
                        / complex_interval_t(2 * ell);
                }
            }
            return;
        }

    public:
        // DEBUG!
        real_t uncertainty(const MatrixXi &X) const {
            real_t a = 0;
            for (const auto &x : X.reshaped())
                a += bmp::width(x);
            return a;
        }
        real_t uncertainty_c(const MatrixXci &X) const {
            real_t a = 0;
            for (const auto &x : X.reshaped())
                a += bmp::width(x.real()) + bmp::width(x.imag());
            return a;
        }

        LSV(const interval_t &gamma = 1) {
            set_gamma(gamma);
        }
        void set_gamma(const interval_t &gamma) {
            gamma_ = gamma;

            // WHY?
            interval_t Rad = interval_t(80) / interval_t(100);

            mlogeps_ = log(interval_t(2)) * PREC_;
            N_ = int(ceil(mlogeps_ / Rad));

            abel_n_ = int(ceil(mlogeps_)) - 1;
            std::cout << "Computing Abel function coeffs and constants...\n";
            std::cout << "abel_n_: " << abel_n_ << "\n";
            compute_abel_stuff();
            std::cout
                << "Nstar_: " << Nstar_ << "\n"
                << "abel_r_: " << abel_r_ << ", "
                << "abel_C0_: " << abel_C0_ << ", "
                << "abel_r_good_: " << abel_r_good_ << "\n"
                << "abel_r1_: " << abel_r1_ << ", "
                << "abel_C1_: " << abel_C1_ << "\n"
                << "abel_delta1_: " << abel_delta1_ << ", "
                << "abel_varkappa1_: " << abel_varkappa1_ << "\n"
                << "abel_am1_minus_C1_: " << abel_am1_minus_C1_ << "\n"
                << "abel_nu_: " << abel_nu_ << "\n"
                << "First Abel coefficients: " << abel_coef_ni_.head(5).transpose()
                << "\n\n";

            // L is the number of derivatives FIXME: what should it be?
            // It used to be
            //      L_ = 1 + Nstar_ / 2;
            // but we need it in the computation of Nstar_, so now we set it there

            // M is the number of points on the circle
            halfM_ = L_;
            M_ = 2 * halfM_;

            assert(M_ >= 2*L_);

            Mhat_ = 2 * M_;

            compute_derivatives_cs();
        }

        int NCheb() const { return N_; };
        int prec() const { return PREC_; };
        int KAbel() const { return abel_n_ + 1; };
        int abel_n() const { return abel_n_; };
        interval_t alpha() const { return 1 / gamma_; };
        interval_t gamma() const { return gamma_; };
        VectorXi abel_coef() const { return abel_coef_; };

        // branches
        Vector2i left(const interval_t &x) const {
            interval_t t = pow(2 * x, gamma_);
            return Vector2i(x * (1 + t), 1 + (gamma_ + 1) * t);
        }
        Vector2i right(const interval_t &x) const {
            return Vector2i(2 * x - 1, 2);
        }

        // full map
        Vector2i map(const interval_t &x) const {
            return (x < 0.5) ? left(x) : right(x);
        }

        // inverses
        Vector2i left_inv(const interval_t &x) const {
            // max derivative is
            interval_t md = left(interval_t(1) / 2)(1);
            // an interval enclosing the root
            interval_t guess(x / md, x);

            auto f = [this, &x] (const interval_t &z) {
                Vector2<interval_t> r = left(z);
                r(0) -= x;
                return r;
            };
            return interval_root_ns::interval_newton(f, guess);
        }

        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> right_inv(const var_t &x) const {
            return Vector2<var_t>((x + var_t(1)) / var_t(2), 0.5);
        }

        // Transfer operator for the induced map, the clever one.
        // This returns an approximation of the induced transfer operator by an
        // N_ by N_ matrix acting on the Chebyshev coefficients.
        MatrixXi Lind() const {
            const interval_cheb_t cheb(interval_t(1) / 2, 1, N_);
            const VectorXi x_nodes = cheb.nodes();

            // values of first N Chebyshev polynomials at the preimage on [1/2,1]
            auto phi = [this, &cheb]<typename var_t>(const var_t &x) {
                var_t y = this->right_inv(x)(0);
                // slow trigonometric evaluation is more precise
                // in interval arithmetic
                VectorX<var_t> r = cheb.basis_values_trig(y, N_);
                return r;
            };

            // constant part of the error in S_small, except for multiplication
            // by |A'(x)|.
            const VectorXi S_small_error = [this] () {
                const interval_t &nu = abel_nu_;
                const real_t &vk = abel_varkappa1_;

                interval_t Lfac = 1;
                for (int ell = 1; ell <= 2 * L_ + 1; ell++)
                    Lfac *= ell;

                interval_t E;

                E = Lfac * nu / (L_ * pow(2 * pi_ * nu * vk, 2 * L_ + 1));

                E += pi_ * pi_ * e_ * nu * vk /
                    ( 6 *
                      (pow((2 * pi_ * e_ * nu * vk) / (2 * L_ - 1), 2) - 1) *
                      (exp(interval_t(M_)) - 1)
                    );

                E /= (gamma_ * abel_am1_minus_C1_ * pow(abel_r1_, 1 + 1 / gamma_));

                // maxima of Chebyshev polynomials in the petal Re(t) > r1
                VectorXi cheb_max_r1(N_);
                {
                    interval_t rr = 2 * pow(abel_r1_, - 1 / gamma_),
                               tt = 1 + rr + sqrt(rr * rr + 2 * rr);
                    for (int i = 0; i < N_; i++)
                        cheb_max_r1(i) = (pow(tt, i) + pow(tt, -i)) / 2;
                }

                VectorXi EE = E * cheb_max_r1;

                for (int i = 0; i < N_; i++) {
                    assert(EE(i) > 0);
                    real_t ei = bmp::upper(EE(i));
                    EE(i) = interval_t(-ei, ei);
                }
                return EE;
            }();


            // for a **small** real x, return a rigorous approximation of S(z) = S(x^\gamma)
            // for the vector-valued observable of first N Chebyshev polynomials
            auto S_small = [this, &cheb, &phi, &S_small_error](const interval_t &x) -> VectorXi {
                //interval_t z = pow(x, gamma_);
                VectorXi r = VectorXi::Zero(N_);
                Vector2i Ax = this->abel(x);

                {   // integral
                    interval_t y = this->right_inv(x)(0);
                    VectorXi bi = cheb.beta_integral(0, y, N_);
                    r += - 2 * Ax(1) * bi;
                }

                // \varphi(z) / 2
                r += phi(x) / real_t(2);

                {   // derivatives
                    const VectorXci &s = derivatives_s_;
                    const VectorXci &c = derivatives_c_;

                    VectorXi der = VectorXi::Zero(N_);
                    for (int m = 1; m <= halfM_; m++) {
                        Vector2ci xm_inv = abel_inv(Ax(0) + s(m-1));
                        der += ( c(m-1) * phi(xm_inv(0)) * xm_inv(1) ).real();
                    }

                    der *= - 2 * Ax(1) / M_;
                    r += der;
                }

                // error
                r += abs(Ax(1)) * S_small_error;

                return r;
            };

            // now compute S when x is not necessarily small
            auto S = [this, &S_small, &phi](const interval_t &x) -> VectorXi const {
                VectorXi r = VectorXi::Zero(N_);

                interval_t xk = x,
                           inv_Jk = 1;

                // add branches naively up to Nstar
                for (int k = 0; k < Nstar_; k++) {
                    r += phi(xk) * inv_Jk;
                    Vector2i y = left_inv(xk);
                    xk = y(0);
                    inv_Jk *= y(1);
                }

                // add the rest when xk is small
                r += S_small(xk) * inv_Jk;
                return r;
            };

            // Compute (L T_n)(x) for n=0,...,N-1 and x in the node points
            MatrixXi L_values(N_, N_);
#pragma omp parallel for schedule(dynamic)
            for (int ix = 0; ix < x_nodes.size(); ix++) {
                interval_t x = x_nodes[ix];
                L_values.col(ix) = S(x) / interval_t(2);
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

            return R;
        }

        interval_cheb_t h() const {
            interval_t a(interval_t(1) / 2),
                       b(1);

            MatrixXi L = Lind();
            int LN = L.cols();

            VectorXi iota(LN);
            VectorXi u = VectorXi::Zero(LN);

            {   // compute values of integrals on [1/2,1] of the basis Chebyshev polynomials
                interval_cheb_t ii(a, b, 1);
                iota = ii.beta_integral(0, b, LN);
            }

            u(0) = 1 / iota(0);

            // Now the invariant density in Chebyshev basis
            // is h = (I - L + u iota)^{-1} u
            MatrixXi S = MatrixXi::Identity(LN,LN) - L + u * iota.transpose();

            VectorXi hv = interval_root_ns::linear_krawczyk(S, u);

            // h as a function
            interval_cheb_t hh(hv, a, b);

            return hh;
        }

        // *** ABEL FUNCTIONS ***
        // * sum for \tA(t), or abel_t, with a provided array of coefficients, without the error term
        template <typename var_t, typename coef_t> requires
            (type_one_of<var_t, real_t, complex_t> && type_one_of<coef_t, VectorXr>) ||
            (type_one_of<var_t, interval_t, complex_interval_t> && type_one_of<coef_t, VectorXi>)
        Vector2<var_t> abel_t_sum(const var_t &t, const coef_t &coef) const {
            assert(coef.size() == abel_n_ + 3);
            var_t A = coef(0) * t + coef(1) * log(t) + coef(2),
                  dA = coef(0) + coef(1) / t;

            const var_t ti = var_t(1) / t;

            var_t A_tail(0), dA_tail(0);

            // Horner's method
            for (int j = abel_n_; j >= 1; j--) {
                A_tail = coef(2 + j) + ti * A_tail;
                dA_tail = var_t(j) * coef(2 + j) + ti * dA_tail;
            }
            A += ti * A_tail;
            dA -= ti * ti * dA_tail;

            return Vector2<var_t>(A, dA);
        }
        // * non-interval version
        template <typename var_t> requires type_one_of<var_t, real_t, complex_t>
        Vector2<var_t> abel_t(const var_t &t) const {
            return abel_t_sum(t, abel_coef_ni_);
        }
        // * interval version
        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> abel_t(const var_t &t) const {
            return abel_t_sum(t, abel_coef_) + abel_t_error(t);
        }
        // * interval A(x)
        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> abel(const var_t &x) const {
            const var_t gamma = var_t(gamma_);
            var_t t = pow(x, -gamma);
            Vector2<var_t> A = abel_t(t);
            return Vector2<var_t>(A(0), -gamma * (t / x) * A(1));
        }

        // inverse
        // * non-interval version in t, Newton method
        template <typename var_t> requires type_one_of<var_t, real_t, complex_t>
        Vector2<var_t> abel_t_inv(const var_t &a) const {
            var_t t = a / abel_coef_ni_(0);
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
                return complex_interval_t(median(x.real()), median(x.imag()));
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
            const var_t & t = Ai(0), dAt = Ai(1);
            var_t x = pow(t, var_t(-1 / gamma_));
            return Vector2<var_t>(x, -(x / t) / gamma_ * dAt);
        }

        // error term for abel_t
        template <typename i_t> requires type_one_of<i_t, interval_t, complex_interval_t>
        Vector2<i_t> abel_t_error(const i_t &t) const {

            // check is we are in the right region where we can compute the
            // Abel function accurately
            if constexpr (type_one_of<i_t, complex_interval_t>)
                assert( bmp::lower(t.real()) >= abel_r_good_ );
            else
                assert( bmp::lower(t) >= abel_r_good_ );

            // ir_t is the underlying real interval type
            using ir_t = std::conditional<
                type_one_of<i_t,complex_interval_t>,  interval_t,  i_t
                >::type;

            ir_t factor0 = abel_C0_ * pow(abs(t), -abel_n_),
                 window0 (-bmp::upper(factor0), bmp::upper(factor0)),
                 factor1 = abel_C0_ * pow(abs(t - i_t(1)), -abel_n_),
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
                    B2[m / 2] = interval_t(A(0)) / interval_t(fact);
            }
            return B2;
        }

};  // class LSV

template <int PREC>
void LSV<PREC>::compute_abel_stuff() {

    using r_t = real_t;
    using i_t = interval_t;

    const int n = abel_n_;
    assert(n > 0);

    abel_coef_.resize(n + 3);
    abel_coef_ni_.resize(n + 3);

    {   // first compute non-constant coefficients (am1, al, a1, a2, ...) of the Abel function
        VectorXi b = VectorXi::Zero(abel_n_ + 2);
        b(0) = -1;

        MatrixXi A = abel_matrix(); // lower triangular

        // brute force solve lower triangular matrix
        // instead of using Eigen that can fail in interval arithmetic:
        //VectorXi x = abel_matrix().template triangularView<Eigen::Lower>().solve(b);
        //VectorXi x = interval_root_ns::linear_krawczyk(A, b);

        VectorXi xx(abel_n_ + 2);
        for (int k = 0; k < n + 2; k++) {
            i_t s = 0;
            for (int j = 0; j < k; j++) {
                s += A(k,j) * xx(j);
            }
            xx(k) = (b(k) - s) / A(k,k);
        }

        abel_coef_(0) = xx(0);
        abel_coef_(1) = xx(1);
        abel_coef_(2) = 0;
        for (int j = 2; j <= n + 1; j++) {
            abel_coef_(j + 1) = xx(j);
        }
    }

    {   // compute r, C0
        i_t q = i_t(1) / 2,
            b = pow(2, gamma_);

        i_t r_m1 = b / q,
            r_m2 = 2 * b * (gamma_ + 1) * (pow(1 - q, -gamma_) - 1) / (gamma_ * q);
        abel_r_ = max(bmp::upper(r_m1), bmp::upper(r_m2));

        i_t R = (1 / abel_r_ + 1 / b) / 2,
            bR = b * R;

        i_t B = 1
            + abs(abel_coef_(0)) / R * (pow(1 - bR, -gamma_) + 1)
            - abs(abel_coef_(1)) * gamma_ * log(1 - bR);
        for (int k = 1; k <= n; k++) {
            B += abs(abel_coef_(2 + k)) * pow(R, k) * (pow(1 + bR, k * gamma_) + 1);
        }

        i_t M = B / pow(R, n + 1);

        i_t I = sqrt(pi_) / 2
            * bmp::tgamma(i_t(n) / 2) / bmp::tgamma(i_t(n + 1) / 2);

        abel_C0_ = bmp::upper(pow(i_t(2), n + 1) * M * I / (gamma_ * b));
    }

    {   // compute a "good" value of r, for which
        // C0 / (r-1)^n is smaller than desired precision
        i_t accuracy = pow(i_t(2), - PREC_);
        abel_r_good_ = 1 + bmp::upper(pow(abel_C0_ / accuracy, i_t(1) / i_t(n)));

        abel_r_good_ = max(abel_r_good_, abel_r_ + r_t(1));
    }

    // now we can call abel_t, although the constant term is still zero

    {   // r1, C1

        // increase r1 until C1 is significantly smaller than a_{-1}
        i_t max_C1 = bmp::lower(i_t(abel_coef_(0) / 16));
        for (abel_r1_ = abel_r_good_ + 1; ; abel_r1_ += 1) {
            i_t iC1 = abs(abel_coef_(1)) / abel_r1_ + abel_C0_ / pow(abel_r1_ - 1, n);
            for (int k = 1; k <= n; k++)
                iC1 += k * abs(abel_coef_(2 + k)) / pow(abel_r1_, k + 1);

            if (bmp::upper(iC1) < max_C1) {
                abel_C1_ = bmp::upper(iC1);
                break;
            }
        }

        i_t C1am1 = abel_C1_ / abel_coef_(0);
        abel_delta1_    = bmp::upper(i_t(pow(pow(C1am1, -2) - 1, i_t(-1) / 2)));
        abel_varkappa1_ = bmp::lower(i_t(pow(1 - pow(C1am1, 2), i_t(1) / 2)));
    }

    i_t t_Nstar;
    {   // Nstar and L
        // (instead of Caroline's Nstar_ = int(ceil(nu + mlogeps_)),
        // we slowly increase Nstar until the constraints are satisfied)
        i_t Ar1 = abel_t(i_t(abel_r1_))(0);
        i_t x = 1, t = 1;

        for (Nstar_ = 1; ; Nstar_++) {
            x = left_inv(x)(0);
            t = pow(x, -gamma_);

            // unremarkable, but we set the number of derivatives L here
            L_ = 1 + Nstar_ / 2;

            if (bmp::lower(t) > abel_r_good_) {
                // +20 is arbitrary: we need \nu suffuciently large
                r_t min_At = 20 + bmp::upper(
                        Ar1 + i_t(L_) / (e_ * pi_ * abel_varkappa1_)
                        );
                r_t At = bmp::lower(abel_t(t)(0));
                if (At > min_At) {
                    t_Nstar = t;
                    break;
                }
            }
        }
    }

    {   // Now set the constant term so that A(1) is approximately 0.
        // Not trying to control the accuracy
        i_t x = 1;
        for (int i=0; i<Nstar_; i++) {
            x = left_inv(x)(0);
        }

        i_t t = pow(x, -gamma_);

        // now, we should have A(t) = Nstar_
        abel_coef_(2) = 0;
        abel_coef_(2) = Nstar_ - bmp::median(abel_t(t)(0));
    }

    // populate the non-interval coefficients
    for (int k = 0; k <= abel_n_ + 2; k++)
        abel_coef_ni_(k) = bmp::median(abel_coef_(k));

    abel_am1_minus_C1_ = bmp::lower(i_t(
                abel_coef_(0) - abel_C1_
                ));
    abel_nu_ = bmp::lower(i_t(
                abel_t(t_Nstar)(0) - abel_t(i_t(abel_r1_))(0)
                ));

    // sanity checks
    assert( abel_r_ > 0 );
    assert( abel_C0_ > 0 );
    assert( abel_r_good_ > 0 );
    assert( abel_r_good_ >= abel_r_ );
    assert( abel_r1_ > abel_r_good_ );
    assert( abel_C1_ > 0 );
    assert( abel_delta1_ > 0 );
    assert( abel_varkappa1_ > 0 );
    assert( abel_am1_minus_C1_ > 0 );
    assert( abel_nu_ > 0 );

    assert(2 * pi_ * e_ * abel_nu_ * abel_varkappa1_ > 2 * L_ - 1);

    return;
}

template <int PREC>
LSV<PREC>::MatrixXi LSV<PREC>::abel_matrix() const {
    // X is a KAbel+1 by KAbel+1 matrix where the columns contain coefficients at 1, 1/t, 1/t^2, ... for:
    // * first column: coefficients of t(left(z))
    // * second column: coefficients of log(t(left(z)))
    // * 2+n column: coefficients of t(left(z))^{-k} - t^{-k}
    MatrixXi X = MatrixXi::Zero(abel_n_ + 2, abel_n_ + 2);
    const interval_t g = gamma_, b = pow(2, g);
    interval_t c;

    c = - g * b;
    for (int j = 0; j <= abel_n_ + 1; j++) {
        X(j,0) = c;
        c *= b * (-g - j - 1) / (j + 2);
    }
    c = - g * b;
    for (int j = 1; j <= abel_n_ + 1; j++) {
        X(j,1) = c / j;
        c *= -b;
    }
    for (int k = 1; k <= abel_n_; k++) {
        c = 1;
        for (int j = k + 1; j <= abel_n_ + 1; j++) {
            c *= (k * g - j + k + 1) * b / (j - k);
            X(j, k+1) = c;
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

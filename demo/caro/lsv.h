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

#include "clenshawcurtis.h"
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

template <int PREC = 16>
class LSV {
    public:
        // PREC is binary target precision;
        // DIGITS is decimal working precision, it should be finer than the target precision
        static constexpr int DIGITS = PREC / 2; // PREC / 3;
        using real_t     = bmp::number<bmp::mpfr_float_backend<DIGITS>>;
        using interval_t = bmp::number<bmp::mpfi_float_backend<DIGITS>>;

        // extra precision
        using real_extra_t     = bmp::number<bmp::mpfr_float_backend<2 * DIGITS>>;
        using interval_extra_t = bmp::number<bmp::mpfi_float_backend<2 * DIGITS>>;

        using complex_t = std::complex<real_t>;
        using complex_interval_t = std::complex<interval_t>;

        template <typename var_t>
            using Vector2 = Eigen::Matrix<var_t,2,1>;
        template <typename var_t>
            using VectorX = Eigen::Matrix<var_t,Eigen::Dynamic,1>;
        template <typename var_t>
            using MatrixX = Eigen::Matrix<var_t,Eigen::Dynamic,Eigen::Dynamic>;
        template <typename var_t>
            using func_t = std::function<var_t (var_t)>;

        typedef VectorX<real_t> VectorXr;
        typedef Vector2<real_t> Vector2r;
        typedef Vector2<interval_t> Vector2i;
        typedef VectorX<interval_t> VectorXi;
        typedef VectorX<interval_extra_t> VectorXix;

        typedef MatrixX<real_t> MatrixXr;
        typedef MatrixX<interval_t> MatrixXi;
        typedef MatrixX<interval_extra_t> MatrixXix;

        typedef func_t<interval_t> interval_func_t;

        typedef Vector2<complex_t> Vector2c;
        typedef Vector2<complex_interval_t> Vector2ci;
        typedef VectorX<complex_interval_t> VectorXci;

        typedef MatrixX<complex_interval_t> MatrixXci;

        typedef func_t<complex_interval_t> complex_func_t;

        typedef Cheb<interval_t> interval_cheb_t;
    private:
        static const interval_t e_;
        static const interval_t pi_;
        static const interval_extra_t pi_x_;
        static const real_t real_eps_;
        static const real_t real_eps_sqrt_;

        interval_t gamma_;
        VectorXi abel_coef_;
        VectorXr  abel_coef_ni_;

        void compute_abel_stuff();
        MatrixXix abel_matrix() const;

        static constexpr int PREC_ = PREC;
        interval_t mlogeps_;
        int N_, Nstar_, L_, halfM_, M_, Mhat_;

        // Constants from the Abel function, with \tA = \tilde{A} and \hA = \hat{A},
        // \tA_n(t) = a_{-1} + a_\ell + a_0 + a_1 t + ... + a_n t^n
        int abel_n_;
        // |\tA(z) - \tA_n(z)| \leq abel_C0_ |t|^{-n}  when  Re(t) \geq abel_r_
        real_t abel_r_, abel_C0_;
        // |\tA'(t) - a_{-1}| \leq abel_C1_ when Re(t) \leq abel_r1_
        real_t abel_r1_, abel_C1_, abel_am1_minus_C1_;
        // min distance the Nstar preimage of 1 in the t-plane to r1
        real_t abel_nu_;
        // derived variables
        real_t abel_delta1_,        // = C1 / sqrt(a_{-1}^2 - C1^2)
               abel_varkappa1_;     // = (1 + delta1^2)^{-1/2}  or  (1 - C_1^2 / a_{-1}^2)^{1/2}

        template <typename i_t> requires type_one_of<i_t, interval_t, interval_extra_t>
        bool abel_t_in_range(const i_t &t) const {
            return bmp::lower(t) >= abel_r_;
        }
        bool abel_t_in_range(const complex_interval_t &t) const {
            return bmp::lower(t.real()) >= abel_r_;
        }
        template <typename i_t> requires type_one_of<i_t, interval_t, interval_extra_t, complex_interval_t>
        i_t abel_t_error(const i_t &t) const {

            // ir_t is the underlying real interval type
            using ir_t = std::conditional<type_one_of<i_t,complex_interval_t>,
                  interval_t,
                  i_t>::type;

            if (!abel_t_in_range(t)) {
                assert(false);
            }

            const int n = abel_n_;
            ir_t factor = abel_C0_ * pow(abs(t), -n),
                 window (-bmp::upper(factor), bmp::upper(factor));

            //DEBUG
            {
                if (factor > 1) {
                    std::cout << "PIZDA!! t: " << t << ", factor: " << factor << "\n";
                    std::cout << "n: " << n << ", abel_C0_: " << abel_C0_ << "\n";
                }
            }

            if constexpr (type_one_of<i_t, complex_interval_t>) {
                return complex_interval_t(window, window);
            } else {
                return window;
            }
        }

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
            // DEBUG debug commented temporary temporarily FIXME:
            //set_gamma(gamma);
        }
        void set_gamma(const interval_t &gamma) {
            gamma_ = gamma;

            // WHY?
            interval_t Rad = interval_t(80) / interval_t(100);

            mlogeps_ = log(interval_t(2)) * PREC_;
            N_ = int(ceil(mlogeps_ / Rad));

            std::cout << "Computing Abel function coeffs and constants...\n";
            abel_n_ = int(ceil(mlogeps_)) - 1;
            compute_abel_stuff();
            std::cout << "   ... done...\n"
                << "abel_n_: " << abel_n_ << ", "
                << "Nstar_: " << Nstar_ << "\n"
                << "abel_r_: " << abel_r_ << ", "
                << "abel_C0_: " << abel_C0_ << "\n"
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
        template <typename i_t> requires type_one_of<i_t, interval_t, interval_extra_t>
        Vector2<i_t> left(const i_t &x) const {
            i_t t = pow(2 * x, gamma_);
            return Vector2<i_t>(x * (1 + t), 1 + (gamma_ + 1) * t);
        }
        Vector2i right(const interval_t &x) const {
            return Vector2i(2 * x - 1, 2);
        }

        // full map
        Vector2i map(const interval_t &x) const {
            return (x < 0.5) ? left(x) : right(x);
        }

        // inverses
        template <typename i_t> requires type_one_of<i_t, interval_t, interval_extra_t>
        Vector2<i_t> left_inv(const i_t &x) const {
            // max derivative is
            i_t md = left(i_t(0.5))(1);
            // an interval enclosing the root
            i_t guess(x / md, x);

            auto f = [this, &x] (const i_t &z) {
                Vector2<i_t> r = left(z);
                r(0) -= x;
                return r;
            };
            return interval_root_ns::interval_newton(f, guess);
        }

        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> right_inv(const var_t &x) const {
            return Vector2<var_t>((x + var_t(1)) / var_t(2), 0.5);
        }

        // transfer operator wrt dx
        interval_t L(const interval_func_t &v, const interval_t &x) const {
            Vector2i y1 = left_inv(x),
                     y2 = right_inv(x);
            return v(y1(0)) * y1(1) + v(y2(0)) * y2(1);
        }
        // transfer operator wrt x^{-gamma} dx
        interval_t Lxgamma (const interval_func_t &v, const interval_t &x) const {
            Vector2i y1 = left_inv(x),
                     y2 = right_inv(x);
            interval_t w1 = y1(1) * (y1(0) <= 0 ? interval_t(1) : pow(x / y1(0), gamma_)),
                   w2 = y2(1) * (y2(0) <= 0 ? interval_t(1) : pow(x / y2(0), gamma_));
            return v(y1(0)) * w1 + v(y2(0)) * w2;
        }

        // transfer operator of induced map with n branches
        interval_t Lind(const interval_func_t &v, const interval_t &x, int n) const {
            interval_t r = 0, z = x, w = 1;
            for (int k = 0; k < n; k++) {
                Vector2i yr = right_inv(z);
                r += v(yr(0)) * w * yr(1);
                Vector2i l = left_inv(z);
                z = l(0);
                w *= l(1);
            }
            return r;
        }

        // transfer operator of induced map, the clever one
        // This returns an approximation of the induced transfer operator by an
        // N_ by N_ matrix acting on the Chebyshev coefficients
        MatrixXi Lind() const {
            interval_cheb_t cheb(0.5, 1.0, N_);

            const VectorXi x_nodes = cheb.nodes();
            assert(x_nodes.size() == N_);

            // values of first N Chebyshev polynomials at a point
            auto v = [&cheb, N = N_]<typename var_t>(const var_t &x) {
                VectorX<var_t> r = cheb.basis_values_trig(x, N);
                r(0) *= 2; // undo halving T_1(x)
                return r;
            };

            // for a **small** real x, return a rigorous approximation of S(z) = S(x^\gamma)
            // for the vector-valued observable of first N Chebyshev polynomials
            auto S = [this, &cheb, &v, N = N_](const interval_t &x) -> VectorXi {
                // FIXME: what happens to the first Chebyshev polynomial, which is either halved or doubled?
                interval_t z = pow(x, gamma_);
                VectorXi r = VectorXi::Zero(N);
                VectorXi Ax = abel(x);

                // integral
                VectorXi bi = cheb.beta_integral(0, (x + 1) / 2, N);
                bi(0) *= 2; // undo halving T_1(x)
                r += - 2 * Ax(1) * bi;

                // \varphi(z) / 2
                r += v(right_inv(x)(0));

                { // derivatives
                    const VectorXci &s = derivatives_s_;
                    const VectorXci &c = derivatives_c_;

                    VectorXi der = VectorXi::Zero(N);
                    for (int m = 1; m <= halfM_; m++) {
                        complex_interval_t xm = abel_inv(Ax(0) + s(m-1))(0);
                        der += ( c(m-1) * v(right_inv(xm)(0)) / abel(xm)(1) ).real();
                    }

                    der *= - 2 * Ax(1) / M_;
                    r += der;
                }

                { // error FIXME: Make this independent of x
                    interval_t nu = Ax(0) - abel_r1_;

                    interval_t Lfac = 1;
                    for (int ell = 1; ell <= 2 * L_ + 1; ell++)
                        Lfac *= ell;

                    interval_t E = 0;
                    const real_t &vk = abel_varkappa1_;

                    E += Lfac * nu / (L_ * pow(2 * pi_ * nu * vk, 2 * L_ + 1));

                    E += pi_ * pi_ * e_ * nu * vk /
                        ( 6 *
                          (pow((2 * pi_ * e_ * nu * vk) / (2 * L_ - 1), 2) - 1) *
                          (exp(interval_t(M_)) - 1)
                        );

                    E *= abs(Ax(1)) / (gamma_ * abel_am1_minus_C1_ * pow(abel_r1_, 1 + 1 / gamma_));

                    // maxima of Chebyshev polynomials in the petal
                    VectorXi TM(N);
                    {
                        interval_t rr = 2 * pow(abel_r1_, - 1 / gamma_),
                                   tt = 1 + rr + sqrt(rr * rr + 2 * rr);
                        for (int i = 0; i < N; i++)
                            TM(i) = (pow(tt, i) + pow(tt, -i)) / 2;
                    }

                    VectorXi EN = E * TM;

                    for (int i = 0; i < N; i++) {
                        assert(EN(i) > 0);
                        real_t ei = bmp::upper(EN(i));
                        r(i) += interval_t(-ei, ei);

                        if (bmp::width(EN(i)) > 0.001 || ei > 0.001)
                            std::cout << "TOLL: " << i
                                << " :: " << bmp::width(EN(i))
                                << " :: " << ei << "\n";
                    }
                }

                return r;
            };

            MatrixXi L_values(N_, N_);
#pragma omp parallel for schedule(dynamic)
            for (int ix = 0; ix < x_nodes.size(); ix++) {
                interval_t x = x_nodes[ix];
                VectorXi r = VectorXi::Zero(N_);

                // add branches naively up to Nstar
                interval_t xk = x, Jk = 1;

                for (int k = 0; k < Nstar_; k++) {
                    std::cout << "TROLL: " << k << " :: "
                        << "xk: " << bmp::width(xk)
                        << ", Jk: " << bmp::width(Jk)
                        << ", v(xk): " << uncertainty(v(xk))
                        << ", r: " << uncertainty(r)
                        << "\n"
                        << "%% ";
                    //for (int t = 0; t < v(xk).size(); t++)
                    //    std::cout << bmp::width(v(xk)(t)) << "  ";
                    //std::cout << "\n";

                    r += v(xk) * Jk;
                    Vector2i y = left_inv(xk);
                    xk = y(0);
                    Jk *= y(1);
                }

                r += S(xk) * Jk;

                L_values.col(ix) = r;
            }

            MatrixXi R(N_, N_);
            {
                interval_cheb_t c(interval_t(1) / 2, 1, N_);
                for (int k = 0; k < N_; k++) {
                    c.set_from_values(L_values.row(k).transpose());
                    VectorXi cc = c.coef();
                    std::cout << "ROCK: " << k
                        << " :: " << uncertainty(L_values.col(k))
                        << "\n";
                    std::cout << "ROLL: " << k
                        << " :: " << uncertainty(L_values.row(k).transpose())
                        << " --vs-- " << uncertainty(cc)
                        << "\n";
                    cc(0) /= 2; // undo halving T_1(x)
                    R.col(k) = cc;
                }
            }

            return R;
        }

        // Approximate Abel function A(x), satisfying A(left(x)) = A(x) - 1.
        // With t = t(x) = x^{-gamma}, and N = KAbel-1, we use an asymptotic approximation
        //      A(x) ≈ am1 t + al log t + a0 + a1 / t + a2 / t^2 + ... + aN / t^N
        // which is valid for sufficiently large Re(t) with a controlled error
        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> abel(const var_t &x) const {
            const var_t gamma = var_t(gamma_);
            var_t t = pow(x, -gamma);
            Vector2<var_t> A = abel_t(t);
            return Vector2<var_t>(A(0), -gamma * (t / x) * A(1));
        }
        template <typename var_t, typename coef_t> requires
            (type_one_of<var_t, real_t, complex_t> && type_one_of<coef_t, VectorXr>) ||
            (type_one_of<var_t, interval_t, complex_interval_t> && type_one_of<coef_t, VectorXi>) ||
            (type_one_of<var_t, interval_extra_t> && type_one_of<coef_t, VectorXix>)
        Vector2<var_t> abel_t(const var_t &t, const coef_t &coef) const {
            assert(coef.size() == abel_n_ + 3);
            var_t A = coef(0) * t + coef(1) * log(t) + coef(2),
                  dA = coef(0) + coef(1) / t;

            const var_t ti = var_t(1) / t;
            var_t tj = ti;
            for (int j = 1; j <= abel_n_; j++) {
                A += coef(2 + j) * tj;
                tj *= ti;
                dA -= var_t(j) * coef(2 + j) * tj;
            }
            if constexpr (type_one_of<var_t, interval_t, interval_extra_t, complex_interval_t>) {
                A += abel_t_error(t);
            }
            return Vector2<var_t>(A, dA);
        }
        template <typename var_t> requires type_one_of<var_t, real_t, complex_t>
        Vector2<var_t> abel_t(const var_t &t) const {
            return abel_t(t, abel_coef_ni_);
        }
        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> abel_t(const var_t &t) const {
            return abel_t(t, abel_coef_);
        }

        // inverse
        // * non-interval version, Newton method
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
        // * real interval version
        Vector2i abel_t_inv(const interval_t &a) const {
            Vector2r g = abel_t_inv<real_t>(bmp::median(a));
            // an interval enclosing the root
            interval_t guess = g(0) + 16 * g(1) * (bmp::width(a) + 0.0001) * interval_t(-1, 1);

            auto f = [this, &a] (const interval_t &t) {
                Vector2i r = abel_t(t);
                r(0) -= a;
                return r;
            };
            return interval_root_ns::interval_newton(f, guess);
        }
        // * complex interval version
        Vector2ci abel_t_inv(const complex_interval_t &a) const {
            auto m = [] (complex_interval_t x) {
                using bmp::median;
                return complex_interval_t(median(x.real()), median(x.imag()));
            };
            Vector2c g = abel_t_inv<complex_t>(m(a));
            complex_interval_t g0 = g(0),
                               g1 = g(1);
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
        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> abel_inv(const var_t &a) const {
            Vector2<var_t> Ai = abel_t_inv(a);
            const var_t & t = Ai(0), dAt = Ai(1);
            var_t x = pow(t, var_t(-1 / gamma_));
            return Vector2<var_t>(x, -(x / t) / gamma_ * dAt);
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

    using rx_t = real_extra_t;
    using ix_t = interval_extra_t;

    abel_coef_.resize(abel_n_ + 3);
    abel_coef_ni_.resize(abel_n_ + 3);

    rx_t r, C0, r1, C1, delta1, varkappa1;

    // We do the computation in extra precision first, then dumb it down
    VectorXix x_coef(abel_n_ + 3);

    { // first compute non-constant coefficients (am1, al, a1, a2, ...) of the Abel function
        VectorXix b = VectorXix::Zero(abel_n_ + 2);
        b(0) = -1;

        MatrixXix A = abel_matrix(); // lower triangular

        // brute force solve lower triangular matrix
        // instead of using Eigen that can fail in interval arithmetic:
        //VectorXi x = abel_matrix().template triangularView<Eigen::Lower>().solve(b);
        //VectorXi x = interval_root_ns::linear_krawczyk(A, b);

        VectorXix xx(abel_n_ + 2);
        for (int k = 0; k < abel_n_ + 2; k++) {
            ix_t s = 0;
            for (int j = 0; j < k; j++) {
                s += A(k,j) * xx(j);
            }
            xx(k) = (b(k) - s) / A(k,k);
        }

        x_coef(0) = xx(0);
        x_coef(1) = xx(1);
        x_coef(2) = 0;
        for (int j = 2; j <= abel_n_ + 1; j++) {
            x_coef(j + 1) = xx(j);
        }
    }

    { // compute r, C0
        const int n = abel_n_;

        ix_t q = ix_t(1) / 2,
             b = pow(2, gamma_);

        r = max(
                bmp::upper(ix_t(
                        b / q
                        )),
                bmp::upper(ix_t(
                        2 * b * (gamma_ + 1) * (pow(1 - q, -gamma_) - 1) / (gamma_ * q)
                        )
                    ));

        ix_t R = (1 / r + 1 / b) / 2,
             bR = b * R;

        ix_t B = 1
            + abs(x_coef(0)) / R * (pow(1 - bR, -gamma_) + 1)
            - abs(x_coef(1)) * gamma_ * log(1 - bR);
        for (int k = 1; k <= n; k++) {
            B += abs(x_coef(2 + k)) * pow(R, k) * (pow(1 + bR, k * gamma_) + 1);
        }

        ix_t M = B / pow(R, n + 1);

        ix_t I = sqrt(pi_) / 2
            * bmp::tgamma(ix_t(n) / 2) / bmp::tgamma(ix_t(n + 1) / 2);

        C0 = bmp::upper(pow(ix_t(2), n + 1) * M * I / (gamma_ * b));
    }

    {   // now r1, C1 ...
        assert(x_coef(0) > 0); // sanity check
        const int n = abel_n_;

        // increase r1 until C1 is significantly smaller than a_{-1}
        ix_t max_C1 = bmp::lower(ix_t(x_coef(0) / 16));
        for (r1 = r + 1; ; r1 += 1) {
            ix_t iC1 = abs(x_coef(1)) / r1 + C0 * pow(r, -n) / (r1 - r);
            for (int k = 1; k <= n; k++)
                iC1 += k * abs(x_coef(2 + k)) / pow(r1, k + 1);

            if (bmp::upper(iC1) < max_C1) {
                C1 = bmp::upper(iC1);
                break;
            }
        }

        ix_t C1am1 = interval_t(C1 / x_coef(0));
        delta1       = bmp::upper(ix_t(pow(pow(C1am1, -2) - 1, real_t(-1) / 2)));
        varkappa1    = bmp::lower(ix_t(pow(1 - pow(C1am1, 2), real_t(1) / 2)));
    }

    ix_t t_Nstar;
    { // set Nstar and L
      // instead of
      //    Nstar_ = int(ceil(nu + mlogeps_));
      // we slowly increase Nstar until the constraints are satisfied.
        ix_t Ar1 = abel_t(ix_t(r1), x_coef)(0);
        ix_t x = 1, t = 1;

        for (Nstar_ = 0; ; Nstar_++) {
            x = left_inv(x)(0);
            t = pow(x, -gamma_);
            rx_t At = bmp::lower(abel_t(t, x_coef)(0));

            // unremarkable, but we set the number of derivatives L here
            L_ = 1 + Nstar_ / 2;

            // FIXME!! Remove +20, figure out what to do
            rx_t min_At = 20 +  bmp::upper(
                    Ar1 + ix_t(L_) / (exp(ix_t(1)) * pi_x_ * varkappa1)
                    );
            if (    abel_t_in_range(t)
                    && bmp::upper(abel_t_error(t)) < pow(real_t(2), -PREC_)
                    && At > min_At   ) {
                t_Nstar = t;
                break;
            }
        }
    }

    { // Now set the constant term so that A(1) is approximately 0.
      // Not trying to control the accuracy
        ix_t x = 1;
        for (int i=0; i<Nstar_; i++) {
            x = left_inv(x)(0);
        }

        ix_t t = pow(x, -gamma_);

        // now, we should have A(t) = Nstar_
        x_coef(2) = 0;
        x_coef(2) = Nstar_ - bmp::median(abel_t(t, x_coef)(0));
    }

    // populate the coefficients
    for (int k = 0; k <= abel_n_ + 2; k++) {
        abel_coef_(k) = interval_t(x_coef(k));
        abel_coef_ni_(k) = real_t(bmp::median(x_coef(k)));
    }

    // we computed the constants in extended precision,
    // we revert to normal precision
    auto round_up = [] (const rx_t &x) {
        return real_t(bmp::upper(ix_t(x)));
    };
    auto round_down = [] (const rx_t &x) {
        return real_t(bmp::lower(ix_t(x)));
    };
    abel_r_  = round_up(r);
    abel_C0_ = round_up(C0);
    abel_r1_ = round_up(r1);
    abel_C1_ = round_up(C1);
    abel_delta1_ = round_up(delta1);
    abel_varkappa1_ = round_down(varkappa1);

    assert( abel_r_ > 0 );
    assert( abel_C0_ > 0 );
    assert( abel_r1_ > abel_r_ );
    assert( abel_C1_ > 0 );
    assert( abel_delta1_ > 0 );
    assert( abel_varkappa1_ > 0 );

    abel_am1_minus_C1_ = bmp::lower(interval_t(
                x_coef(0) - C1
                ));
    abel_nu_ = bmp::lower(interval_t(
                abel_t(t_Nstar, x_coef)(0) - abel_t(ix_t(r1), x_coef)(0)
                ));

    assert( abel_am1_minus_C1_ > 0 );
    assert( abel_nu_ > 0 );

    return;
}

template <int PREC>
LSV<PREC>::MatrixXix LSV<PREC>::abel_matrix() const {
    // X is a KAbel+1 by KAbel+1 matrix where the columns contain coefficients at 1, 1/t, 1/t^2, ... for:
    // * first column: coefficients of t(left(z))
    // * second column: coefficients of log(t(left(z)))
    // * 2+n column: coefficients of t(left(z))^{-k} - t^{-k}
    MatrixXix X = MatrixXix::Zero(abel_n_ + 2, abel_n_ + 2);
    const interval_extra_t g = gamma_, b = pow(2, g);
    interval_extra_t c;

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
const LSV<PREC>::interval_extra_t LSV<PREC>::pi_x_ = 4 * atan(interval_extra_t(1));

template <int PREC>
const LSV<PREC>::real_t LSV<PREC>::real_eps_ = std::numeric_limits<real_t>::epsilon();

template <int PREC>
const LSV<PREC>::real_t LSV<PREC>::real_eps_sqrt_ = sqrt(std::numeric_limits<real_t>::epsilon());

}  // namespace lsv_ns

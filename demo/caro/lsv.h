#pragma once

#include <limits>
#include <vector>
#include <complex>
#include <cmath> // for pow, abs, log, ... for standard types
#include <type_traits>
#include <Eigen/Core>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpfi.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/eigen.hpp>

// for arbitrary precision integers (for Bernoulli)
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

        typedef Cheb<interval_t> real_cheb_t;
    private:
        static const interval_t pi_;
        static const real_t real_eps_;
        static const real_t real_eps_sqrt_;

        interval_t gamma_;
        VectorXi abel_coef_;
        VectorXr  abel_coef_ni_;

        void compute_abel_stuff();
        MatrixXix abel_matrix() const;

        static constexpr int PREC_ = PREC;
        interval_t kappa_, mlogeps_, tau_, W_;
        int N_, Nstar_, L_, KAbel_, halfM_, M_, Mhat_, P_;

        // Constants from the Abel function, with \tA = \tilde{A},
        // |\tA(z) - \tA_n(z)| \leq abel_C0_ |t|^{-n}  when  Re(t) \geq abel_r_
        real_t abel_r_, abel_C0_;
        // |\tA'(t) - a_{-1}| \leq abel_C1_ when Re(t) \leq abel_r1_
        real_t abel_r1_, abel_C1_;
        // derived variables
        real_t abel_delta1_,        // = abel_C1_ / sqrt(a_{-1}^2 - abel_C1_^2)
               abel_varkappa1_;     // = (1 + abel_delta1_^2)^{-1/2}  or  (1 - C_1^2 / a_{-1}^2)^{1/2}

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

            const int n = KAbel_ - 1;
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

        VectorXci taylorpoints_, taylorpointweights_;
        MatrixXi integralccpoints_, integralccweights_;

        void compute_integral_derivative_weights() {
            // check that required constants are set
            assert(halfM_ > 0 && Nstar_ > 0 && tau_ > 0 && M_ > 0 && L_ > 0 && P_ > 0 && W_ > 0);

            {  // derivatives
                std::cout << "Computing derivatives, halfM_: " << halfM_ << "...\n";
                taylorpoints_.resize(halfM_);
                for (int k = 0; k < halfM_; k++)
                    taylorpoints_(k) = complex_interval_t(Nstar_, 0)
                        + tau_ * exp(complex_interval_t(0, (2 * k + 1) * pi_ / M_));
                MatrixXci FFTweights(halfM_, L_);
#pragma omp parallel for collapse(2) schedule(dynamic)
                for (int m = 0; m < halfM_; m++)
                    for (int k = 0; k < L_; k++)
                        FFTweights(m, k) =
                            exp(complex_interval_t(0, - (pi_ * (2 * k + 1) * (2 * m + 1)) / M_))
                            / complex_interval_t(halfM_);
                std::cout << "   ... computing bernoulli2k(" << L_ << ") ...\n";
                VectorXi EMderivativeweights = -bernoulli2k(L_).tail(L_);
                std::cout << "   ... bernoulli2k done ...\n";
                for (int k = 0; k < L_; k++)
                    EMderivativeweights(k) /= pow(tau_, 2*k + 1) * 2 * (k + 1);
                taylorpointweights_ = FFTweights * EMderivativeweights;
                { // DEBUG
                    std::cout << "error in taylor point weights: " << uncertainty_c(taylorpointweights_) << "\n"
                        << "error in FFTweights: " << uncertainty_c(FFTweights) << "\n"
                        << "error in EMderivativeweights: " << uncertainty_c(EMderivativeweights) << "\n"
                        << "error in bernoulli: " << uncertainty_c(bernoulli2k(L_)) << "\n";
                }
                std::cout << "   ... done\n";
            }

            {  // integral
                std::cout << "Computing integral...\n";
                std::cout << "   ... computing CCQ(" << Mhat_ << ") ...\n";
                ClenshawCurtisQuadrature<interval_t> cc(Mhat_);
                std::cout << "   ... computing weights...\n";
                integralccpoints_.resize(Mhat_, P_);
                integralccweights_.resize(Mhat_, P_);
#pragma omp parallel for collapse(2) schedule(dynamic)
                for (int m = 0; m < Mhat_; m++) {
                    for (int p = 0; p < P_; p++) {
                        integralccpoints_(m, p) = Nstar_ * exp(2 * W_ * (cc.nodes()(m) + p));
                        integralccweights_(m, p) = integralccpoints_(m, p) * 2 * W_ * cc.weights()(m);
                    }
                }
                std::cout << "   ... done\n";
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

            // Caro's BASIC CONSTANTS
            interval_t nu = 10,
                   Rad = interval_t(80) / interval_t(100);

            kappa_ = interval_t(1) / gamma_;
            mlogeps_ = log(interval_t(2)) * PREC_;
            N_ = int(ceil(mlogeps_ / Rad));
            KAbel_ = int(ceil(mlogeps_));

            std::cout << "KAbel_: " << KAbel_ << "\n";
            std::cout << "Computing Abel function coeffs and errors...\n";
            compute_abel_stuff();
            std::cout << "   ... done...\n";
            std::cout << "Nstar_: " << Nstar_ << "\n";
            std::cout << "\nError constants: "
                << "  abel_r_: " << abel_r_ << ",   abel_C0_: " << abel_C0_ << "\n";
            std::cout << "First Abel coefficients: " << abel_coef_.head(5).transpose() << "\n\n";


            // L is the number of derivatives FIXME: what should it be?
            L_ = 1 + Nstar_ / 2;
            // M is the number of points on the circle
            halfM_ = L_;
            M_ = 2 * halfM_;
            // tau is the radius of the circle
            tau_ = (Nstar_ - nu) * exp(-interval_t(1));

            assert(M_ >= 2*L_);
            assert(exp(interval_t(1)) * 2 * pi_ * (Nstar_ - nu) > 2 * L_);

            Mhat_ = 2 * M_;
            W_ = 1;
            P_ = int(ceil(mlogeps_ / W_ / kappa_));

            compute_integral_derivative_weights();

        }

        int NCheb() const { return N_; };
        int prec() const { return PREC_; };
        int KAbel() const { return KAbel_; };
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
            real_cheb_t cheb(0.5, 1.0, N_);

            const VectorXi x_nodes = cheb.nodes();
            assert(x_nodes.size() == N_);

            MatrixXi L_values(N_, N_);
#pragma omp parallel for shared(L_values)
            for (int ix = 0; ix < x_nodes.size(); ix++) {
                interval_t x = x_nodes[ix];
                VectorXi r = VectorXi::Zero(N_);

                // evaluation of all chebs
                auto v = [&cheb, N = N_]<typename var_t>(const var_t &x) {
                    return cheb.basis_values_trig(x, N);
                };

                auto evaluate_branch_real = [&v, this](const interval_t &v0xn, const interval_t &dv0xn,
                                                       const interval_t &weight) -> VectorXi {
                    Vector2i y = right_inv(v0xn);
                    return v(y(0)) * (weight * dv0xn * y(1));
                };
                auto evaluate_branch_cplx = [&v, this](const complex_interval_t &v0xn, const complex_interval_t &dv0xn,
                                                       const complex_interval_t &weight) -> VectorXci {
                    Vector2ci y = right_inv(v0xn);
                    return v(y(0)) * (weight * dv0xn * y(1));
                };

                // add branches naively up to Nstar
                interval_t v0xn = x, dv0xn = 1;

                for (int n = 0; n < Nstar_; n++) {
                    r += evaluate_branch_real(v0xn, dv0xn, 1);
                    Vector2i y = left_inv(v0xn);
                    v0xn = y(0);
                    dv0xn *= y(1);
                }

                // 1/2 sum with the Nstar
                r += evaluate_branch_real(v0xn, dv0xn, 0.5);
                { // DEBUG
                    auto ee = uncertainty(r);
                    if (ee > 1) {
                        std::cout << "Accuracy at " << __LINE__ << ": " << ee << "\n";
                    }
                }

                // figure out Abel function of x and its derivative
                Vector2i Av0xn = abel(v0xn);
                interval_t Ax = Av0xn(0) - Nstar_,
                           dAx = Av0xn(1) * dv0xn;
                if (false) {
                    std::cout << "v0xn: " << v0xn << ", " << dv0xn
                        << ", Av0xn: " << Av0xn.transpose()
                        << ", Ax: " << Ax << ", dAx: " << dAx
                        << "\n        *badness*: " << bmp::width(Ax) + bmp::width(dAx)
                        << "\n";
                }

                // add the derivatives (all in one go)
                for (int m = 0; m < halfM_; m++) {
                    const complex_interval_t &n_cplx = taylorpoints_(m);
                    const complex_interval_t &weight = taylorpointweights_(m);

                    Vector2ci Aixn = abel_inv(complex_interval_t(Ax) + n_cplx);
                    complex_interval_t dv0xn = dAx * Aixn(1);

                    { // DEBUG
                        auto ee = uncertainty_c(evaluate_branch_cplx(Aixn(0), dv0xn, weight));
                        auto width = [] (const complex_interval_t &x) {
                            return bmp::width(x.real()) + bmp::width(x.imag());
                        };
                        if (ee > 1) {
                            Vector2ci ri = right_inv(Aixn(0));
                            std::cout << "Accuracy at " << __LINE__ << ": " << ee << "\n"
                                << "Aixn(0): " << Aixn(0) << ", width: " << width(Aixn(0)) << "\n"
                                << "dv0xn: " << dv0xn << ", width: " << width(dv0xn) << "\n"
                                << "weight: " << weight << ", width: " << width(weight) << "\n"
                                << "ri(0):  " << ri(0)
                                << ", width: " << width(ri(0)) + width(ri(1)) << "\n"
                                << "N_: " << N_ << "\n";
                            for (int j=0; j<N_; j++) {
                                complex_interval_t c1, c2;
                                c1 = cheb.basis_value(ri(0), j);
                                c2 = cheb.basis_value_trig(ri(0), j);

                                interval_t a = abs((c1 - c2).real()) + abs((c1 - c2).imag());

                                real_t rr = bmp::upper(a);

                                if (rr > 1e-8) {
                                    std::cout << "j: " << j << ", c1: " << c1 << ", c2: " << c2
                                        << ", widths: " << width(c1) << "  vs  " << width(c2)
                                        << "\n";
                                } else {
                                    std::cout << "j: " << j << " ok!\n";
                                }

                            }
                            assert(false);
                        }
                    }

                    r += evaluate_branch_cplx(Aixn(0), dv0xn, weight).real();
                }
                { // DEBUG
                    auto ee = uncertainty(r);
                    if (ee > 1) {
                        std::cout << "Accuracy at " << __LINE__ << ": " << ee << "\n";
                    }
                }

                // do the integral, again as a linear combination of point evaluations
                // we use an exp transformation plus Clenshaw-Curtis
                for (int m = 0; m < Mhat_; m++) {
                    for (int p = 0; p < P_; p++) {
                        const interval_t &n = integralccpoints_(m, p);
                        const interval_t &weight = integralccweights_(m, p);

                        Vector2i Aixn = abel_inv<interval_t>(Ax + n);
                        interval_t dv0xn = dAx * Aixn(1);

                        { // DEBUG
                            auto ee = uncertainty(evaluate_branch_real(Aixn(0), dv0xn, weight));
                            if (ee > 1) {
                                std::cout << "Accuracy at " << __LINE__ << ": " << ee << "\n"
                                    << "Aixn(0): " << Aixn(0) << ", width: " << bmp::width(Aixn(0)) << "\n"
                                    << "dv0xn: " << dv0xn << ", width: " << bmp::width(dv0xn) << "\n"
                                    << "weight: " << weight << ", width: " << bmp::width(weight) << "\n";

                                assert(false);
                            }
                        }
                        r += evaluate_branch_real(Aixn(0), dv0xn, weight);
                    }
                }
#pragma omp critical
                L_values.col(ix) = r;
            }

            MatrixXi R(N_, N_);
            for (int k = 0; k < N_; k++) {
                real_cheb_t c(0.5, 1.0, N_);
                c.set_from_values(L_values.row(k).transpose());
                R.col(k) = c.coef();
            }

            return R;
        }

        // Approximate Abel function A(x), satisfying A(left(x)) = A(x) - 1.
        // With t = t(x) = x^{-gamma}, and N = KAbel-1, we use an asymptotic approximation
        //      A(x) ≈ am1 t + al log t + a0 + a1 / t + a2 / t^2 + ... + aN / t^N
        // which is valid for sufficiently large Re(t) with a controlled error
        template <typename var_t> requires type_one_of<var_t, interval_t, complex_interval_t>
        Vector2<var_t> abel(const var_t &x) const {
            var_t t = pow(x, -gamma_);
            Vector2<var_t> A = abel_t(t);
            return Vector2<var_t>(A(0), -gamma_ * (t / x) * A(1));
        }
        template <typename var_t, typename coef_t> requires
            (type_one_of<var_t, real_t, complex_t> && type_one_of<coef_t, VectorXr>) ||
            (type_one_of<var_t, interval_t, complex_interval_t> && type_one_of<coef_t, VectorXi>) ||
            (type_one_of<var_t, interval_extra_t> && type_one_of<coef_t, VectorXix>)
        Vector2<var_t> abel_t(const var_t &t, const coef_t &coef) const {
            assert(coef.size() == KAbel_ + 2);
            var_t A = coef(0) * t + coef(1) * log(t) + coef(2),
                  dA = coef(0) + coef(1) / t;

            const var_t ti = var_t(1) / t;
            var_t tj = ti;
            for (int j = 1; j <= KAbel_ - 1; j++) {
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

        // Akiyama–Tanigawa algorithm for second Bernoulli numbers B+n
        static VectorXi bernoulli2k(int p) {
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

    using ix_t = interval_extra_t;

    abel_coef_.resize(KAbel_ + 2);
    abel_coef_ni_.resize(KAbel_ + 2);

    // We do the computation in extra precision first, then dumb it down
    VectorXix x_coef(KAbel_ + 2);

    { // first compute non-constant coefficients (am1, al, a1, a2, ...) of the Abel function
        VectorXix b = VectorXix::Zero(KAbel_ + 1);
        b(0) = -1;

        MatrixXix A = abel_matrix(); // lower triangular

        // brute force solve lower triangular matrix
        // instead of using Eigen that can fail in interval arithmetic:
        //VectorXi x = abel_matrix().template triangularView<Eigen::Lower>().solve(b);
        //VectorXi x = interval_root_ns::linear_krawczyk(A, b);

        VectorXix xx(KAbel_ + 1);
        for (int k = 0; k < KAbel_ + 1; k++) {
            ix_t s = 0;
            for (int j = 0; j < k; j++) {
                s += A(k,j) * xx(j);
            }
            xx(k) = (b(k) - s) / A(k,k);
        }

        x_coef(0) = xx(0);
        x_coef(1) = xx(1);
        x_coef(2) = 0;
        for (int j = 2; j <= KAbel_; j++) {
            x_coef(j + 1) = xx(j);
        }
    }

    { // compute abel error constants
        int n = KAbel_ - 1;

        ix_t q = ix_t(1) / 2,
             b = pow(2, gamma_);

        ix_t r = max(
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

        // OUTPUT:
        abel_r_ = bmp::upper(interval_t(
                    r
                    ));
        abel_C0_ = bmp::upper(interval_t(
                    pow(interval_t(2), n + 1) * M * I / (gamma_ * b)
                    ));
        std::cout << "abel_r_: " << abel_r_ << "\n"
            << "abel_C0_: " << abel_C0_ << "\n";

        // now abel_r1_, abel_C1_ ...
        assert(x_coef(0) > 0); // sanity check
        ix_t C1, r1;
        r1 = abel_r_;

        for (;;) {
            r1 += 1;
            C1 = abs(x_coef(1)) / r1 + abel_C0_ * pow(abel_r_, -n) / (r1 - abel_r_);
            for (int k = 1; k <= n; k++)
                C1 += k * abs(x_coef(2 + k)) / pow(r1, k + 1);
            if (bmp::upper(ix_t(C1 / x_coef(0))) < ix_t(0.25)) { // arbitrary small constant
                break;
            }
        }
        abel_r1_ = bmp::upper(interval_t(r1));
        abel_C1_ = bmp::upper(interval_t(C1));

        real_t C1am1 = bmp::upper(interval_t(C1 / x_coef(0)));
        abel_delta1_ = pow(pow(C1am1, -2) - 1, real_t(-1) / 2);
        abel_varkappa1_ = pow(1 - pow(C1am1, 2), real_t(1) / 2);

        std::cout << "abel_r1_: " << abel_r1_ << "\n"
            << "abel_C1_: " << abel_C1_ << "\n"
            << "abel_delta1_: " << abel_delta1_ << "\n"
            << "abel_varkappa1_: " << abel_varkappa1_ << "\n";
    }

    { // set Nstar empirically, instead of
      // Nstar_ = int(ceil(nu + mlogeps_));
        ix_t x = 1, t = 1;
        for (Nstar_ = 0; ; Nstar_++) {
            x = left_inv(x)(0);
            t = pow(x, -gamma_);
            if (abel_t_in_range(t) && bmp::upper(abel_t_error(t)) < pow(real_t(2), -PREC_)) {
                std::cout << "Nstar_ = " << Nstar_ << ", t: " << t << "\n";
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
    for (int k = 0; k <= KAbel_ + 1; k++) {
        abel_coef_(k) = interval_t(x_coef(k));
        abel_coef_ni_(k) = real_t(bmp::median(x_coef(k)));
    }
}

template <int PREC>
LSV<PREC>::MatrixXix LSV<PREC>::abel_matrix() const {
    // X is a KAbel+1 by KAbel+1 matrix where the columns contain coefficients at 1, 1/t, 1/t^2, ... for:
    // * first column: coefficients of t(left(z))
    // * second column: coefficients of log(t(left(z)))
    // * 2+n column: coefficients of t(left(z))^{-k} - t^{-k}
    MatrixXix X = MatrixXix::Zero(KAbel_ + 1, KAbel_ + 1);
    const interval_extra_t g = gamma_, b = pow(2, g);
    interval_extra_t c;

    c = - g * b;
    for (int j = 0; j <= KAbel_; j++) {
        X(j,0) = c;
        c *= b * (-g - j - 1) / (j + 2);
    }
    c = - g * b;
    for (int j = 1; j <= KAbel_; j++) {
        X(j,1) = c / j;
        c *= -b;
    }
    for (int k = 1; k <= KAbel_ - 1; k++) {
        c = 1;
        for (int j = k + 1; j <= KAbel_; j++) {
            c *= (k * g - j + k + 1) * b / (j - k);
            X(j, k+1) = c;
        }
    }

    return X;
}

template <int PREC>
const LSV<PREC>::interval_t LSV<PREC>::pi_ = 4 * atan(interval_t(1));

template <int PREC>
const LSV<PREC>::real_t LSV<PREC>::real_eps_ = std::numeric_limits<real_t>::epsilon();

template <int PREC>
const LSV<PREC>::real_t LSV<PREC>::real_eps_sqrt_ = sqrt(std::numeric_limits<real_t>::epsilon());

}  // namespace lsv_ns

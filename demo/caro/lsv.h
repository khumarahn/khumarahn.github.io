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

#include "clenshawcurtis.h"
#include "cheb.h"
#include "interval-root.h"

// Range checking is disabled if NDEBUG or EIGEN_NO_DEBUG is defined

namespace lsv_ns {

namespace bmp = boost::multiprecision;

template <typename T, typename var_t>
concept real_or_complex = std::is_same_v<T, var_t> || std::is_same_v<T, std::complex<var_t>>;

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

        interval_t alpha_, gamma_;
        VectorXi abel_coef_;
        VectorXr  abel_coef_ni_;

        void compute_abel_coef();
    public: // TEMPORARY
        MatrixXix abel_matrix() const;
    private: // REMOVE

        // make public for debug
        static constexpr int PREC_ = PREC;
        interval_t kappa_, mlogeps_, tau_, W_;
        int N_, Nstar_, K_, KAbel_, halfM_, M_, Mhat_, P_;

        VectorXci taylorpoints_, taylorpointweights_;
        MatrixXi integralccpoints_, integralccweights_;

    public: // Newton methods aare TEMPORARY PUBLIC
        // interval version of Newton method. f is assumed convex, and the root is in [a,b]
        //Vector2i interval_newton(const std::function<Vector2i (interval_t)>& f,
        //        const interval_t &guess) const;
        // ... and a Krawczyk method for complex f
        Vector2ci complex_krawczyk(const std::function<Vector2ci (complex_interval_t)>& f,
                const complex_interval_t &guess) const;
        // .. and a Krawczyk method for a linear equation Ax=b with invertible A
        VectorXi linear_krawczyk(const MatrixXi &A, const VectorXi &b) const;
    public:
        LSV(const interval_t &gamma = 1) {
            set_gamma(gamma);
        }
        void set_gamma(const interval_t &gamma) {
            gamma_ = gamma;
            alpha_ = interval_t(1) / gamma_;

            // Caro's BASIC CONSTANTS
            interval_t nu = 10,
                   Rad = interval_t(80) / interval_t(100);

            kappa_ = alpha_;
            mlogeps_ = log(interval_t(2)) * PREC_;
            N_ = int(ceil(mlogeps_ / Rad));
            Nstar_ = int(ceil(nu + mlogeps_));
            K_ = int(floor(((Nstar_ - nu) - 1) * 0.5));
            KAbel_ = int(ceil(mlogeps_));
            halfM_ = K_;
            M_ = 2 * halfM_;
            tau_ = (Nstar_ - nu) * exp(-interval_t(1));
            assert(M_ >= 2*K_);
            assert(exp(interval_t(1)) * 2 * pi_ * (Nstar_ - nu) > 2 * K_);

            Mhat_ = 2 * M_;
            W_ = 1;
            P_ = int(ceil(mlogeps_ / W_ / kappa_));

            {  // derivatives
                taylorpoints_.resize(halfM_);
                for (int k = 0; k < halfM_; k++)
                    taylorpoints_(k) = complex_interval_t(Nstar_, 0) + tau_ * exp(complex_interval_t(0, (2 * k + 1) * pi_ / M_));
                MatrixXci FFTweights(halfM_, K_);
                for (int m = 0; m < halfM_; m++)
                    for (int k = 0; k < K_; k++)
                        FFTweights(m, k) =
                            exp(complex_interval_t(0, -pi_ * (2 * k + 1) * (2 * m + 1) / interval_t(M_))) / complex_interval_t(halfM_);
                VectorXi EMderivativeweights = -bernoulli2k(K_).tail(K_);
                for (int k = 0; k < K_; k++)
                    EMderivativeweights(k) /= pow(tau_, 2*k + 1) * 2 * (k + 1);
                taylorpointweights_ = FFTweights * EMderivativeweights;
            }

            {  // integral
               ClenshawCurtisQuadrature<interval_t> cc(Mhat_);
               integralccpoints_.resize(Mhat_, P_);
               integralccweights_.resize(Mhat_, P_);
               for (int m = 0; m < Mhat_; m++) {
                   for (int p = 0; p < P_; p++) {
                       integralccpoints_(m, p) = Nstar_ * exp(2 * W_ * (cc.nodes()(m) + p));
                       integralccweights_(m, p) = integralccpoints_(m, p) * 2 * W_ * cc.weights()(m);
                   }
               }
            }

            compute_abel_coef();
        }

        int NCheb() const { return N_; };
        int prec() const { return PREC_; };
        int KAbel() const { return KAbel_; };
        interval_t alpha() const { return alpha_; };
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
            interval_t md = left(0.5)(1);
            // an interval enclosing the root
            interval_t guess(x / md, x);

            auto f = [this, &x] (const interval_t &z) {
                Vector2i r = left(z);
                r(0) -= x;
                return r;
            };
            return interval_root_ns::interval_newton(f, guess);
        }

        template <real_or_complex<interval_t> var_t>
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
//#pragma omp parallel for shared(L_values)
            for (int ix = 0; ix < x_nodes.size(); ix++) {
                interval_t x = x_nodes[ix];
                VectorXi r = VectorXi::Zero(N_);

                // evaluation of all chebs
                auto v = [&cheb, N = N_]<typename var_t>(const var_t &x) {
                    return cheb.basis_values(x, N);
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
                //std::cout << "v0xn: " << v0xn << ", width: " << bmp::width(v0xn) << "\n";

                // 1/2 sum with the Nstar
                r += evaluate_branch_real(v0xn, dv0xn, 0.5);

                // figure out Abel function of x and its derivative
                Vector2i Av0xn = abel(v0xn);
                //std::cout << "Av0xn: " << Av0xn.transpose() << "\n";
                //std::cout << "Nstar_: " << Nstar_ << "\n";
                interval_t Ax = Av0xn(0) - Nstar_,
                           dAx = Av0xn(1) * dv0xn;
                //std::cout << "Ax: " << Ax << "\n";

                // add the derivatives (all in one go)
                for (int m = 0; m < halfM_; m++) {
                    const complex_interval_t &n_cplx = taylorpoints_(m);
                    const complex_interval_t &weight = taylorpointweights_(m);

                    Vector2ci Aixn = abel_inv(complex_interval_t(Ax) + n_cplx);
                    complex_interval_t dv0xn = dAx * Aixn(1);

                    r += evaluate_branch_cplx(Aixn(0), dv0xn, weight).real();
                }

                // do the integral, again as a linear combination of point evaluations
                // we use an exp transformation plus Clenshaw-Curtis
                for (int m = 0; m < Mhat_; m++) {
                    for (int p = 0; p < P_; p++) {
                        const interval_t &n = integralccpoints_(m, p);
                        const interval_t &weight = integralccweights_(m, p);

                        Vector2i Aixn = abel_inv<interval_t>(Ax + n);
                        interval_t dv0xn = dAx * Aixn(1);

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

        // Abel function
        // with t = t(x) = x^{-gamma}
        // A(x) = am1 t + al log t + a0 + a1 / t + a2 / t^2 + ...
        // satisfies A(left(x)) = A(x) - 1
        // We compute it with KAbel + 2 coefficients am1, al, a0, ..., a{KAbel-1}
        template <real_or_complex<interval_t> var_t>
        Vector2<var_t> abel(const var_t &x) const {
            var_t t = pow(x, -gamma_);
            Vector2<var_t> A = abel_t(t);
            return Vector2<var_t>(A(0), -gamma_ * (t / x) * A(1));
        }
        template <real_or_complex<interval_t> var_t>
        Vector2<var_t> abel_t(const var_t &t) const {
            //std::cout << "called interval abel_t(" << t << ")\n";
            var_t A = abel_coef_(0) * t + abel_coef_(1) * log(t) + abel_coef_(2),
                  dA = abel_coef_(0) + abel_coef_(1) / t;

            const var_t ti = var_t(1) / t;
            var_t tj = ti;
            for (int j = 1; j <= KAbel_ - 1; j++) {
                A += abel_coef_(2 + j) * tj;
                if (0) if constexpr (std::is_same_v<var_t, complex_interval_t>) {
                    std::cout << "    in interval abel_t: increment: " << abel_coef_(2 + j) * tj
                        << " of width: " << bmp::width((abel_coef_(2 + j) * tj).real())
                        << " x " << bmp::width((abel_coef_(2 + j) * tj).imag())
                        << "\n";
                } else if constexpr (std::is_same_v<var_t, interval_t>) {
                    std::cout << "    in interval abel_t: increment: " << abel_coef_(2 + j) * tj
                        << " of width: " << bmp::width(interval_t(abel_coef_(2 + j) * tj))
                        << "\n";
                }
                tj *= ti;
                dA -= var_t(j) * abel_coef_(2 + j) * tj;
            }
            return Vector2<var_t>(A, dA);
        }
        template <real_or_complex<real_t> var_t>
        Vector2<var_t> abel_t(const var_t &t) const {
            //std::cout << "called abel_t(" << t << ")\n";
            var_t A = abel_coef_ni_(0) * t + abel_coef_ni_(1) * log(t) + abel_coef_ni_(2),
                  dA = abel_coef_ni_(0) + abel_coef_ni_(1) / t;

            const var_t ti = var_t(1) / t;
            var_t tj = ti;
            for (int j = 1; j <= KAbel_ - 1; j++) {
                A += abel_coef_ni_(2 + j) * tj;
                tj *= ti;
                dA -= var_t(j) * abel_coef_ni_(2 + j) * tj;
            }
            return Vector2<var_t>(A, dA);
        }

        // inverse
        // * non-interval version, Newton method
        template <real_or_complex<real_t> var_t>
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
            //std::cout << "called abel_t_inv(" << a << ")\n";
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
        template <real_or_complex<interval_t> var_t>
        Vector2<var_t> abel_inv(const var_t &a) const {
            //std::cout << "called abel_inv(" << a << ")\n";
            Vector2<var_t> Ai = abel_t_inv(a);
            const var_t & t = Ai(0), dAt = Ai(1);
            var_t x = pow(t, var_t(-1 / gamma_));
            return Vector2<var_t>(x, -(x / t) / gamma_ * dAt);
        }

        // Akiyama–Tanigawa algorithm for second Bernoulli numbers B+n
        static VectorXi bernoulli2k(int p) {
            VectorXi B2(p+1),
                     A(2*p+1);
            interval_t fact = 1;
            for (int m = 0; m <= 2*p; m++) {
                A(m) = fact;
                fact *= m + 1;
                for (int j = m; j >= 1; j--)
                    A(j - 1) = j * ((m + 1) * A(j - 1) - A(j));
                if (m % 2 == 0)
                    B2[m / 2] = A(0) / fact;
            }
            return B2;
        }

};  // class LSV

template <int PREC>
void LSV<PREC>::compute_abel_coef() {
    abel_coef_.resize(KAbel_ + 2);

    // first compute non-constant coefficients (am1, al, a1, a2, ...) of the Abel function
    {
        VectorXix b = VectorXix::Zero(KAbel_ + 1);
        b(0) = -1;

        MatrixXix A = abel_matrix(); // lower triangular

        // brute force solve lower triangular matrix
        // instead of using Eigen that can fail in interval arithmetic:
        //VectorXi x = abel_matrix().template triangularView<Eigen::Lower>().solve(b);
        
        VectorXix xx(KAbel_ + 1);
        for (int k = 0; k < KAbel_ + 1; k++) {
            interval_extra_t s = 0;
            for (int j = 0; j < k; j++) {
                s += A(k,j) * xx(j);
            }
            xx(k) = (b(k) - s) / A(k,k);
        }
        
        //VectorXi x = interval_root_ns::linear_krawczyk(A, b);

        abel_coef_(0) = interval_t(xx(0));
        abel_coef_(1) = interval_t(xx(1));
        for (int j = 2; j<=KAbel_; j++) {
            abel_coef_(j+1) = interval_t(xx(j));
        }
    }

    // Now set the constant term so that A(1) is approximately 0.
    // This is rather dubious: say, 1e+07 is approximately 0
    {
        interval_t x = 1.0;
        int NAbel = 5 * KAbel_;
        for (int i = 1; i <= NAbel; i++) {
            x = left_inv(x)(0);
        }
        //
        interval_t t = pow(x, -gamma_);

        // now, A(t) = NAbel
        abel_coef_(2) = 0;
        abel_coef_(2) = NAbel - bmp::median(abel_t(t)(0));
    }

    // populate the non-interval coefficients
    abel_coef_ni_.resize(KAbel_ + 2);
    for (int k = 0; k < KAbel_ + 1; k++) {
        abel_coef_ni_(k) = bmp::median(abel_coef_(k));
    }
}

template <int PREC>
LSV<PREC>::MatrixXix LSV<PREC>::abel_matrix() const {
    // X is a KAbel+1 by KAbel+1 matrix where the columns contain coefficients at 1, 1/t, 1/t^2, ... for:
    // * first column: coefficients of t(left(z))
    // * second column: coefficients of log(t(left(z)))
    // * 2+n column: coefficients of t(left(z))^{-k} - t^{-k}
    MatrixXix X = MatrixXix::Zero(KAbel_ + 1, KAbel_ + 1);
    interval_extra_t twopowgamma_ = pow(2, interval_extra_t(gamma_)),
                     c = - interval_extra_t(gamma_) * twopowgamma_;
    for (int j = 0; j <= KAbel_; j++) {
        X(j,0) = c;
        c *= twopowgamma_ * (-gamma_ - j - 1) / (j + 2);
    }
    c = - gamma_ * twopowgamma_;
    for (int j = 1; j <= KAbel_; j++) {
        X(j,1) = c / j;
        c *= -twopowgamma_;
    }
    for (int k = 1; k <= KAbel_ - 1; k++) {
        c = 1;
        for (int j = k + 1; j <= KAbel_; j++) {
            c *= (k * gamma_ - j + k + 1) * twopowgamma_ / (j - k);
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

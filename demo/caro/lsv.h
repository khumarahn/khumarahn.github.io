#pragma once

#include <limits>
#include <vector>
#include <complex>
#include <cmath> // for pow, abs, log, ... for standard types
#include <type_traits>
#include <Eigen/Core>

#include "clenshawcurtis.h"
#include "cheb.h"

// Range checking is disabled if NDEBUG or EIGEN_NO_DEBUG is defined

namespace lsv_ns {

template <typename T, typename var_t>
concept real_or_complex = std::is_same_v<T, var_t> || std::is_same_v<T, std::complex<var_t>>;

using std::pow;
using std::log;
using std::ceil;
using std::floor;
using std::exp;

using cheb_ns::Cheb;

template <typename real_t, int PREC = 0>
class LSV {
    public:
        typedef std::complex<real_t> complex_t;

        template <typename var_t>
            using Vector2 = Eigen::Matrix<var_t,2,1>;
        template <typename var_t>
            using VectorX = Eigen::Matrix<var_t,Eigen::Dynamic,1>;
        template <typename var_t>
            using MatrixX = Eigen::Matrix<var_t,Eigen::Dynamic,Eigen::Dynamic>;
        template <typename var_t>
            using func_t = std::function<var_t (var_t)>;

        typedef Vector2<real_t> Vector2r;
        typedef VectorX<real_t> VectorXr;
        typedef MatrixX<real_t> MatrixXr;
        typedef func_t<real_t> real_func_t;

        typedef Vector2<complex_t> Vector2c;
        typedef VectorX<complex_t> VectorXc;
        typedef MatrixX<complex_t> MatrixXc;
        typedef func_t<complex_t> complex_func_t;

        typedef Cheb<real_t> real_cheb_t;
    private:
        static const real_t pi_;
        static const real_t real_eps_;

        real_t alpha_, gamma_, twopowgamma_;
        VectorXr abel_coef_;

        void compute_abel_coef();
        MatrixXr abel_matrix() const;

        // make public for debug
        static constexpr int PREC_ = (PREC > 0) ? PREC : std::numeric_limits<real_t>::digits; // precision in bits
        real_t kappa_, mlogeps_, tau_, W_;
        int N_, Nstar_, K_, KAbel_, halfM_, M_, Mhat_, P_;

        VectorXc taylorpoints_, taylorpointweights_;
        MatrixXr integralccpoints_, integralccweights_;
    public:
        LSV(const real_t &gamma = 1) {
            set_gamma(gamma);
        }
        void set_gamma(const real_t &gamma) {
            gamma_ = gamma;
            alpha_ = real_t(1) / gamma_;
            twopowgamma_ = pow(2, gamma_);

            // Caro's BASIC CONSTANTS
            real_t nu = 10,
                   Rad = real_t(80) / real_t(100);

            kappa_ = 1 / gamma_;
            mlogeps_ = log(real_t(2)) * PREC_;
            N_ = int(ceil(mlogeps_ / Rad));
            Nstar_ = int(ceil(nu + mlogeps_));
            K_ = int(floor(((Nstar_ - nu) - 1) * 0.5));
            KAbel_ = int(ceil(mlogeps_));
            halfM_ = K_;
            M_ = 2 * halfM_;
            tau_ = (Nstar_ - nu) * exp(-real_t(1));
            assert(M_ >= 2*K_);
            assert(exp(real_t(1)) * 2 * pi_ * (Nstar_ - nu) > 2 * K_);

            Mhat_ = 2 * M_;
            W_ = 1;
            P_ = int(ceil(mlogeps_ / W_ / kappa_));

            {  // derivatives
                taylorpoints_.resize(halfM_);
                for (int k = 0; k < halfM_; k++)
                    taylorpoints_(k) = complex_t(Nstar_, 0) + tau_ * exp(complex_t(0, (2 * k + 1) * pi_ / M_));
                MatrixXc FFTweights(halfM_, K_);
                for (int m = 0; m < halfM_; m++)
                    for (int k = 0; k < K_; k++)
                        FFTweights(m, k) =
                            exp(complex_t(0, -pi_ * (2 * k + 1) * (2 * m + 1) / real_t(M_))) / complex_t(halfM_);
                VectorXr EMderivativeweights = -bernoulli2k(K_).tail(K_);
                for (int k = 0; k < K_; k++)
                    EMderivativeweights(k) /= pow(tau_, 2*k + 1) * 2 * (k + 1);
                taylorpointweights_ = FFTweights * EMderivativeweights;
            }

            {  // integral
               ClenshawCurtisQuadrature<real_t> cc(Mhat_);
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
        real_t alpha() const { return alpha_; };
        real_t gamma() const { return gamma_; };
        VectorXr abel_coef() const { return abel_coef_; };

        // branches
        Vector2r left(const real_t &x) const {
            real_t t = pow(2 * x, gamma_);
            return Vector2r(x * (1 + t), 1 + (gamma_ + 1) * t);
        }
        Vector2r right(const real_t &x) const {
            return Vector2r(2 * x - 1, 2);
        }

        // full map
        Vector2r map(const real_t &x) const {
            return (x < 0.5) ? left(x) : right(x);
        }

        // inverses
        Vector2r left_inv(const real_t &x) const {
            // convex, so Newton should work very well
            real_t y = x,
                   p;
            for (;;) {
                Vector2r l = left(y);
                real_t f = l(0) - x;
                p = l(1);
                real_t z = y - f / p;
                if (z >= y)
                    break;
                y = z;
            }
            return Vector2r(y, 1 / p);
        }
        template <real_or_complex<real_t> var_t>
        Vector2<var_t> right_inv(const var_t &x) const {
            return Vector2<var_t>((x + var_t(1)) / var_t(2), 0.5);
        }

        // transfer operator wrt dx
        real_t L(const real_func_t &v, const real_t &x) const {
            Vector2r y1 = left_inv(x),
                     y2 = right_inv(x);
            return v(y1(0)) * y1(1) + v(y2(0)) * y2(1);
        }
        // transfer operator wrt x^{-gamma} dx
        real_t Lxgamma (const real_func_t &v, const real_t &x) const {
            Vector2r y1 = left_inv(x),
                     y2 = right_inv(x);
            real_t w1 = y1(1) * (y1(0) <= 0 ? real_t(1) : pow(x / y1(0), gamma_)),
                   w2 = y2(1) * (y2(0) <= 0 ? real_t(1) : pow(x / y2(0), gamma_));
            return v(y1(0)) * w1 + v(y2(0)) * w2;
        }

        // transfer operator of induced map with n branches
        real_t Lind(const real_func_t &v, const real_t &x, int n) const {
            real_t r = 0, z = x, w = 1;
            for (int k = 0; k < n; k++) {
                Vector2r yr = right_inv(z);
                r += v(yr(0)) * w * yr(1);
                Vector2r l = left_inv(z);
                z = l(0);
                w *= l(1);
            }
            return r;
        }

        // transfer operator of induced map, the clever one
        // This returns an approximation of the induced transfer operator by an
        // N_ by N_ matrix acting on the Chebyshev coefficients
        MatrixXr Lind() const {
            std::vector<real_cheb_t> chebs(N_);
            for (int k = 0; k < N_; k++) {
                VectorXr c = VectorXr::Zero(N_);
                c(k) = 1;
                chebs[k] = real_cheb_t(c, 0.5, 1.0);
            }

            const VectorXr x_nodes = chebs[0].nodes();
            assert(x_nodes.size() == N_);

            MatrixXr L_values(N_, N_);
#pragma omp parallel for shared(L_values)
            for (int ix = 0; ix < x_nodes.size(); ix++) {
                real_t x = x_nodes[ix];
                VectorXr r = VectorXr::Zero(N_);

                // evaluation of all chebs
                auto v = [&chebs, this]<typename var_t>(var_t x) {
                    VectorX<var_t> ret(N_);
                    for (int k = 0; k < N_; k++)
                        ret(k) = chebs[k](x);
                    return ret;
                };

                auto evaluate_branch_real = [&v, this](const real_t &v0xn, const real_t &dv0xn,
                                                       const real_t &weight) -> VectorXr {
                    Vector2r y = right_inv(v0xn);
                    return v(y(0)) * (weight * dv0xn * y(1));
                };
                auto evaluate_branch_cplx = [&v, this](const complex_t &v0xn, const complex_t &dv0xn,
                                                       const complex_t &weight) -> VectorXc {
                    Vector2c y = right_inv(v0xn);
                    return v(y(0)) * (weight * dv0xn * y(1));
                };

                // add branches naively up to Nstar
                real_t v0xn = x, dv0xn = 1;

                for (int n = 0; n < Nstar_; n++) {
                    r += evaluate_branch_real(v0xn, dv0xn, 1);
                    Vector2r y = left_inv(v0xn);
                    v0xn = y(0);
                    dv0xn *= y(1);
                }

                // 1/2 sum with the Nstar
                r += evaluate_branch_real(v0xn, dv0xn, 0.5);

                // figure out Abel function of x and its derivative
                Vector2r Av0xn = abel(v0xn);
                real_t Ax = Av0xn(0) - Nstar_, dAx = Av0xn(1) * dv0xn;

                // add the derivatives (all in one go)
                for (int m = 0; m < halfM_; m++) {
                    const complex_t &n_cplx = taylorpoints_(m);
                    const complex_t &weight = taylorpointweights_(m);

                    Vector2c Aixn = abel_inv(complex_t(Ax) + n_cplx);
                    complex_t dv0xn = dAx * Aixn(1);

                    r += evaluate_branch_cplx(Aixn(0), dv0xn, weight).real();
                }

                // do the integral, again as a linear combination of point evaluations
                // we use an exp transformation plus Clenshaw-Curtis
                for (int m = 0; m < Mhat_; m++) {
                    for (int p = 0; p < P_; p++) {
                        const real_t &n = integralccpoints_(m, p);
                        const real_t &weight = integralccweights_(m, p);

                        Vector2r Aixn = abel_inv<real_t>(Ax + n);
                        real_t dv0xn = dAx * Aixn(1);

                        r += evaluate_branch_real(Aixn(0), dv0xn, weight);
                    }
                }
#pragma omp critical
                L_values.col(ix) = r;
            }

            MatrixXr R(N_, N_);
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
        template <real_or_complex<real_t> var_t>
        Vector2<var_t> abel(const var_t &x) const {
            var_t t = pow(x, -gamma_);
            Vector2<var_t> A = abel_t(t);
            return Vector2<var_t>(A(0), -gamma_ * (t / x) * A(1));
        }
        template <real_or_complex<real_t> var_t>
        Vector2<var_t> abel_t(const var_t &t) const {
            var_t A = abel_coef_(0) * t + abel_coef_(1) * log(t) + abel_coef_(2),
                  dA = abel_coef_(0) + abel_coef_(1) / t;

            const var_t ti = var_t(1) / t;
            var_t tj = ti;
            for (int j = 1; j <= KAbel_ - 1; j++) {
                A += abel_coef_(2 + j) * tj;
                tj *= ti;
                dA -= var_t(j) * abel_coef_(2 + j) * tj;
            }
            return Vector2<var_t>(A, dA);
        }

        // inverse
        template <real_or_complex<real_t> var_t>
        Vector2<var_t> abel_t_inv(const var_t &a) const {
            var_t t = a / abel_coef_(0);
            Vector2<var_t> A;
            for (int i=0; i < 1000; i++) {
                A = abel_t(t);
                if (abs(A(0) - a) < 100 * real_eps_) {
                    assert(i < 100);
                    break;
                }
                t -= (A(0) - a) / A(1);
            }
            return Vector2<var_t>(t, var_t(1) / A(1));
        }
        template <real_or_complex<real_t> var_t>
        Vector2<var_t> abel_inv(const var_t &a) const {
            Vector2<var_t> Ai = abel_t_inv(a);
            const var_t & t = Ai(0), dAt = Ai(1);
            var_t x = pow(t, var_t(-1 / gamma_));
            return Vector2<var_t>(x, -(x / t) / gamma_ * dAt);
        }

        // Akiyama–Tanigawa algorithm for second Bernoulli numbers B+n
        static VectorXr bernoulli2k(int p) {
            VectorXr B2(p+1),
                     A(2*p+1);
            real_t fact = 1;
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

template <typename real_t, int PREC>
void LSV<real_t, PREC>::compute_abel_coef() {
    abel_coef_.resize(KAbel_ + 2);

    // first compute non-constant coefficients (am1, al, a1, a2, ...) of the Abel function
    {
        VectorXr b = VectorXr::Zero(KAbel_ + 1);
        b(0) = -1;

        VectorXr x = abel_matrix().template triangularView<Eigen::Lower>().solve(b);
        abel_coef_(0) = x(0);
        abel_coef_(1) = x(1);
        for (int j = 2; j<=KAbel_; j++) {
            abel_coef_(j+1) = x(j);
        }
    }

    // now the constant term so that A(1) = 0
    {
        real_t x = 1.0;
        int NAbel = 5 * KAbel_;
        for (int i = 1; i <= NAbel; i++) {
            x = left_inv(x)(0);
        }
        //
        real_t t = pow(x, -gamma_);

        // now, A(t) = NAbel
        abel_coef_(2) = 0;
        abel_coef_(2) = NAbel - abel_t(t)(0);
    }
}

template <typename real_t, int PREC>
LSV<real_t, PREC>::MatrixXr LSV<real_t, PREC>::abel_matrix() const {
    // X is a KAbel+1 by KAbel+1 matrix where the columns contain coefficients at 1, 1/t, 1/t^2, ... for:
    // * first column: coefficients of t(left(z))
    // * second column: coefficients of log(t(left(z)))
    // * 2+n column: coefficients of t(left(z))^{-k} - t^{-k}
    MatrixXr X = MatrixXr::Zero(KAbel_ + 1, KAbel_ + 1);
    real_t c = - gamma_ * twopowgamma_;
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

template <typename real_t, int PREC>
const real_t LSV<real_t, PREC>::pi_ = 4 * atan(real_t(1));

template <typename real_t, int PREC>
const real_t LSV<real_t, PREC>::real_eps_ = std::numeric_limits<real_t>::epsilon();

}  // namespace lsv_ns

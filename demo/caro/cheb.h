#pragma once

#include <cmath> // for sin and atan for standard types
#include <type_traits>
#include <unsupported/Eigen/FFT>
#include <type_traits>

namespace cheb_ns {

using std::atan;
using std::cos;

template <typename T, typename var_t>
concept real_or_complex = std::is_same_v<T, var_t> || std::is_same_v<T, std::complex<var_t>>;

// home made class for basic Chebyshev approximation
template <typename real_t>
class Cheb {
    public:
        typedef Eigen::Matrix<real_t,Eigen::Dynamic,1> VectorXr;
    private:
        static const real_t pi_;

        int N_;
        real_t a_, b_;
        VectorXr coef_;

        // precomputed things
        real_t bpa_, bpa2_, bma_, bmai_, bma2_, bma2i_, Ni_;
        void set_abN(const real_t &a, const real_t &b, int N) {
            a_ = a;
            b_ = b;
            N_ = N;
            assert ( N_ > 0 && a_ < b_);

            bpa_ = b_ + a_;
            bpa2_ = bpa_ / 2;
            bma_ = b_ - a_;
            bmai_ = 1 / bma_;
            bma2_ = bma_ / 2;
            bma2i_ = 1 / bma2_;
            Ni_ = real_t(1) / real_t(N_);
        };
    public:
        Cheb() {};
        Cheb(const std::function<real_t (real_t)> &f, const real_t &a, const real_t &b, int N) {
            set_abN(a, b, N);

            // compute a vector of values of f
            Cheb::VectorXr f_val(N_);
            for (int k = 0; k < N_; k++)
                f_val(k) = f(cos(pi_ * (k + 0.5) * Ni_) * bma2_ + bpa2_);

            set_from_values(f_val);
        };
        Cheb(const VectorXr &coef, const real_t &a, const real_t &b) {
            set_abN(a, b, int(coef.size()));
            coef_ = coef;
        };
        Cheb(const real_t &a, const real_t &b, int N) {
            set_abN(a, b, N);
            coef_ = VectorXr::Zero(N);
        };

        VectorXr coef() const { return coef_; };
        real_t a() const {return a_; };
        real_t b() const {return b_; };
        int N() const { return N_; };

        VectorXr nodes() const {
            VectorXr no;
            no.resize(N_);
            for (int k = 0; k < N_; k++)
                no(k) = cos(pi_ * (k + 0.5) * Ni_) * bma2_ + bpa2_;
            return no;
        }

        void set_from_values(const VectorXr &values) {
            assert(values.size() == N_);

            coef_ = kissDCT2(values) * 2 * Ni_;
        }

        template <real_or_complex<real_t> var_t>
        var_t value(const var_t &x) const {
            // Clenshaw recurrence
            var_t d(0),
                  dd(0),
                  y = (var_t(2) * x - var_t(bpa_)) * var_t(bmai_),
                  y2 = var_t(2) * y;
            for (int j = N_ - 1; j > 0; j--) {
                var_t sv = d;
                d = y2 * d - dd + coef_(j);
                dd = sv;
            }
            return y * d - dd + coef_(0) / var_t(2);
        }
        template <real_or_complex<real_t> var_t>
        var_t operator()(const var_t &x) const {
            return value(x);
        }

        Cheb derivative() const {
            VectorXr d_coef(N_);
            if (N_ > 0)
                d_coef(N_ - 1) = 0;
            if (N_ > 1)
                d_coef(N_ - 2) = 2 * (N_ - 1) * coef_(N_ - 1);
            for (int j = N_ - 2; j > 0; j--)
                d_coef(j - 1) = d_coef(j + 1) + 2 * j * coef_(j);

            d_coef *= bma2i_;

            Cheb r(d_coef, a_, b_);
            return r;
        }

        Cheb integral() const {
            real_t sum = 0.0, fac = 1.0, con = (b_ - a_) / 4;
            VectorXr cint(N_);
            for (int j = 1; j < N_ - 1; j++) {
                cint[j] = con * (coef_[j - 1] - coef_[j + 1]) / j;
                sum += fac * cint[j];
                fac = -fac;
            }
            cint[N_ - 1] = con * coef_[N_ - 2] / (N_ - 1);
            sum += fac * cint[N_ - 1];
            cint[0] = 2.0 * sum;
            return Cheb(cint, a_, b_);
        }

        // discrete cosine transform type 2: eigen's kissfft
        // and a reference naive implementation
        static VectorXr kissDCT2(const VectorXr &v) {
            long N = v.size();
            VectorXr vd = VectorXr::Zero(4 * N);
            for (long i = 0; i < N; i++) {
                vd(2 * i + 1) = v(i);
                vd(4 * N - 2 * i - 1) = v(i);
            }

            // dynamic complex vector
            Eigen::Matrix<std::complex<real_t>,Eigen::Dynamic,1> rc;

            Eigen::FFT<real_t> fft;
            fft.fwd(rc, vd);

            VectorXr r(N);
            for (long i = 0; i < N; i++) {
                r(i) = std::real(rc(i)) / 2;
            }

            return r;
        }
        static VectorXr naiveDCT2(const VectorXr &v) {
            int N = int(v.size());
            real_t Ni = real_t(1) / real_t(N);
            VectorXr r(N);
            for (int k = 0; k < N; k++) {
                r(k) = 0;
                for (int n = 0; n < N; n++) {
                    r(k) += v(n) * cos(pi_ * (n + 0.5) * k * Ni);
                }
            }
            return r;
        }
};

template <typename real_t>
const real_t Cheb<real_t>::pi_ = 4 * atan(real_t(1));

} // namespace cheb_ns

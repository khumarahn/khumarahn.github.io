#pragma once

#include <Eigen/Core>
#include <cmath> // for cos for standard types

using std::cos;

// loosely translated from Julia's QuadratureRules
template <typename real_t>
class ClenshawCurtisQuadrature {
    public:
        typedef Eigen::Matrix<real_t,Eigen::Dynamic,1> VectorXr;
    private:
        static const real_t pi_;
        VectorXr nodes_, weights_;
    public:
        ClenshawCurtisQuadrature(int s) {
            auto b = [](int j, int n) { return (n == 2 * j) ? 1 : 2; };
            auto c = [](int k, int n) -> real_t { return (k % n == 0) ? 1 : 2; };
            auto theta = [](int k, int n) { return (2 * pi_ * k ) / n; };

            auto cctrm = [&b, &theta](int k, int j, int n) {
                long J = j;
                return b(j, n) * cos(j * theta(k,n)) / real_t(4 * J * J - 1);
            };
            auto ccsum = [&c, &cctrm](int k, int n) {
                real_t S = 0;
                for (int j = 1; j <= n / 2; j++)
                    S += cctrm(k, j, n);
                return c(k, n) / real_t(2 * n) * (real_t(1) - S);
            };

            nodes_.resize(s);
            weights_.resize(s);
            real_t si = 1 / real_t(s - 1);
            for (int k = 0; k < s; k++) {
                nodes_(k) = (1 - cos(pi_ * k * si)) / 2 ;
                weights_(k) = ccsum(k, s-1);
            }
        }
        VectorXr nodes() const { return nodes_; }
        VectorXr weights() const { return weights_; }
};

template <typename real_t>
const real_t ClenshawCurtisQuadrature<real_t>::pi_ = 4 * atan(real_t(1));

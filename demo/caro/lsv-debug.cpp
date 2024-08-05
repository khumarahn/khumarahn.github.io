//#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <numeric>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <iostream>

#include "lsv.h"

typedef double real_t;
const int PREC = 42;

typedef Eigen::Matrix<real_t,3,1>                           Vector3r;
typedef Eigen::Matrix<real_t,Eigen::Dynamic,1>              VectorXr;
typedef Eigen::Matrix<real_t,3,3>                           Matrix3r;
typedef Eigen::Matrix<real_t,Eigen::Dynamic,Eigen::Dynamic> MatrixXr;

typedef cheb_ns::Cheb<real_t> real_cheb_t;

typedef std::function<real_cheb_t(real_cheb_t)> cheb_operator_t;

template<typename var_t> var_t SQR(var_t x) {
    return x * x;
}
template<typename var_t> var_t QUB(var_t x) {
    return x * x * x;
}

typedef lsv_ns::LSV<real_t, PREC> BaseLSV;

class LSV : public BaseLSV {
    private:
        real_cheb_t h_cheb_, h_cheb_p_, h_cheb_pp_;
        MatrixXr R_;
    public:
        void set_gamma(double gamma) {
            BaseLSV::set_gamma(gamma);

            real_t a = 0.5, b = 1.0;

            R_ = Lind();
            Eigen::EigenSolver<MatrixXr> eigensolver(R_);

            VectorXr ev_abs = eigensolver.eigenvalues().cwiseAbs();

            std::vector<int> idx(R_.rows());
            std::iota(idx.begin(), idx.end(), 0);
            std::sort(idx.begin(), idx.end(),
                    [&ev_abs] (int i, int j) { return ev_abs(i) > ev_abs(j); });

            // approximation of invariant function
            real_cheb_t hn(eigensolver.eigenvectors().col(idx[0]).real(), a, b);
            // normalize so that integral on [0.5, 1.0] is one
            real_cheb_t hi = hn.integral();
            real_t norm = hi(real_t(1.0)) - hi(real_t(0.5));

            h_cheb_ = real_cheb_t(hn.coef() / norm, a, b);
            h_cheb_p_ = h_cheb_.derivative();
            h_cheb_pp_ = h_cheb_p_.derivative();
        }
        LSV() {
            set_gamma(1.0);
        }
        double gamma() const {
            return BaseLSV::gamma();
        }
        double h(real_t x) const {
            if (x < 1./128 )
                return 0;
            if (x > 1)
                return h(1);
            if (x >= 0.5)
                return h_cheb_(x);

            Vector2r z = left(x);
            Vector2r w = right_inv(z(0));

            return (h(z(0)) - h(w(0)) * w(1)) * z(1);
        }
        // inverse derivatives
        Vector3r www(real_t x) const {
            real_t tt = (gamma() + 1) * pow(2 * x, gamma()),
                   fp = 1 + tt,
                   fpp = gamma() * tt / x,
                   fppp = (gamma() - 1) * fpp / x,
                   w0 = 1 / fp,
                   w1 = - fpp / (fp*fp),
                   w2 = - fppp / (fp*fp) + 2 * (fpp*fpp) / (fp*fp*fp);
            return Vector3r(w0, w1, w2);
        }
        Vector3r h_full(const real_t &x) const {
            if (x < 1./128 - 1./1024) {
                return Vector3r(0,0,0);
            } else {
                Matrix3r factor = Matrix3r::Identity();
                Vector3r sum = Vector3r::Zero();

                real_t xx = x;

                while (xx < 0.5) {

                    Vector3r wx = www(xx);
                    real_t wy = 0.5;

                    real_t tx = left(xx)(0),
                           y = (tx + 1) / 2;

                    Vector3r hy = h_full(y);

                    Matrix3r A = Matrix3r::Zero();
                    Vector3r b = Vector3r::Zero();

                    A(0, 0) = 1 / wx(0);
                    b(0) = -wy / wx(0) * hy(0);

                    A(1, 1) = 1 / SQR(wx(0));
                    A(1, 0) = -wx(1) / wx(0) * A(0, 0);
                    b(1) = -SQR(wy) / SQR(wx(0)) * hy(1) - wx(1) / wx(0) * b(0);

                    A(2, 2) = 1 / QUB(wx(0));
                    A(2, 1) = -3 * wx(1) / wx(0) * A(1, 1);
                    A(2, 0) = -3 * wx(1) / wx(0) * A(1, 0) - (wx(2) / wx(0) + SQR(wx(1)) / SQR(wx(0))) * A(0, 0);
                    b(2) = -QUB(wy) / QUB(wx(0)) * hy(2) - 3 * wx(1) / wx(0) * b(1) -
                        (wx(2) / wx(0) + SQR(wx(1)) / SQR(wx(0))) * b(0);

                    // now h_full(xx) = A * h_full(tx) + b
                    sum += factor * b;
                    factor = factor * A;

                    xx = tx;
                }

                Vector3r hx(h_cheb_(xx), h_cheb_p_(xx), h_cheb_pp_(xx));
                return factor * hx + sum;
            }
        }
        real_t rho(real_t x) {
            return h_full(x)(0);
        }
        real_t rho_p(real_t x) {
            return h_full(x)(1);
        }
        real_t rho_pp(real_t x) {
            return h_full(x)(2);
        }
};

int main() {
    LSV lsv;

    using std::cout;

    for (double alpha = 0.25; alpha <= 4.0; alpha += 0.125) {
        lsv.set_gamma(1. / alpha);
        cout << "gamma: " << lsv.gamma()
            << "\n"
            << "  h(1/128) / h(1): " << (lsv.h_full(1./128).array() / lsv.h_full(1.).array()).transpose()  << "\n"
            << "\n";
    }

    return 0;
}

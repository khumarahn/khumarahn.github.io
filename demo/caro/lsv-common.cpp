// This extends the base LSV class with the computation of the invariant
// density and other things of interest

#include <numeric>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

typedef Eigen::Matrix<real_t,3,1>                           Vector3r;
typedef Eigen::Matrix<real_t,Eigen::Dynamic,1>              VectorXr;
typedef Eigen::Matrix<real_t,3,3>                           Matrix3r;
typedef Eigen::Matrix<real_t,Eigen::Dynamic,Eigen::Dynamic> MatrixXr;

typedef std::complex<real_t> complex_t;
typedef Eigen::Matrix<complex_t,Eigen::Dynamic,1>           VectorXc;

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
        VectorXc R_evalues_;
    public:
        void set_gamma(real_t gamma) {
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

            // save eigenvalues sorted in decreasing order
            R_evalues_.resize(R_.rows());
            for (int j  = 0; j < R_.rows(); j++)
                R_evalues_(j) = eigensolver.eigenvalues()(idx[j]);
        }
        LSV() {
            set_gamma(1.0);
        }
        real_t gamma() const {
            return BaseLSV::gamma();
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
        //
        Vector3r f_abel(real_t x) {
            const int n = Nstar();

            // Let y = f^{-n}(x) on the left branch.
            // p = (f^{-n})'(x) and q = (f^{-n})''(x)
            real_t y = x, p = 1.0, q = 0.0;

            for (int j = 0; j < n; j++) {
                Vector2r yy = left_inv(y);
                real_t w = yy(0);       // new y
                real_t h_p = yy(1);     // h'(old_y) = 1 / f'(w)

                // f''(w) = (gamma + 1) * gamma * (2w)^(gamma - 1) * 2
                real_t f_pp = (gamma() + 1) * gamma() * pow(2 * w, gamma() - 1) * 2;

                // h''(old_y) = -f''(w) / (f'(w))^3 = -f''(w) * (h_p)^3
                real_t h_pp = -f_pp * h_p * h_p * h_p;

                // Chain rule for the next step:
                // y_{new}'' = h''(y_{old}) * (y_{old}')^2 + h'(y_{old}) * y_{old}''
                q = h_pp * p * p + h_p * q;

                // y_{new}' = h'(y_{old}) * y_{old}'
                p *= h_p;

                y = w;
            }

            Vector3r A = abel(y);

            return Vector3r(
                    A(0) - n,
                    A(1) * p,
                    A(2) * p * p + A(1) * q
            );
        }

        real_t full_abel(real_t x) {
            return f_abel(x)(0);
        }
        real_t full_abel_p(real_t x) {
            return f_abel(x)(1);
        }
        real_t full_abel_pp(real_t x) {
            return f_abel(x)(2);
        }
        //
        Vector2r f_abel_inv(real_t z) {
            const int n = Nstar();
            real_t zn = z + n;

            Vector2r AI = abel_inv(zn);

            real_t y = AI(0), p = 1;
            for (int j = 0; j < n; j++) {
                Vector2r yy = left(y);
                y = yy(0);
                p *= yy(1);
            }

            return Vector2r(
                    y,
                    AI(1) * p
                    );
        }
        real_t full_abel_inv(real_t x) {
            return f_abel_inv(x)(0);
        }
        real_t full_abel_inv_p(real_t x) {

            return f_abel_inv(x)(1);
        }
        //
        // Combined method for w(x) and w'(x) to reuse Abel evaluations
        Vector2r w_full(real_t x) {
            // y = A^{-1}(x)
            real_t y = full_abel_inv(x);
            // z = 2y - 1
            real_t z = 2 * y - 1;

            // Evaluate A, A', A'' at y and z
            Vector3r Ay = f_abel(y);
            Vector3r Az = f_abel(z);

            real_t Ay_p = Ay(1);
            real_t Ay_pp = Ay(2);

            real_t Az_p = Az(1);
            real_t Az_pp = Az(2);

            // w(x) = A'(y) / (2 * A'(z))
            real_t w_val = Ay_p / (2.0 * Az_p);

            // w'(x) = (A''(y)*A'(z) - 2*A'(y)*A''(z)) / (2 * A'(y) * (A'(z))^2)
            real_t w_deriv = (Ay_pp * Az_p - 2 * Ay_p * Az_pp) / (2 * Ay_p * Az_p * Az_p);

            return Vector2r(w_val, w_deriv);
        }

        real_t w(real_t x) {
            return w_full(x)(0);
        }

        real_t w_p(real_t x) {
            return w_full(x)(1);
        }
        real_t h(real_t x) const {
            return h_full(x)(0);
        }
        real_t h_p(real_t x) const {
            return h_full(x)(1);
        }
        real_t h_pp(real_t x) const {
            return h_full(x)(2);
        }
        real_t h_coef(int k) const {
            return h_cheb_.coef(k);
        }
        int R_cols() const {
            return R_.cols();
        }
        real_t R_coef(int i, int j) const {
            return R_(i, j);
        }
        complex_t R_evalues(int i) const {
            return R_evalues_(i);
        }
        real_t R_evalues_real(int i) const {
            return R_evalues_(i).real();
        }
        real_t R_evalues_imag(int i) const {
            return R_evalues_(i).imag();
        }
};

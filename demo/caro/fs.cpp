//#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>
#include <format>
#include <numeric>

#include <Eigen/Dense>

#include <boost/math/tools/roots.hpp>
#include <tuple>

#include "cheb.h"
#include "clenshawcurtis.h"


using real_t = long double;
const int PREC = 64;
const int real_digits = std::numeric_limits<real_t>::digits;

using real_func_t = std::function<real_t (real_t)>;
using VectorXr = Eigen::Matrix<real_t, Eigen::Dynamic, 1>;
using MatrixXr = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

using Cheb = cheb_ns::Cheb<real_t>;

const int CHEB_N = 128;

const ClenshawCurtisQuadrature<real_t> clenshaw(CHEB_N);

const real_t pi = 4 * std::atan(real_t(1));


// find a root of f(y) = value, when f and its derivative are given by Cheb
// using Newton-Raphson in interval [low,high]
real_t cheb_root(const Cheb &f, const Cheb &fp,
        const real_t &low, const real_t &high, const real_t &value) {
    auto FF = [&f, &fp, &value] (real_t y) {
        double fy = f.value(y) - value;
        double dfy = fp.value(y);
        return std::make_tuple(fy, dfy);
    };

    real_t guess = (low + high) / 2;
    real_t root = boost::math::tools::newton_raphson_iterate(
            FF, guess, low, high, real_digits);
    return root;
}

// transfer operator of a map f:[0,1] -> [0,2] (assumed we do mod 1)
// when f and its derivative are Cheb
real_t L(const Cheb &f, const Cheb &fp, const real_func_t &v, const real_t &y) {
    real_t s = cheb_root(f, fp, 0.0, 1.0, y);
    real_t t = cheb_root(f, fp, 0.0, 1.0, y + 1);
    return v(s) / std::abs(fp.value(s)) + v(t) / std::abs(fp.value(t));
}
real_t L(const Cheb &f, const Cheb &fp, const Cheb &v, const real_t &y) {
    auto vf = [&v] (const real_t y) {
        return v.value(y);
    };
    return L(f, fp, vf, y);
}

// ... a matrix representation of the transfer operator
MatrixXr L_mat(const Cheb &f, const Cheb &fp) {
    Cheb cheb(0.0, 1.0, CHEB_N);

    const VectorXr y_nodes = cheb.nodes();


    MatrixXr L_values(CHEB_N, CHEB_N);
    for (int iy = 0; iy < y_nodes.size(); iy++) {
        real_t y = y_nodes[iy];
        VectorXr r = VectorXr::Zero(CHEB_N);

        // evaluation of all chebs
        auto v = [&cheb, N = CHEB_N] (const real_t &y) {
            return cheb.basis_values(y, N);
        };

        // loop over the two branches
        for (int b = 0; b <= 1; b++) {
            real_t z = cheb_root(f, fp, 0.0, 2.0, y + b);
            real_t fpz = std::abs(fp.value(z));

            r += v(z) / fpz;
        }
        L_values.col(iy) = r;
    }

    MatrixXr R(CHEB_N, CHEB_N);
    for (int k = 0; k < CHEB_N; k++) {
        Cheb c(0.0, 1.0, CHEB_N);
        c.set_from_values(L_values.row(k).transpose());
        R.col(k) = c.coef();
    }

    return R;
}

// ... density of invariant probability measure
Cheb inv_density(const Cheb &f, const Cheb &fp) {
    auto R = L_mat(f, fp);
    Eigen::EigenSolver<MatrixXr> eigensolver(R);

    VectorXr ev_abs = eigensolver.eigenvalues().cwiseAbs();

    std::vector<int> idx(R.rows());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
            [&ev_abs] (int i, int j) { return ev_abs(i) > ev_abs(j); });

    // approximation of invariant function
    Cheb hn(eigensolver.eigenvectors().col(idx[0]).real(), 0.0, 1.0);
    // normalize so that the integral is 1
    Cheb hi = hn.integral();
    real_t norm = hi(real_t(1.0)) - hi(real_t(0.0));

    Cheb h(hn.coef() / norm, 0.0, 1.0);

    return h;
}

Cheb inv_density(real_func_t f) {
    Cheb f_cheb(f, 0.0, 1.0, CHEB_N);
    return inv_density(f_cheb, f_cheb.derivative());
}

// integrate a function on [0,1] with Clenshaw-Curtis
real_t integrate(real_func_t v) {
    const VectorXr &nodes = clenshaw.nodes();
    const VectorXr &weights = clenshaw.weights();
    real_t r = 0;
    for (int k = 0; k < CHEB_N; k++) {
        r += v(nodes[k]) * weights[k];
    }
    return r;
}

class fs_t {
    private:
        // this is v before centering
        // it has to be periodic in x with period [-pi,pi],
        // so that we can reliably compute and subtract averages
        real_t v_(real_t x, real_t y) const {
            return std::sin(x + y);
        }
        Cheb v_mean_;
    public:
        fs_t() {
            auto v_m = [this] (real_t x) {
                Cheb id = inv_density([this, x](real_t y) { return f(x, y); } );
                return integrate([this, &x, &id](real_t y) { return v_(x, y) * id.value(y); });
            };
            v_mean_ = Cheb(v_m, -pi, pi, CHEB_N);
        }
        real_t v(real_t x, real_t y) const {
            real_t xx = std::fmod(x + pi, 2 * pi) - pi;
            return v_(xx, y) - v_mean_.value(xx);
        }
        // f has to be smooth, expanding, have TWO branches,
        // necessarily f(*,0) = 0, f(*,1) = 2)
        real_t f(real_t x, real_t y) const {
            return 2 * y + 0.1 * y * (1 - y) * std::sin(x);
        }

};

class std_pair_t {
    public:
        Cheb x, rho;
        std_pair_t(real_func_t x_, real_func_t rho_) {
            x = Cheb(x_, 0.0, 1.0, CHEB_N);
            rho = Cheb(rho_, 0.0, 1.0, CHEB_N);
        }
        std_pair_t() : std_pair_t([](real_t y) { return 0; }, [](real_t y) { return 1; }) {}
};

using std_fam_t = std::vector<std_pair_t>;

std_fam_t fam_evolve(const std_fam_t &fam, const fs_t &fs, const real_t &eps) {
    std_fam_t r;
    const real_t epsq = std::sqrt(eps);
    for (const auto &p : fam) {
        auto ff = [&p, &fs] (real_t y) {
            return fs.f(p.x.value(y), y);
        };
        Cheb f(ff, 0.0, 1.0, CHEB_N);
        Cheb fp = f.derivative();
        auto vv = [&p, &fs] (real_t y) {
            return fs.v(p.x.value(y), y);
        };

        // loop over branches of f
        for (int b = 0; b <= 1; b++) {
            auto finv = [&f, &fp, &b] (real_t y) {
                return cheb_root(f, fp, 0.0, 1.0, y + b);
            };
            auto x_ = [&p, &vv, &finv, &epsq] (real_t y) {
                real_t yi = finv(y);
                return p.x.value(yi) + epsq * vv(yi);
            };
            auto rho_ = [&p, &fp, &finv] (real_t y) {
                real_t yi = finv(y);
                return p.rho.value(yi) / fp(yi);
            };

            r.push_back(std_pair_t(x_, rho_));
        }
    }
    return r;
}

std_fam_t fam_evolve(const std_pair_t &p, const fs_t &fs, const real_t &eps) {
    std_fam_t fam;
    fam.push_back(p);

    return fam_evolve(fam, fs, eps);
}

int main() {
    // construct a Chebyshev approximation and try to invert it!

    using std::cout;

    auto f = [] (real_t x) -> real_t {
        return 1.5 * x + 0.5 * x*x;
    };

    cheb_ns::Cheb<real_t> f_cheb(f, 0.0, 1.0, CHEB_N);
    auto fp_cheb = f_cheb.derivative();

    real_t max_err = 0.0;
    for (real_t val = 0.5; val < 0.6; val += 0.0001) {
        real_t root = cheb_root(f_cheb, fp_cheb, 0.0, 1.0, val);
        max_err = std::max(max_err, std::abs(f(root) - val));
    }

    cout << "root max err: " << max_err
        << "\n";

    auto h = inv_density(f_cheb, fp_cheb);

    max_err = 0.0;
    for (real_t y = 0.0; y <= 1.0; y += 1./1024.) {
        max_err = std::max(max_err, std::abs(h.value(y) - L(f_cheb, fp_cheb, h, y) ));
    }
    cout << "invariant density max err: " << max_err
        << "\n";

    real_t idi = integrate([&h](real_t y) { return h.value(y); });
    cout << "integral of the invariant density - 1: " << idi - 1.0
        << "\n";

    fs_t fs;
    /*
    max_err = 0;
    for (real_t x = 4.0; x <= 6.0; x += 1. / 128.) {
        auto ds = [&fs, x](real_t y) {
            return fs.f(x, y);
        };
        auto hh = inv_density(ds);
        real_t vi = integrate([&hh, &x, &fs](real_t y) { return hh.value(y) * fs.v(x,y); });
        max_err = std::max(max_err, std::abs(vi));
    }
    cout << "error of averages of v in fs: " << max_err
        << "\n";
    */

    for (real_t eps = 1./8.; eps >= 1./1000000.; eps *= 0.75) {
        std_fam_t fam(1);
        for (int k = 0; k < 8; k++) {
            fam = fam_evolve(fam, fs, eps);
        }
        real_t W = 0, D = 0;
        for (const auto &p : fam) {
            auto w = [&p] (real_t y) {
                return p.rho.value(y);
            };
            auto d = [&p] (real_t y) {
                return p.rho.value(y) * p.x.value(y);
            };
            W += integrate(w);
            D += integrate(d);
        }
        cout << std::format("eps: {:11.8f}, W-1: {:11.8f}, D: {:11.8f}, D/eps: {:11.8f}\n",
                eps, W-1, D, D/eps);
    }


    return 0;
}

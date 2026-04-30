// Interval root finding routines
//
// Designed to work with these basic types
//      using real_t     = boost::multiprecision::number<bmp::mpfr_float_backend<DIGITS>>;
//      using interval_t = boost::multiprecision::number<bmp::mpfi_float_backend<DIGITS>>;

#pragma once

#include <complex>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpfi.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/eigen.hpp>


namespace interval_root_ns {

namespace bmp = boost::multiprecision;

// TYPES

template <unsigned d10>
    using real_t     = boost::multiprecision::number<bmp::mpfr_float_backend<d10>>;
template <unsigned d10>
    using interval_t = boost::multiprecision::number<bmp::mpfi_float_backend<d10>>;

template <unsigned d10>
    using complex_t           = std::complex<real_t<d10>>;
template <unsigned d10>
    using complex_interval_t  = std::complex<interval_t<d10>>;

template <unsigned d10>
    using Vector2r = Eigen::Matrix<real_t<d10>, 2, 1>;
template <unsigned d10>
    using Vector2i = Eigen::Matrix<interval_t<d10>, 2, 1>;
template <unsigned d10>
    using Vector2c = Eigen::Matrix<complex_t<d10>, 2, 1>;
template <unsigned d10>
    using Vector2ci = Eigen::Matrix<complex_interval_t<d10>, 2, 1>;
template <unsigned d10>
    using VectorXr = Eigen::Matrix<real_t<d10>, Eigen::Dynamic, 1>;
template <unsigned d10>
    using VectorXi = Eigen::Matrix<interval_t<d10>, Eigen::Dynamic, 1>;
template <unsigned d10>
    using MatrixXr = Eigen::Matrix<real_t<d10>, Eigen::Dynamic, Eigen::Dynamic>;
template <unsigned d10>
    using MatrixXi = Eigen::Matrix<interval_t<d10>, Eigen::Dynamic, Eigen::Dynamic>;

// HELPER FUNCTIONS FOR COMPLEX INTERVALS

template <unsigned d10>
complex_t<d10> median(const complex_interval_t<d10> &x) {
    return complex_t<d10>(bmp::median(x.real()), bmp::median(x.imag()));
}

template <unsigned d10>
complex_interval_t<d10> intersect(const complex_interval_t<d10> &a, const complex_interval_t<d10> &b) {
    return complex_interval_t<d10>(
            bmp::intersect(a.real(), b.real()),
            bmp::intersect(a.imag(), b.imag())
            );
}

template <unsigned d10>
bool overlap(const complex_interval_t<d10> &a, const complex_interval_t<d10> &b) {
    return bmp::overlap(a.real(), b.real()) && bmp::overlap(a.imag(), b.imag());
}

template <unsigned d10>
bool subset(const complex_interval_t<d10> &a, const complex_interval_t<d10> &b) {
    return bmp::subset(a.real(), b.real()) && bmp::subset(a.imag(), b.imag());
}

template <unsigned d10>
bool proper_subset(const complex_interval_t<d10> &a, const complex_interval_t<d10> &b) {
    return
        ( bmp::proper_subset(a.real(), b.real()) && bmp::subset(a.imag(), b.imag()) ) ||
        ( bmp::subset(a.real(), b.real()) && bmp::proper_subset(a.imag(), b.imag()) ) ;
}

// HELPER FUNCTIONS FOR INTERVAL VECTORS

template <unsigned d10>
VectorXi<d10> intersect(const VectorXi<d10> &a, const VectorXi<d10> &b) {
    assert(a.size() == b.size());
    const int N = a.size();
    VectorXi<d10> r(N);
    for (int i=0; i<N; i++) {
        r(i) = bmp::intersect(a(i),b(i));
    }
    return r;
}

template <unsigned d10>
bool subset(const VectorXi<d10> &a, const VectorXi<d10> &b) {
    assert(a.size() == b.size());
    const int N = a.size();
    for (int i=0; i<N; i++) {
        if (!bmp::subset(a(i),b(i))) {
            return false;
        }
    }
    return true;
}

template <unsigned d10>
bool proper_subset(const VectorXi<d10> &a, const VectorXi<d10> &b) {
    assert(a.size() == b.size());
    bool proper = false;
    const int N = a.size();
    for (int i=0; i<N; i++) {
        if (!bmp::subset(a(i),b(i))) {
            return false;
        }
        if (!proper && bmp::proper_subset(a(i),b(i))) {
            proper = true;
        }
    }
    return proper;
}

// Univariate interval Newton method. Requires correct input.
// The first argument is a function
//      interval_t -> Vector2<interval_t>
// where the 0th component of the vector is the function value,
// and the 1st component is the derivative.
// The second argument is an interval containing a root.
//
// We terminate when the interval stops changing. We return the solution
// and the inverse derivative.
template <unsigned d10, typename f_t>
Vector2i<d10> interval_newton(const f_t &f, const interval_t<d10> &guess) {
    using i_t = interval_t<d10>;
    using V2i_t = Vector2i<d10>;

    i_t Y = guess,
        P;

    for (;;) {
        V2i_t fY = f(Y);
        P = fY(1);

        i_t m = bmp::median(Y);
        i_t Z = m - f(m)(0) / P;
        i_t ZZ = bmp::intersect(Z, Y); // is this needed?

        if (!bmp::proper_subset(ZZ,Y))
            break;

        Y = ZZ;
    }

    // Sanity check: image of the returned interval contains zero
    assert(bmp::subset(i_t(0), f(Y)(0)));

    return V2i_t(Y, 1 / P);
}

// Krawczyk method for a complex function complex_krawczyk(f, guess).
// Here f is a function
//      complex_interval_t -> Vector2<complex_interval_t>
// where the 0th component of the vector is the function value,
// and the 1st component is the derivative.
//
// guess is an interval containing a root.
//
// We terminate when the interval stops changing. We return the solution
// and the inverse derivative.
template <unsigned d10, typename f_t>
Vector2ci<d10> complex_krawczyk(const f_t &f, const complex_interval_t<d10> &guess) {
    using c_t = complex_t<d10>;
    using ci_t = complex_interval_t<d10>;
    using V2ci_t = Vector2ci<d10>;

    ci_t X = guess;

    V2ci_t fX;
    for (;;) {
        fX = f(X);
        ci_t y = median(X),
             Z = X - y,
             fpX = fX(1),
             Y = c_t(1) / median(fpX), // can be optimized
             K = y - Y * f(y)(0) + (ci_t(1) - Y * fpX) * Z;

        assert(overlap(X, K));
        ci_t XX = intersect(X, K);

        if (!proper_subset(XX,X)) {
            break;
        }

        X = XX;
    }

    assert(subset(ci_t(0), fX(0)));

    return V2ci_t(X, ci_t(1) / fX(1));
}

// Krawczyk method for solving a real linear equation Ax = b
//
// We terminate when the interval stops changing.
template <unsigned d10>
VectorXi<d10> linear_krawczyk(const MatrixXi<d10> &A, const VectorXi<d10> &b) {
    using r_t = real_t<d10>;
    using i_t = interval_t<d10>;
    using M_t = MatrixXr<d10>;
    using Mi_t = MatrixXi<d10>;
    using Vi_t = VectorXi<d10>;

    assert(A.rows() == A.cols());
    const int N = A.cols();

    // norm
    auto norm = [] (const Mi_t &U) {
        // abs for intervals
        auto iabs = [] (const auto &x) {
            return max(abs(bmp::lower(x)), abs(bmp::upper(x)));
        };
        r_t no = 0;
        for (int i=0; i<U.rows(); i++) {
            r_t s = 0;
            for (int j=0; j<U.cols(); j++) {
                s += iabs(U(i,j));
            }
            no = max(no,s);
        }
        return no;
    };

    // Y is an approximate inverse of the median of A
    Mi_t Y(N,N);
    {
        M_t Y2(N,N);
        for (int j=0; j<N; j++)
            for (int i=0; i<N; i++)
                Y2(i,j) = bmp::median(A(i,j));
        M_t Y1 = Y2.inverse();
        for (int j=0; j<N; j++)
            for (int i=0; i<N; i++)
                Y(i,j) = Y1(i,j);
    }

    Mi_t E = Mi_t::Identity(N,N) - Y * A;
    r_t normE = norm(E);
    assert(normE < 1);

    i_t x0 = i_t(-1,1) * norm(Y * b) / (1 - normE);
    Vi_t X = Vi_t::Constant(N, x0);

    // check that the initial interval X makes sense
    assert(subset(b, Vi_t(A * X)));

    for (;;) {
        Vi_t XX = Y * b + E * X,
             XXX = intersect(X, XX);

        if (!proper_subset(XXX,X)) {
            break;
        }
        X = XXX;
    }

    //std::cout << "Linear Kraw, error norm: " << norm(A * X - b) << "\n";

    // check that the initial interval X makes sense
    assert(subset(b, Vi_t(A * X)));

    return X;
}

// Given a matrix, return an interval enclosing all singular values, if only_top = false,
// and only the top singular value if only_top = true
template <unsigned d10>
interval_t<d10> matrix_sigma_bounds(const MatrixXi<d10> &W, bool only_top = false) {
    using i_t = interval_t<d10>;
    using Mi_t = MatrixXi<d10>;

    assert(W.rows() == W.cols());
    const int n = W.rows();

    Mi_t WtW = W.transpose() * W;
    i_t eig;

    for (int i = 0; i < n; i++) {
        i_t center = WtW(i,i);
        i_t radius = 0;
        for (int j = 0; j < n; j++)
            if (i != j)
                radius += bmp::upper(i_t(abs(WtW(i,j))));
        i_t e = center + i_t(-bmp::upper(radius), bmp::upper(radius));

        if (i == 0) {
            eig = e;
        } else if (!only_top) {
            eig = bmp::hull(eig, e);
        } else if (only_top) {
            eig = i_t(
                    max(bmp::lower(eig), bmp::lower(center)),
                    max(bmp::upper(eig), bmp::upper(e))
                    );
        }
    }

    if (bmp::lower(eig) < 0)
        eig = i_t(0, bmp::upper(eig));

    return sqrt(eig);
}

// A basic rigorous estimate or the top singular value of a square
// interval matrix, using a floating-point SVD and Gershgorin discs
template <unsigned d10>
interval_t<d10> matrix_L2_norm(const MatrixXi<d10> &Q) {
    using i_t = interval_t<d10>;
    using M_t = MatrixXr<d10>;
    using Mi_t = MatrixXi<d10>;

    assert(Q.rows() == Q.cols());
    const int n = Q.rows();

    // median
    M_t M(n, n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            M(i,j) = bmp::median(Q(i,j));

    // floating point SVD
    Eigen::JacobiSVD<M_t> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // M = U D V^*
    M_t U_fp = svd.matrixU();
    M_t V_fp = svd.matrixV();
    //std::cout << "****\n"
    //    << (U_fp.transpose() * M * V_fp).topLeftCorner(5,5)
    //    << "\n****\n";

    // interval versions of U, V
    Mi_t U(n, n), V(n, n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            U(i,j) = i_t(U_fp(i,j));
            V(i,j) = i_t(V_fp(i,j));
        }
    }

    // an approximate SVD in interval arithmetic:
    // Q contained in (U^*)^{-1} S V^{-1}, where S = U^* Q V
    Mi_t S = U.transpose() * Q * V;

    // now, the max singular value of Q is in the interval
    // sigma_U^{-1} * top_sigma_S * sigma_V^{-1},
    // where sigma_U, sigma_V are intervals enclosing singular values of U,V,
    // and top_sigma_S encloses the top singular value of S, all
    // computed simply with Gershgorin discs

    i_t sigma_U = matrix_sigma_bounds(U),
          sigma_V = matrix_sigma_bounds(V),
          sigma_S = matrix_sigma_bounds(S, true);

    assert(sigma_U > 0 && sigma_V > 0);

    //std::cout << "sigma_U: " << sigma_U << ", sigma_V: " << sigma_V << ", sigma_S: " << sigma_S << "\n"
    //    << S.topLeftCorner(5,5) << "\n";

    return sigma_S / (sigma_U * sigma_V);
}

} // namespace

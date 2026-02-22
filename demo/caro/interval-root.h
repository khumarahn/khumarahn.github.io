// Interval root finding routines
//
// Designed to work with these basic types
//      using real_t     = boost::multiprecision::number<bmp::mpfr_float_backend<DIGITS>>;
//      using interval_t = boost::multiprecision::number<bmp::mpfi_float_backend<DIGITS>>;

#pragma once

#include <complex>
#include <Eigen/Core>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpfi.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/eigen.hpp>


namespace interval_root_ns {

namespace bmp = boost::multiprecision;

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

template <unsigned d10>
complex_t<d10> median(const complex_interval_t<d10> &x) {
    return complex_t<d10>(bmp::median(x.real()), bmp::median(x.imag()));
}

template <unsigned d10>
real_t<d10> width(const complex_interval_t<d10> &x) {
    return max(bmp::width(x.real()), bmp::width(x.imag()));
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
    bool proper = false;
    assert(a.size() == b.size());
    const int N = a.size();
    for (int i=0; i<N; i++) {
        if (!bmp::subset(a(i),b(i))) {
            return false;
        }
        if (bmp::proper_subset(a(i),b(i))) {
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
// The second component is an interval containing a root.
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

        if (!bmp::proper_subset(Z,Y))
            break;

        Y = Z;
    }

    // Sanity check: image of the returned interval contains zero
    assert(bmp::subset(i_t(0), f(Y)(0)));

    return V2i_t(Y, 1 / P);
}

// Krawczyk for a complex function
template <unsigned d10, typename f_t>
Vector2ci<d10> complex_krawczyk(const f_t &f, const complex_interval_t<d10> &guess) {
    using c_t = complex_t<d10>;
    using ci_t = complex_interval_t<d10>;
    using V2ci_t = Vector2ci<d10>;

    ci_t X = guess;

    int k;
    V2ci_t fX;
    for (k = 0; k < 1024; k++) {
        fX = f(X);
        ci_t y = median(X),
             Z = X - y,
             fpX = fX(1),
             Y = c_t(1) / median(fpX), // can be optimized
             K = y - Y * f(y)(0) + (ci_t(1) - Y * fpX) * Z;

        assert(overlap(X, K));
        complex_interval_t XX = intersect(X, K);

        if (!proper_subset(XX,X)) {
            break;
        }

        X = XX;
    }

    assert(k < 64);
    assert(subset(ci_t(0), fX(0)));
    //assert(width(X) < real_eps_sqrt_);

    return V2ci_t(X, ci_t(1) / fX(1));
}

// Krawczyk method for solving a real linear equation Ax = b
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
        real_t<d10> no = 0;
        for (int i=0; i<U.rows(); i++) {
            real_t<d10> s = 0;
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

    for (int k = 0; k < 1024; k++) {
        Vi_t XX = Y * b + E * X;

        X = intersect(X, XX);

        if (!proper_subset(XX,X)) {
            break;
        }
        X = XX;

        assert(k < 64);
    }

    //std::cout << "Linear Kraw, error norm: " << norm(A * X - b) << "\n";

    // check that the initial interval X makes sense
    assert(subset(b, Vi_t(A * X)));

    return X;
}


}

// black magic which we need when compiling for the web

// bug in boost
#define variable_precision_optionss variable_precision_options

#include <iostream>

namespace std {
    template <typename T> requires std::is_class_v<T> T abs(const T& x);
    template <typename T> requires std::is_class_v<T> T sqrt(const T& x);
    template <typename T> requires std::is_class_v<T> T exp(const T& x);
    template <typename T> requires std::is_class_v<T> T log(const T& x);
    template <typename T> requires std::is_class_v<T> T sin(const T& x);
    template <typename T> requires std::is_class_v<T> T cos(const T& x);
    template <typename T> requires std::is_class_v<T> T sinh(const T& x);
    template <typename T> requires std::is_class_v<T> T cosh(const T& x);
    template <typename T> requires std::is_class_v<T> T atan2(const T& y, const T& x);
    template <typename T> requires std::is_class_v<T> T hypot(const T& x, const T& y);
    template <typename T> requires std::is_class_v<T> T copysign(const T& x, const T& y);
    template <typename T> requires std::is_class_v<T> bool isfinite(const T& x);
    template <typename T> requires std::is_class_v<T> bool isinf(const T& x);
    template <typename T> requires std::is_class_v<T> bool isnan(const T& x);
    template <typename T> requires std::is_class_v<T> bool signbit(const T& x);
}

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpfi.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <Eigen/Dense>

namespace std {
    template <typename T> requires std::is_class_v<T>
    T sqrt(const T& x) { return boost::multiprecision::sqrt(x); }

    template <typename T> requires std::is_class_v<T>
    T abs(const T& x) { return boost::multiprecision::abs(x); }

    template <typename T> requires std::is_class_v<T>
    T exp(const T& x) { return boost::multiprecision::exp(x); }

    template <typename T> requires std::is_class_v<T>
    T log(const T& x) { return boost::multiprecision::log(x); }

    template <typename T> requires std::is_class_v<T>
    T sin(const T& x) { return boost::multiprecision::sin(x); }

    template <typename T> requires std::is_class_v<T>
    T cos(const T& x) { return boost::multiprecision::cos(x); }

    template <typename T> requires std::is_class_v<T>
    T sinh(const T& x) { return boost::multiprecision::sinh(x); }

    template <typename T> requires std::is_class_v<T>
    T cosh(const T& x) { return boost::multiprecision::cosh(x); }

    template <typename T> requires std::is_class_v<T>
    T atan2(const T& y, const T& x) { return boost::multiprecision::atan2(y, x); }

    template <typename T> requires std::is_class_v<T>
    T hypot(const T& x, const T& y) { return boost::multiprecision::sqrt(x*x + y*y); }

    template <typename T> requires std::is_class_v<T>
    T copysign(const T& x, const T& y) { return x; }

    template <typename T> requires std::is_class_v<T>
    bool isfinite(const T&) { return true; }

    template <typename T> requires std::is_class_v<T>
    bool isinf(const T&) { return false; }

    template <typename T> requires std::is_class_v<T>
    bool isnan(const T&) { return false; }

    template <typename T> requires std::is_class_v<T>
    bool signbit(const T&) { return false; }
}

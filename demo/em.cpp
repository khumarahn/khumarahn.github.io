#include <vector>
#include <emscripten/bind.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using std::vector;
using namespace Eigen;
using namespace emscripten;

vector<double> do_something_mad(double p) {
    // Arnold's cat map
    Matrix2d m {
        {2, 1},
        {1, 1}
    };
    
    // compute mp = m^p
    MatrixPower<Matrix2d> m_pow(m);
    Matrix2d mp = m_pow(p);

    // convert to a vector and return
    vector<double> mp_as_vector(mp.data(), mp.data() + mp.size());
    return mp_as_vector;
}

EMSCRIPTEN_BINDINGS(my_module) {
    function("do_something_mad", &do_something_mad);
    register_vector<double>("vector<double>");
}

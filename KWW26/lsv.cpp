#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <iostream>
#include <string>
#include <stdexcept>

#include "web.h"

#include <emscripten/bind.h>

#include "lsv.h"

const int PREC = 32;

#include "lsv-common.cpp"

EMSCRIPTEN_BINDINGS(my_module) {
    emscripten::class_<LSV>("LSV")
        .constructor()
        .function("set_alpha", &LSV::set_alpha)
        .function("compute_L", &LSV::compute_L)
        .function("compute_h_meta", &LSV::compute_h_meta)
        .function("compute_h_cheb", &LSV::compute_h_cheb)
        .function("compute_F", &LSV::compute_F)
        .function("compute_derivative_signs_right", &LSV::compute_derivative_signs_right)
        .function("compute_derivative_bounds", &LSV::compute_derivative_bounds)
        .function("compute_tau", &LSV::compute_tau)
        .function("compute_lambda", &LSV::compute_lambda)
        .function("oracle", &LSV::oracle)
        .function("double_gamma", &LSV::double_gamma)
        .function("double_alpha", &LSV::double_alpha);
}

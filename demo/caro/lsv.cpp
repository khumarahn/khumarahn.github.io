#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <emscripten/bind.h>

#include "lsv.h"

typedef double real_t;
const int PREC = 52;

#include "lsv-common.cpp"

EMSCRIPTEN_BINDINGS(my_module) {
    emscripten::class_<LSV>("LSV")
        .constructor()
        .function("gamma", &LSV::gamma)
        .function("set_gamma", &LSV::set_gamma)
        .function("h", &LSV::h)
        .function("h_p", &LSV::h_p)
        .function("h_pp", &LSV::h_pp)
        .function("R_evalues_real", &LSV::R_evalues_real)
        .function("R_evalues_imag", &LSV::R_evalues_imag)
        .function("R_size", &LSV::R_size)
        .function("R_coef", &LSV::R_coef)
        .function("h_coef", &LSV::h_coef);
}

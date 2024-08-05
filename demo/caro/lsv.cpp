#define NDEBUG // disable Eigen's range checking, asserts and other good things

#include <emscripten/bind.h>

#include "lsv.h"

typedef double real_t;
const int PREC = 42;

#include "lsv-common.cpp"

EMSCRIPTEN_BINDINGS(my_module) {
    emscripten::class_<LSV>("LSV")
        .constructor()
        .function("gamma", &LSV::gamma)
        .function("set_gamma", &LSV::set_gamma)
        .function("h", &LSV::h)
        .function("h_p", &LSV::h_p)
        .function("h_pp", &LSV::h_pp);
}

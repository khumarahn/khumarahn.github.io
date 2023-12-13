/*
    To be compiled into javascript/wasm and ran on a webpage
*/

#include <iostream>
#include <string>
#include <sstream>
#include <emscripten/bind.h>
#include "hardballs.h"

using namespace emscripten;

using hardballs::f_type;
using hardballs::f_rand;
using hardballs::f_rand_norm;
using hardballs::table_c;
using hardballs::long_time_t;
using hardballs::SQR;

table_c<hardballs::TORUS> table;

std::stringstream table_log_;

std::string table_log() {
    return table_log_.str();
}

typedef struct {
    double x,y,vx,vy,r,m;
    int infected;
} particle;

std::vector<int> infecta;

f_type bottom_wall = 0.25;
f_type top_wall = 0.75;

std::vector<particle> hist_po_;
int hist_po_N() {
    return hist_po_.size();
}
particle hist_po(int n) {
    return hist_po_[n];
}

long num_col_pp = 0;
long table_num_col_pp() {
    return num_col_pp;
}
long num_col_po = 0;
long table_num_col_po() {
    return num_col_po;
}
long table_num_col() {
    return num_col_pp + num_col_po;
}

//#define heavy0
#ifdef heavy0
f_type test_c_r = 0.025,
       test_c_x = 0.5,
       test_c_y = 0.5;
#endif

void table_init(int N, double r, char col_rule_bottom, char col_rule_top, double gravity_x, double gravity_y, double alpha_bottom, double beta_bottom, double alpha_top, double beta_top) {

    table.SetTime(long_time_t(0.));
    num_col_pp = 0;
    num_col_po = 0;

    if (std::abs(gravity_x) + std::abs(gravity_y) < 2 * hardballs::f_type_large_epsilon) {
        table.SetStraight();
    } else {
        table.SetGravity(gravity_x, gravity_y);
    }
    
    table.AddObstacleHLine(bottom_wall, col_rule_bottom, alpha_bottom, beta_bottom);
    table.AddObstacleHLine(top_wall, col_rule_top, alpha_top, beta_top);
   
    table.SetParticlesNumber(N);
    table.SetRadars((std::sqrt(N) * 15) / 2);

    
    {
        // infected status
        infecta.resize(N);
        infecta.assign(N,0);
        infecta[0] = 1;

        // normal velocities with average energy 1
        std::vector<f_type> vvx(N, 0.);
        std::vector<f_type> vvy(N, 0.);
        f_type total_energy2 = 0.;
        for (int k=0; k<N; k++) {
            vvx[k] = f_rand_norm(0,1);
            vvy[k] = f_rand_norm(0,1);
            total_energy2 += (vvx[k]*vvx[k] + vvy[k]*vvy[k]);
        }
        f_type energy_correction = std::sqrt(N / total_energy2);


        f_type cr = 1.01 * r;
        f_type real_r = r;
        int max_N = 0;

        for (int test = 0; test < 2; test++) {
            f_type x = 0., y = bottom_wall + cr;
            int k = 0, layer = 0;
            while((test == 0 || (test == 1 && k < N)) && y < top_wall - cr) {
                if (
#ifdef heavy0
                std::sqrt(SQR(x-test_c_x) + SQR(y-test_c_y)) > test_c_r + cr
#else
                true
#endif
                ) {
                    if (test == 0) {
                        max_N ++;
                    } else {
                        table.SetParticlePosition(k, x, y);
                        f_type vx = vvx[k] * energy_correction;
                        f_type vy = vvy[k] * energy_correction;
                        table.SetParticleVelocity(k, vx, vy);
                        table.SetParticleRadius(k, real_r);
                        table.SetParticleMass(k, 1);
                        k ++;
                    }
                }
                
                if (test == 0) {
                    x += 2 * cr;
                } else {
                    x += (max_N / N) * 2 * cr;
                }
               
                if (((layer % 2 == 1) && x > 1. - cr) || ((layer % 2 == 0) && x > 1. - 2. * cr)) {
                    layer ++;
                    x = (layer % 2) ? cr : 0.;
                    y += std::sqrt(3.) * cr;
                }
            }
            if (test == 0 && max_N < N) {
                real_r = 0.0;
                cr = 0.5 / N;
                max_N = N;
            }
        }
#ifdef heavy0
        { 
            f_type m0 = (real_r > 0) ? SQR(test_c_r / real_r) : 1.0;
            f_type vx0 = table.vx(0) / std::sqrt(m0);
            f_type vy0 = table.vy(0) / std::sqrt(m0);
            table.SetParticleRadius(0, test_c_r);
            table.SetParticleMass(0, m0);
            table.SetParticlePosition(0, test_c_x, test_c_y);
            table.SetParticleVelocity(0, vx0, vy0);
        }
#endif
    }

    table.SetReady();
}

void table_adjust(char col_rule_bottom, char col_rule_top, double gravity_x, double gravity_y, double alpha_bottom, double beta_bottom, double alpha_top, double beta_top) {
    if (std::abs(gravity_x) + std::abs(gravity_y) < 2 * hardballs::f_type_large_epsilon) {
        table.SetStraight();
    } else {
        table.SetGravity(gravity_x, gravity_y);
    }
    table.RemoveObstacles();
    table.AddObstacleHLine(bottom_wall, col_rule_bottom, alpha_bottom, beta_bottom);
    table.AddObstacleHLine(top_wall, col_rule_top, alpha_top, beta_top);
    table.SetReady();
}

double table_time() {
    return double(table.T().f_val());
}

int table_N() {
    return table.N();
}

particle table_particle(int n) {
    particle p;
    p.x = table.x(n);
    p.y = table.y(n);
    p.vx = table.vx(n);
    p.vy = table.vy(n);
    p.r = table.r(n);
    p.m = table.m(n);
    p.infected = infecta[n];
    return p;
}

void table_live(double t) {
    // live until time t, record all collisions with obstacles
    hist_po_.resize(0);
    hardballs::event_t ev;
    do {
        table.Live(long_time_t(t), hardballs::EV_PO_COLL | hardballs::EV_PP_COLL);
        ev = table.CurrentEvent();
        if (ev.kind == hardballs::EV_PO_COLL) {
            num_col_po ++;
            hist_po_.push_back(table_particle(ev.particle));
        } else if (ev.kind == hardballs::EV_PP_COLL) {
            num_col_pp ++;
            size_t p1 = ev.particle,
                   p2 = ev.object;
            if (infecta[p1] == 1 || infecta[p2] == 1) {
                infecta[p1] = infecta[p2] = 1;
            }
        }
    } while (ev.kind != hardballs::EV_STOP);
}

EMSCRIPTEN_BINDINGS(my_module) {
    value_object<particle>("particle")
        .field("x", &particle::x)
        .field("y", &particle::y)
        .field("vx", &particle::vx)
        .field("vy", &particle::vy)
        .field("r", &particle::r)
        .field("m", &particle::m)
        .field("infected", &particle::infected);
    function("table_init", &table_init);
    function("table_adjust", &table_adjust);
    function("table_live", &table_live);
    function("table_time", &table_time);
    function("table_particle", &table_particle);
    function("table_N", &table_N);
    function("table_num_col_pp", &table_num_col_pp);
    function("table_num_col_po", &table_num_col_po);
    function("table_num_col", &table_num_col);
    function("table_log", &table_log);
    function("hist_po_N", &hist_po_N);
    function("hist_po", &hist_po);
}

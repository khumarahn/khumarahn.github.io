/*
   Standard table: compute an unfolded trajectory
*/

#include <time.h>
#include <emscripten/bind.h>
#include <vector>
#include "hardballs.h"
using namespace std;
using namespace hardballs;
using namespace emscripten;

const size_t max_orbit_size = 32*1024;

struct xy_t {
    double x;
    double y;
    xy_t(f_type xx, f_type yy) : x(xx), y(yy) {};
    xy_t() {};
};

class orbit_c {
    private:
        table_c<TORUS> table_;
        vector<xy_t> xy_;
        long long step_size_;
        long long step_count_;
    public:
        orbit_c() {
            srand (time(NULL));
            
            table_.SetStraight();

            table_.AddObstacleCircle(0.0, 0.0, 0.4);
            table_.AddObstacleCircle(0.5, 0.5, 0.2);

            table_.SetParticlesNumber(1);

            table_.SetParticlePosition(0, 0.3, 0.3);
            // velocities masses radii
            {
                f_type vvx = f_rand_norm(1, 0.4),
                       vvy = 0.99,
                       vvv = sqrt(vvx*vvx + vvy*vvy);
                table_.SetParticleVelocity(0, vvx / vvv, vvy / vvv);
                table_.SetParticleMass(0, 1.);
                table_.SetParticleRadius(0, 0);
            }

            table_.SetReady();

            step_size_ = 1;
            step_count_ = 0;

            xy_.reserve(2 * max_orbit_size);
            
            xy_.push_back(xy_t(table_.X(0), table_.Y(0)));
        };
        void step(double dt) {
            event_kind_t stop_events = EV_PO_COLL;
            for (;;) {
                table_.Live(table_.T() + dt, stop_events);
                if (table_.CurrentEvent().kind == stop_events) {
                    step_count_++;
                    if (step_count_ == step_size_) {
                        xy_.push_back(xy_t(table_.X(0), table_.Y(0)));
                        step_count_ = 0;
                        if (xy_.size() > max_orbit_size) {
                            for (size_t i=1; i*2 < xy_.size(); i++) {
                                xy_[i] = xy_[i*2];
                            }
                            xy_.resize((xy_.size() + 1)  / 2);
                            step_size_ *= 2;
                        }
                    }
                } else {
                    break;
                }
            }

        }
        xy_t xy(int n) {
            if (n < xy_.size()) {
                return xy_[n]; 
            } else {
                return xy_t(table_.X(0), table_.Y(0));
            }
        };
        long size() { return xy_.size(); };
        double T() { return double(table_.T().f_val()); };
};

// remove sASSERTIONS

EMSCRIPTEN_BINDINGS(my_module) {
    class_<orbit_c>("orbit_c")
    .constructor()
    .function("step", &orbit_c::step)
    .function("size", &orbit_c::size)
    .function("T", &orbit_c::T)
    .function("xy", &orbit_c::xy);
    value_object<xy_t>("xy_t")
        .field("x",     &xy_t::x)
        .field("y",     &xy_t::y);
    //register_vector<xy_t>("vector<xy_t>");
}

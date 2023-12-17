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
            orbit_c(1);
        }
        orbit_c(int t) {
            srand(time(NULL));

            table_.SetStraight();
            table_.SetParticlesNumber(1);
            if (t==0) { // standard layout
                // obstacles
                table_.AddObstacleCircle(0.0, 0.0, 0.4);
                table_.AddObstacleCircle(0.5, 0.5, 0.2);


                table_.SetParticlePosition(0, 0.3, 0.3);
            } else { // skewed layout
                // obstacles : long arcs
                for (int ii=-1; ii<=1; ii+=2) {
                    f_type d = 16.*ii,
                           x1 = +1.25,
                           x2 = -1.25,
                           y1 = 3./16.,
                           y2 = 13./16.,
                           xm = 0.5*(x1+x2),
                           ym = 0.5*(y1+y2),
                           n = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)),
                           xc = xm+(y1-y2)/n*d,
                           yc = ym-(x1-x2)/n*d,
                           R = sqrt(d*d+n*n/4.),
                           theta = (ii==1) ? (atan2(y1-yc, x1-xc)): (atan2(y2-yc, x2-xc)),
                           omega = 2*asin(n/2./R);
                    for (int i=0; i<4; i++) {
                        f_type x1x = xc + R*cos(theta+omega*i/4.),
                               y1y = yc + R*sin(theta+omega*i/4.);
                        table_.AddObstacleCircleArc(xc-floor(x1x),yc-floor(y1y),R,theta+omega*i/4.,omega/4.);
                    }
                }

                // obstacles: short arcs
                {
                    f_type d = 1,
                           R = sqrt(0.1875*0.1875+d*d),
                           theta = asin(0.1875/R),
                           omega = 2*theta;
                    table_.AddObstacleCircleArc(0.5+d,0.,R,PI-theta,omega);
                    table_.AddObstacleCircleArc(0.5-d,1.,R,-theta,omega);
                }
                // particle
                table_.SetParticlePosition(0, 0.5, 0.5);
            } // skewed layout

            table_.SetParticleRadius(0, 0);        
            table_.SetParticleMass(0, 1.0);
            
            // velocities masses radii
            {
                f_type vvx = f_rand_norm(1, 0.4),
                       vvy = 0.99,
                       vvv = sqrt(vvx*vvx + vvy*vvy);
                table_.SetParticleVelocity(0, vvx / vvv, vvy / vvv);
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
    .constructor<int>()
    .function("step", &orbit_c::step)
    .function("size", &orbit_c::size)
    .function("T", &orbit_c::T)
    .function("xy", &orbit_c::xy);
    value_object<xy_t>("xy_t")
        .field("x",     &xy_t::x)
        .field("y",     &xy_t::y);
    //register_vector<xy_t>("vector<xy_t>");
}

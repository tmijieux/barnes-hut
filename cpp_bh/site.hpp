#ifndef TDP_SITE_H
#define TDP_SITE_H

#include <cstdbool>
#include <cstdint>
#include <memory>
#include "util.hpp"

#include "particle.hpp"

namespace barnes_hut {

class site;

typedef std::shared_ptr<site> site_ptr;

class site {
    bool _is_leaf, _is_empty;
    particle *_particle;

    double _min_x, _min_y, _max_x, _max_y;
    double _x_width, _y_height;
    particle _mass_center;

    site_ptr _up_left;
    site_ptr _up_right;
    site_ptr _down_left;
    site_ptr _down_right;
    void subdivide_site();

public:
    site(double min_x, double min_y, double max_x, double max_y);
    void insert_particle(particle *p);
    particle_vec tree_gen(int64_t particles_count);
    void compute_bh_force(particle *p);
    void update_center(particle *p);

    static site_ptr make_shared(double min_x, double min_y, double max_x, double max_y);
};

#define THRESHOLD (0.707)
#define GRAVITATIONAL_CONSTANT (6.67408e-11)
#define COMPUTE_FORCE(SRC_, DST_)                                       \
do {                                                                    \
    const double _a = (SRC_)->x - (DST_)->x;                            \
    const double _b = (SRC_)->y - (DST_)->y;                            \
    const double u_n2 = _a*_a + _b*_b;                                  \
    const double f_n = ((GRAVITATIONAL_CONSTANT                         \
                         * (DST_)->mass * (SRC_)->mass) / u_n2);        \
    const double u_n = sqrt(u_n2);                                      \
    (DST_)->fx += f_n * (_a / u_n);                                     \
    (DST_)->fy += f_n * (_b / u_n);                                     \
} while (0)

}; // end namespace barnes_hut

#endif // TDP_SITE_H

#ifndef TDP_PARTICLE_H
#define TDP_PARTICLE_H

#include "util.hpp"

namespace barnes_hut {

struct particle;
typedef std::vector<particle>  particle_vec;

struct particle {
    double x, y;
    double mass;
    double fx, fy;
    particle();
    particle(const particle&);
    particle(double min_x, double min_y, double max_x, double max_y);
    bool operator<(const particle &);
};

}; // end namespace barnes_hut

#endif // TDP_PARTICLE_H

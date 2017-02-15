#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "particle.hpp"

using namespace barnes_hut;

static inline double
drand_(double min, double max)
{
    return ((double) rand() / RAND_MAX) * (max - min) + min;
}

static inline int64_t interleave(uint64_t a, uint64_t b)
{
    static const uint64_t B[] = {
        0x5555555555555555,
        0x3333333333333333,
        0x0F0F0F0F0F0F0F0F,
        0x00FF00FF00FF00FF,
        0x0000FFFF0000FFFF
    };
    static const unsigned int S[] = {1, 2, 4, 8, 16};

    a = (a | (a << S[4])) & B[4];
    a = (a | (a << S[3])) & B[3];
    a = (a | (a << S[2])) & B[2];
    a = (a | (a << S[1])) & B[1];
    a = (a | (a << S[0])) & B[0];

    b = (b | (b << S[4])) & B[4];
    b = (b | (b << S[3])) & B[3];
    b = (b | (b << S[2])) & B[2];
    b = (b | (b << S[1])) & B[1];
    b = (b | (b << S[0])) & B[0];

    return  a | (b << 1);
}

particle::particle(double min_x, double min_y,
                   double max_x, double max_y, int64_t nb_particles):
    x(drand_(min_x, max_x)),
    y(drand_(min_y, max_y)),
    mass(drand_(2000.0, 300000.0)),
    fx(0.0), fy(0.0),
    z(0)
{
    uint64_t xx = ((x-min_x) / (max_x-min_x)) * 2*nb_particles;
    uint64_t yy = ((y-min_y) / (max_y-min_y)) * 2*nb_particles;
    z = interleave(xx, yy);
}

particle::particle():
    x(0.0), y(0.0),  mass(0.0), fx(0.0), fy(0.0), z(0)
{
}

particle::particle(const particle& o):
    x(o.x), y(o.y), mass(o.mass), fx(0.0), fy(0.0), z(o.z)
{
}

bool particle::operator<(const particle &o)
{
    return z < o.z;
}

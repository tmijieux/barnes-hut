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

particle::particle(double min_x, double min_y,
                   double max_x, double max_y):
    x(drand_(min_x, max_x)),
    y(drand_(min_y, max_y)),
    mass(drand_(2000.0, 300000.0)),
    fx(0.0), fy(0.0)
{
}

particle::particle():
    x(0.0), y(0.0),  mass(0.0), fx(0.0), fy(0.0)
{
}

particle::particle(const particle& o):
    x(o.x), y(o.y), mass(o.mass), fx(0.0), fy(0.0)
{
}

bool particle::operator<(const particle &o)
{
    if (x < o.x)
        return true;
    if (x > o.x)
        return false;
    if (y < o.y)
        return true;
    return false;
}

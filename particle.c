#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "particle.h"

static inline double
drand_(double min, double max)
{
    return ((double) rand() / RAND_MAX) * (max - min) + min;
}

void tdp_particle_init_random(tdp_particle *p,
                              double min_x, double min_y,
			      double max_x, double max_y)
{
    memset(p, 0, sizeof*p);
    p->x = drand_(min_x, max_x);
    p->y = drand_(min_y, max_y);
    p->mass = drand_(2000.0, 300000.0);
    p->fx = 0.0;
    p->fy = 0.0;
}

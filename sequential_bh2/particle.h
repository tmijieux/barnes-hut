#ifndef TDP_PARTICLE_H
#define TDP_PARTICLE_H

typedef struct tdp_particle_ tdp_particle;

struct tdp_particle_ {
  double mass;
  double x, y;
  double fx, fy;
};

void tdp_particle_init_random(
    tdp_particle *p, double min_x, double min_y, double max_x, double max_y);

#endif // TDP_PARTICLE_H

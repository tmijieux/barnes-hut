#ifndef TDP_SITE_H
#define TDP_SITE_H

#include <stdbool.h>
#include <stdint.h>
#include "util.h"

#include "particle.h"

typedef struct tdp_site_ tdp_site;

struct tdp_site_ {
    bool is_leaf, is_empty;

    tdp_particle *particle;
    tdp_particle mass_center;

    double min_x, min_y, max_x, max_y;
    double x_width, y_height;

    tdp_site *up_left;
    tdp_site *up_right;
    tdp_site *down_left;
    tdp_site *down_right;
};

void tdp_tree_alloc(uint64_t N);
void tdp_site_init(tdp_site *site,
                   double min_x, double min_y, double max_x, double max_y);
tdp_site *tdp_site_new(double min_x, double min_y, double max_x, double max_y);
void tdp_site_insert_particle(tdp_site *tree, tdp_particle *p);
tdp_particle *tdp_site_tree_gen(tdp_site *tree, int64_t particles_count);

void tdp_site_compute_bh_force(tdp_site *site, tdp_particle *p);
void tdp_site_update_center(tdp_site *tree, tdp_particle *p);
void tdp_leaves_copy(tdp_particle *dest, tdp_particle *src);
void tdp_site_free(tdp_site *tree);

#define THRESHOLD (0.707)
//#define THRESHOLD (0.05)

#define GRAVITATIONAL_CONSTANT (6.67408e-11)
#define COMPUTE_FORCE(SRC_, DST_)					\
    do {								\
	double f_n, u_x, u_y, u_n;					\
	u_n = hypot((DST_)->x - (SRC_)->x, (DST_)->y - (SRC_)->y);	\
	f_n = ((GRAVITATIONAL_CONSTANT					\
		* (DST_)->mass * (SRC_)->mass) / (u_n*u_n));		\
	u_x = ((SRC_)->x - (DST_)->x) / u_n;				\
	u_y = ((SRC_)->y - (DST_)->y) / u_n;				\
	(DST_)->fx += f_n * u_x;					\
	(DST_)->fy += f_n * u_y;					\
    } while (0)

#endif // TDP_SITE_H

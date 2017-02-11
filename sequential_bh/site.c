#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "site.h"
#include "util.h"

void tdp_site_init(tdp_site *site,
                   double min_x, double min_y,
                   double max_x, double max_y)
{
    memset(site, 0, sizeof*site);

    site->is_empty = true;
    site->is_leaf = true;
    site->particle = NULL;

    site->min_x = min(min_x, max_x);
    site->min_y = min(min_y, max_y);
    site->max_x = max(min_x, max_x);
    site->max_y = max(min_y, max_y);

    site->x_width = site->max_x - site->min_x;
    site->y_height = site->max_y - site->min_y;
}

tdp_site *tdp_site_new(double min_x, double min_y,
		       double max_x, double max_y)
{
    tdp_site *site = malloc(sizeof*site);
    tdp_site_init(site, min_x, min_y, max_x, max_y);
    return site;
}

static void subdivide_site(tdp_site *site)
{
    site->up_left = tdp_site_new(site->min_x,
                                 site->min_y + (site->y_height/2.0),
                                 site->min_x + (site->x_width/2.0),
                                 site->max_y);

    site->up_right = tdp_site_new(site->min_x + (site->x_width/2.0),
                                  site->min_y + (site->y_height/2.0),
                                  site->max_x,
                                  site->max_y);

    site->down_left = tdp_site_new(site->min_x,
                                   site->min_y,
                                   site->min_x + (site->x_width/2.0),
                                   site->min_y + (site->y_height/2.0));

    site->down_right = tdp_site_new(site->min_x + (site->x_width/2.0),
                                    site->min_y,
                                    site->max_x,
                                    site->min_y + (site->y_height/2.0));
}

void tdp_site_insert_particle(tdp_site *tree, tdp_particle *p)
{
    if (tree->is_leaf){
        if (tree->is_empty) { // does not contain a body
            tree->is_empty = false;
            tree->particle = p;
        } else {
            // two bodies in the same region
            // Subdividing the region further by creating four children

            tdp_particle *b = tree->particle;
            assert( b != NULL );

            tree->is_leaf = false;
            tree->is_empty = false;
            tree->particle = NULL;

            subdivide_site(tree);

            tdp_site_insert_particle(tree, p);
            tdp_site_insert_particle(tree, b);
        }
    } else { // intern (not leaf) node
        tdp_site_update_center(tree, p);

        //recursively inserting body in appropriate quadrant
        if (p->x > tree->min_x + (tree->x_width/2.0))
        {
	    if (p->y > tree->min_y + (tree->y_height/2.0))
                tdp_site_insert_particle(tree->up_right, p);
	    else
                tdp_site_insert_particle(tree->down_right, p);
        }
        else
        {
	    if (p->y > tree->min_y + (tree->y_height/2.0))
                tdp_site_insert_particle(tree->up_left, p);
	    else
                tdp_site_insert_particle(tree->down_left, p);
        }
    }
}

void tdp_site_update_center(tdp_site *tree, tdp_particle *p)
{
    tree->mass_center.x = tree->mass_center.x*tree->mass_center.mass + p->mass*p->x;
    tree->mass_center.y = tree->mass_center.y*tree->mass_center.mass + p->mass*p->y;

    tree->mass_center.mass += p->mass;

    tree->mass_center.x /= tree->mass_center.mass;
    tree->mass_center.y /= tree->mass_center.mass;
}

tdp_particle *tdp_site_tree_gen(tdp_site *site, int64_t particles_count)
{
    tdp_particle *buf = malloc(sizeof*buf * particles_count);
    for (int64_t i = 0; i < particles_count; ++i) {
        tdp_particle_init_random(
            buf+i, site->min_x, site->min_y, site->max_x, site->max_y);
        tdp_site_insert_particle(site, buf+i);
    }
    return buf;
}

void tdp_site_compute_bh_force(tdp_site *site, tdp_particle *p)
{
    if (site->is_leaf) {
        if (site->is_empty || site->particle == p)
            return;
        COMPUTE_FORCE(site->particle, p);
    } else {
        double distance, width;
        width = max(site->x_width, site->y_height);
        distance = hypot(p->x-site->mass_center.x, p->y-site->mass_center.y);
        
        if ((width / distance) < THRESHOLD)   // Under threshold: approximation
            COMPUTE_FORCE(&site->mass_center, p);
        else {
            tdp_site_compute_bh_force(site->up_left, p);
            tdp_site_compute_bh_force(site->up_right, p);
            tdp_site_compute_bh_force(site->down_left, p);
            tdp_site_compute_bh_force(site->down_right, p);
        }
    }
}

void tdp_leaves_copy(tdp_particle *dest, tdp_particle *src)
{
    dest->mass = src->mass;
    dest->x = src->x;
    dest->y = src->y;
    dest->fx = 0.0;
    dest->fy = 0.0;
}

void tdp_site_free(tdp_site *tree)
{
    if (!tree->is_leaf) {
        tdp_site_free(tree->up_right);
        tdp_site_free(tree->up_left);
        tdp_site_free(tree->down_right);
        tdp_site_free(tree->down_left);
    }
    free(tree);
}

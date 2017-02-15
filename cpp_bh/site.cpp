#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>

#include "util.hpp"
#include "site.hpp"
#include "util.hpp"

using namespace barnes_hut;

site::site( double min_x, double min_y, double max_x, double max_y):
    _is_leaf(true), _is_empty(true),
    _particle(NULL),
    _min_x(std::min(min_x, max_x)),
    _min_y(std::min(min_y, max_y)),
    _max_x(std::max(min_x, max_x)),
    _max_y(std::max(min_y, max_y)),
    _x_width(_max_x-_min_x),
    _y_height(_max_y-_min_y)
{
}

site_ptr site::make_shared(double min_x, double min_y, double max_x, double max_y)
{
    return std::make_shared<site>(min_x, min_y, max_x, max_y);
}

void site::subdivide_site()
{
    _up_left = make_shared(_min_x,
                           _min_y + (_y_height/2.0),
                           _min_x + (_x_width/2.0),
                           _max_y);

    _up_right = make_shared(_min_x + (_x_width/2.0),
                            _min_y + (_y_height/2.0),
                            _max_x,
                            _max_y);

    _down_left = make_shared(_min_x,
                             _min_y,
                             _min_x + (_x_width/2.0),
                             _min_y + (_y_height/2.0));

    _down_right = make_shared(_min_x + (_x_width/2.0),
                              _min_y,
                              _max_x,
                              _min_y + (_y_height/2.0));
}

void site::insert_particle(particle *p)
{
    if (_is_leaf){
        if (_is_empty) { // does not contain a body
            _is_empty = false;
            _particle = p;
        } else {
            // two bodies in the same region
            // Subdividing the region further by creating four children

            particle *b = _particle;
            assert( b != NULL );

            _is_leaf = false;
            _is_empty = false;
            _particle = NULL;

            subdivide_site();

            insert_particle(p);
            insert_particle(b);
        }
    } else { // intern (not leaf) node
        update_center(p);

        //recursively inserting body in appropriate quadrant
        if (p->x > _min_x + (_x_width/2.0)) {
	    if (p->y > _min_y + (_y_height/2.0))
                _up_right->insert_particle(p);
	    else
                _down_right->insert_particle(p);
        } else {
	    if (p->y > _min_y + (_y_height/2.0))
                _up_left->insert_particle(p);
	    else
                _down_left->insert_particle(p);
        }
    }
}

void site::update_center(particle *p)
{
    _mass_center.x = _mass_center.x*_mass_center.mass + p->mass*p->x;
    _mass_center.y = _mass_center.y*_mass_center.mass + p->mass*p->y;

    _mass_center.mass += p->mass;

    _mass_center.x /= _mass_center.mass;
    _mass_center.y /= _mass_center.mass;
}

particle_vec site::tree_gen(int64_t particles_count)
{
    particle_vec v;
    v.reserve(particles_count);

    for (int64_t i = 0; i < particles_count; ++i)
        v.emplace_back( _min_x, _min_y, _max_x, _max_y, particles_count);

    std::sort(v.begin(), v.end());

    for (auto &p : v)
        insert_particle(&p);

    return v;
}

void site::compute_bh_force(particle *p)
{
    if (_is_leaf) {
        if (_is_empty || _particle == p)
            return;
        COMPUTE_FORCE(_particle, p);
    } else {
        double distance, width;
        width = std::max(_x_width, _y_height);
        distance = hypot(p->x-_mass_center.x, p->y-_mass_center.y);

        if ((width / distance) < THRESHOLD)   // Under threshold: approximation
            COMPUTE_FORCE(&_mass_center, p);
        else {
            _up_left->compute_bh_force(p);
            _up_right->compute_bh_force(p);
            _down_left->compute_bh_force(p);
            _down_right->compute_bh_force(p);
        }
    }
}

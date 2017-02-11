#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "../perf/perf.h"
#include "site.h"
#include "util.h"

#define TIME( instr, msg ) ({                           \
            perf_t p1_, p2_;                            \
            perf(&p1_);                                 \
            { instr; }                                  \
            perf(&p2_);                                 \
            perf_diff(&p1_, &p2_);                      \
            double T_ = perf_get_seconds(&p2_);         \
            fprintf(stderr, "time " msg " %g s\n", T_); \
            T_;                                         \
        })

static int compar_double(const void *a, const void *b)
{
    if (*(double*)a < *(double*)b)
        return -1;
    if (*(double*)a > *(double*)b)
        return 1;
    return 0;
}

static void compute_difference(
    int64_t particle_count, tdp_particle *particles, tdp_particle *particles_copy)
{
    int64_t p_count = 0, nan_count = 0, inf_count = 0;
    double total_err = 0.0;
    double *errors = malloc(sizeof*errors * particle_count);
    double err_max = 0.0;

    for (int64_t i = 0; i < particle_count; ++i) {
        tdp_particle *p_bh = particles+i;
        tdp_particle *p_exact = particles_copy+i;

        double f_diff = hypot(p_exact->fx-p_bh->fx, p_exact->fy-p_bh->fy);
        double f_exact = hypot(p_exact->fx, p_exact->fy);

        double err = f_diff / f_exact;
        nan_count += !!isnan(err);
        inf_count += !!isinf(err);

        if (isnan(err) || isinf(err))
            continue;
        errors[p_count] = err;
        total_err += err;
        err_max = max(err, err_max);
        ++p_count;
    }

    total_err /= p_count;

    qsort(errors, p_count, sizeof errors[0], compar_double);

    printf("\n\n\nnombre total de particules: %ld\n", particle_count);
    printf("erreur maximale: %g%%\n", err_max*100);
    printf("erreur moyenne: %g%%\n", total_err*100);
    printf("1er quartile: %g%%\n", errors[(int)((double)p_count/4.)]*100);
    printf("erreur médiane: %g%%\n", errors[(int)((double)p_count/2.)]*100);
    printf("3ème quartile: %g%%\n", errors[(int)(3.*p_count/4.)]*100);

    printf("%ld erreur de calcul\n", particle_count - p_count);
    printf("%ld nan values\n", nan_count);
    printf("%ld inf values\n\n\n", inf_count);

    free(errors);
}

static void compute_barnes_hut_force(
    tdp_site *site, tdp_particle *tree_leaves, int64_t nb_particles)
{
    printf("COMPUTE barnes-hut forces >>\n");

    #pragma omp parallel for
    for (int64_t i = 0; i < nb_particles; ++i)
        tdp_site_compute_bh_force(site, tree_leaves+i);
    puts("<<END");
}

/*
 * ""exact"" n^2 method
 */
static void compute_exact_force(tdp_particle *particles, int64_t particle_count)
{
    puts("COMPUTE exact forces>>");

    #pragma omp parallel for collapse(2)
    for (int64_t i = 0; i < particle_count; ++i) {
        for (int64_t j = 0; j < particle_count; ++j) {
            if (j != i)
                COMPUTE_FORCE(particles+j, particles+i);
        }
    }
    puts("<<END");
}

static void barnes_hut(const int64_t N)
{
    tdp_tree_alloc(N);
    
    tdp_site *tree = tdp_site_new(0.0, 0.0, 200.0, 200.0);
    tdp_particle *tree_leaves;

    TIME( tree_leaves = tdp_site_tree_gen(tree, N), "tree_gen" ) ; // génération de l'arbre

    tdp_particle *leaves_copy = malloc(sizeof(tdp_particle) * N);
    for (int64_t i = 0; i<N; ++i)
        tdp_leaves_copy(leaves_copy+i, tree_leaves+i);

    double t1, t2;
    t1 = TIME( compute_barnes_hut_force(tree, tree_leaves, N), "Barnes-Hut Forces" );
    t2 = TIME( compute_exact_force(leaves_copy, N), "exact Forces" );
    printf("Barnes-Hut SPEEDUP: %g\n", t2/t1);

    compute_difference(N, tree_leaves, leaves_copy);

    free(tree_leaves);
    free(leaves_copy);
    tdp_site_free(tree);
}

int main(int argc, char *argv[])
{
    srand(time(NULL)+(long)&argc);
    int64_t N = 20000;

    if (argc >= 2)
        N = strtol(argv[1], NULL, 10);
    ASSERT_MSG( N >= 2,
                "interactions en PLUSIEURS corps, N doit être supérieur ou égal à 2");

    barnes_hut(N);
    return EXIT_SUCCESS;
}

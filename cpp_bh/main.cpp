#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <ctime>

#include "../perf/perf.h"
#include "site.hpp"
#include "util.hpp"
#include "proc.hpp"

#define TIME( instr, msg ) ({                                           \
            perf_t p1_, p2_;                                            \
            perf(&p1_);                                                 \
            { instr; }                                                  \
            perf(&p2_);                                                 \
            perf_diff(&p1_, &p2_);                                      \
            double T_ = perf_get_seconds(&p2_);                         \
            std::cout << "time " msg " "<< T_<<" s"<< std::endl;        \
            T_;                                                         \
        })

using namespace barnes_hut;

static void compute_difference(
    particle_vec &particles_bh, particle_vec &particles_exact)
{
    const uint64_t N = particles_bh.size();
    uint64_t nan_count = 0, inf_count = 0;
    double total_err = 0.0, total_N_err = 0.0, total_N = 0.0, err_max = 0.0;
    std::vector<double> errors;

    errors.reserve(N);

    for (uint64_t i = 0; i < N; ++i) {
        particle *p_bh = &particles_bh[i], *p_exact = &particles_exact[i];
        double f_diff = hypot(p_exact->fx-p_bh->fx, p_exact->fy-p_bh->fy);
        double f_exact = hypot(p_exact->fx, p_exact->fy);

        double err = f_diff / f_exact;
        nan_count += !!std::isnan(err);
        inf_count += !!std::isinf(err);

        if (std::isnan(err) || std::isinf(err))
            continue;
        errors.push_back(err);
        total_err += err;
        total_N_err += err * f_exact;
        total_N += f_exact;
        err_max = std::max(err, err_max);
    }

    int p_count = errors.size();
    total_err /= p_count;
    total_N_err /= total_N;

    std::sort(errors.begin(), errors.end());

    using std::endl;
    std::cout
        << endl << endl
        << "nombre total de particules: " << N << endl
        << "erreur maximale: " << err_max*100 << "%" << endl
        << "erreur moyenne: " << total_err*100 << "%" << endl
        << "erreur moyenne pondéré: " << total_N_err*100 <<"%"<< endl

        << "1er quartile: " << errors[(int)((double)p_count/4.)]*100 << "%" << endl
        << "erreur médiane: " << errors[(int)((double)p_count/2.)]*100 <<"%" << endl
        << "3ème quartile: " << errors[(int)(3.*p_count/4.)]*100 << "%" << endl

        << N - p_count << " erreur de calcul" << endl
        << nan_count << " nan values" << endl
        << inf_count << " inf values" << endl << endl << endl;
}

static void
compute_barnes_hut_force(site_ptr site, particle_vec &v)
{
    std::cout<< "COMPUTE barnes-hut forces >>" << std::endl;
    const uint64_t N = v.size();

    #pragma omp parallel for
    for (uint64_t i = 0; i < N; ++i)
        site->compute_bh_force(&v[i]);
    std::cout<< "<<END" << std::endl;
}

/*
 * ""exact"" n^2 method
 */
static void compute_exact_force(particle_vec &v)
{
    puts("COMPUTE exact forces>>");
    const uint64_t N = v.size();

    #pragma omp parallel for
    for (uint64_t i = 0; i < N; ++i) {
        for (uint64_t j = 0; j < N; ++j) {
            if (j != i)
                COMPUTE_FORCE(&v[j], &v[i]);
        }
    }
    puts("<<END");
}

static void do_barnes_hut(const int64_t N)
{
    site_ptr tree = site::make_shared(0.0, 0.0, 200.0, 200.0);
    particle_vec leaves_bh, leaves_exact;

    TIME( leaves_bh = tree->tree_gen(N), "tree_gen" ) ; // génération de l'arbre

    leaves_exact = leaves_bh;

    double t1, t2;
    t1 = TIME( compute_barnes_hut_force(tree, leaves_bh), "Barnes-Hut Forces" );
    t2 = TIME( compute_exact_force(leaves_exact), "exact Forces" );
    std::cout << "Barnes-Hut SPEEDUP: " << t2/t1 << std::endl;

    compute_difference(leaves_bh, leaves_exact);
}

int main(int argc, char *argv[])
{
    srand(time(NULL)+(long)&argc);
    int64_t N = 20000;

    if (argc >= 2)
        N = strtol(argv[1], NULL, 10);
    ASSERT_MSG( N >= 2,
                "interactions en PLUSIEURS corps, N doit être supérieur ou égal à 2");

    do_barnes_hut(N);

    return EXIT_SUCCESS;
}

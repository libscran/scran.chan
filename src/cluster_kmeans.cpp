#include "config.h"

#include "Rcpp.h"
#include "knncolle.h"
#include "cluster_kmeans.h"
#ifdef _OPENMP
#include "omp.h"
#endif

//[[Rcpp::export(rng=false)]]
SEXP cluster_kmeans(Rcpp::NumericMatrix data, int nclusters, int init_method, int nthreads) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    Kmeans runner(data, nclusters, init_method);
    runner.run();
    return runner.yield();
}

#include "config.h"

#include "Rcpp.h"
#include "knncolle.h"
#include "cluster_snn_graph.h"
#ifdef _OPENMP
#include "omp.h"
#endif

//[[Rcpp::export(rng=false)]]
SEXP build_graph(SEXP nnptr, int k, std::string method, double resolution, int nthreads) {
    KnncollePtr nns(nnptr);
    scran::BuildSNNGraph builder;
    builder.set_neighbors(k);

#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    return Rcpp::XPtr<SNNGraph>(new SNNGraph(nns->nobs(), builder.run(nns.get()), method, resolution));
}

//[[Rcpp::export(rng=false)]]
SEXP cluster_graph(SEXP ptr) {
    Rcpp::XPtr<SNNGraph> x(ptr);
    x->run();
    return x->yield();
}

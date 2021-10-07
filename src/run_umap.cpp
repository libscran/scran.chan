#include "Rcpp.h"
#include "knncolle.h"
#include "run_umap.h"
#ifdef _OPENMP
#include "omp.h"
#endif

//[[Rcpp::export(rng=false)]]
SEXP initialize_umap(SEXP nnptr, int num_neighbors, int nthreads) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    KnncollePtr nns(nnptr);

    umap runner;
    runner.set_num_neighbors(num_neighbors);
    std::vector<float> embedding(2 * nns->nobs());
    auto status = runner.initialize(nns.get(), 2, embedding.data());

    return Rcpp::XPtr<InitializedUmap>(new InitializedUmap(std::move(runner), std::move(status), std::move(embedding), nns->nobs()));
}

//[[Rcpp::export(rng=false)]]
SEXP run_umap(SEXP init) {
    Rcpp::XPtr<InitializedUmap> ptr(init);
    ptr->run();
    return ptr->yield();
}


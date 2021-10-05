#include "Rcpp.h"
#include "knncolle.h"
#include "umappp/umappp.hpp"
#ifdef _OPENMP
#include "omp.h"
#endif

typedef umappp::Umap<float> umap;

struct InitializedUmap {
    InitializedUmap(umap r, umap::Status s, std::vector<float> e, int n) : runner(std::move(r)), status(std::move(s)), embedding(std::move(e)), nobs(n) {}
    umap runner;
    umap::Status status;
    std::vector<float> embedding;
    int nobs;
};

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
    auto t = ptr->runner;
    t.run(ptr->status, 2, ptr->embedding.data());

    Rcpp::NumericMatrix output(2, ptr->nobs);
    std::copy(ptr->embedding.begin(), ptr->embedding.end(), output.begin());
    return output;
}


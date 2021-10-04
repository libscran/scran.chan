#include "Rcpp.h"
#include "knncolle.h"
#include "qdtsne/qdtsne.hpp"
#ifdef _OPENMP
#include "omp.h"
#endif

typedef qdtsne::Tsne<2, float> tsne;

struct InitializedTsne {
    InitializedTsne(tsne r, tsne::Status<int> s, int n) : runner(std::move(r)), status(std::move(s)), nobs(n) {}
    tsne runner;
    tsne::Status<int> status;
    int nobs;
};

//[[Rcpp::export(rng=false)]]
SEXP initialize_tsne(SEXP nnptr, double perplexity, int nthreads) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    KnncollePtr nns(nnptr);

    tsne runner;
    runner.set_perplexity(perplexity).set_max_depth(7).set_interpolation(100);
    auto status = runner.initialize(nns.get());
   
    return Rcpp::XPtr<InitializedTsne>(new InitializedTsne(std::move(runner), std::move(status), nns->nobs()));
}

//[[Rcpp::export(rng=false)]]
SEXP run_tsne(SEXP init, int nthreads) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif

    Rcpp::XPtr<InitializedTsne> ptr(init);
    std::vector<float> embedding(2 * ptr->nobs);
    qdtsne::initialize_random(embedding.data(), ptr->nobs);

    auto t = ptr->runner;
    t.run(ptr->status, embedding.data());

    Rcpp::NumericMatrix output(2, ptr->nobs);
    std::copy(embedding.begin(), embedding.end(), output.begin());
    return output;
}


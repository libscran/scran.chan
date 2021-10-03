#include "Rcpp.h"
#include "knncolle.h"
#include "qdtsne/qdtsne.hpp"
#ifdef _OPENMP
#include "omp.h"
#endif

struct InitializedTsne {
    InitializedTsne(qdtsne::Tsne<> r, qdtsne::Tsne<>::Status<int> s, int n) : runner(std::move(r)), status(std::move(s)), nobs(n) {}
    qdtsne::Tsne<> runner;
    qdtsne::Tsne<>::Status<int> status;
    int nobs;
};

//[[Rcpp::export(rng=false)]]
SEXP initialize_tsne(SEXP nnptr, double perplexity, int nthreads) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    KnncollePtr nns(nnptr);

    qdtsne::Tsne<> runner;
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
    Rcpp::NumericMatrix output(2, ptr->nobs);
    auto optr = static_cast<double*>(output.begin());
    qdtsne::initialize_random(optr, ptr->nobs);

    auto t = ptr->runner;
    t.run(ptr->status, optr);
    return output;
}


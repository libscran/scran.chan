#include "knncolle.h"
#include "knncolle/Annoy/Annoy.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

// [[Rcpp::export(rng=false)]]
SEXP build_nn_index(Rcpp::NumericMatrix data) {
    size_t nr = data.nrow(), nc = data.ncol();
    auto ptr = static_cast<const double*>(data.begin());
    return KnncollePtr(new knncolle::AnnoyEuclidean<int, float>(nr, nc, ptr));
}

// [[Rcpp::export(rng=false)]]
Rcpp::List find_nearest_neighbors(SEXP index, int k, int nthreads) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif

    KnncollePtr nns(index);
    auto nnptr = nns.get();
    size_t nobs = nnptr->nobs();

    Rcpp::NumericMatrix distmat(k, nobs);
    double* dptr = static_cast<double*>(distmat.begin());
    Rcpp::IntegerMatrix idxmat(k, nobs);
    int* iptr = static_cast<int*>(idxmat.begin());

    #pragma omp parallel for
    for (size_t o = 0; o < nobs; ++o) {
        auto res = nnptr->find_nearest_neighbors(o, k);
        auto dcopy = dptr + o * k;
        auto icopy = iptr + o * k;
        for (auto rIt = res.begin(); rIt != res.end(); ++rIt, ++dcopy, ++icopy) {
            *icopy = rIt->first;
            *dcopy = rIt->second;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("index") = idxmat,
        Rcpp::Named("distance") = distmat
    );
}

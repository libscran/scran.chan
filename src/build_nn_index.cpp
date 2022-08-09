#include "config.h"

#include "Rcpp.h"
#include "knncolle.h"
#include "knncolle/knncolle.hpp"

// [[Rcpp::export(rng=false)]]
SEXP build_nn_index(Rcpp::NumericMatrix data, bool approximate) {
    size_t nr = data.nrow(), nc = data.ncol();
    auto ptr = static_cast<const double*>(data.begin());
    if (approximate) {
        return KnncollePtr(new knncolle::AnnoyEuclidean<int, float>(nr, nc, ptr));
    } else {
        return KnncollePtr(new knncolle::VpTreeEuclidean<int, float>(nr, nc, ptr));
    }
}

// [[Rcpp::export(rng=false)]]
Rcpp::List find_nearest_neighbors(SEXP index, int k, int nthreads) {
    KnncollePtr nns(index);
    auto nnptr = nns.get();
    size_t nobs = nnptr->nobs();

    Rcpp::NumericMatrix distmat(k, nobs);
    double* dptr = static_cast<double*>(distmat.begin());
    Rcpp::IntegerMatrix idxmat(k, nobs);
    int* iptr = static_cast<int*>(idxmat.begin());

    #pragma omp parallel for num_threads(nthreads)
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

// [[Rcpp::export(rng=false)]]
Rcpp::List find_nearest_neighbor_indices(SEXP index, int k, int nthreads) {
    KnncollePtr nns(index);
    auto nnptr = nns.get();
    size_t nobs = nnptr->nobs();

    Rcpp::IntegerMatrix idxmat(k, nobs);
    int* iptr = static_cast<int*>(idxmat.begin());

    #pragma omp parallel for num_threads(nthreads)
    for (size_t o = 0; o < nobs; ++o) {
        auto res = nnptr->find_nearest_neighbors(o, k);
        auto icopy = iptr + o * k;
        for (auto rIt = res.begin(); rIt != res.end(); ++rIt, ++icopy) {
            *icopy = rIt->first;
        }
    }

    return Rcpp::List::create(Rcpp::Named("index") = idxmat);
}

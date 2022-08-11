#ifndef UTILS_KNNCOLLE_H
#define UTILS_KNNCOLLE_H

#include "knncolle/utils/Base.hpp"
#include "Rcpp.h"

typedef Rcpp::XPtr<knncolle::Base<int, double> > KnncollePtr;

#include <vector>

template<typename Idx = int, class Dist = double>
std::vector<std::vector<std::pair<Idx, Dist> > > unpack_neighbors(Rcpp::IntegerMatrix nnidx, Rcpp::NumericMatrix nndist) {
    size_t nobs = nnidx.cols();
    size_t nneighbors = nnidx.rows();
    const int* iptr = static_cast<const int*>(nnidx.begin());
    const double* dptr = static_cast<const double*>(nndist.begin());

    std::vector<std::vector<std::pair<Idx, Dist> > > neighbors(nobs);
    for (auto& current : neighbors) {
        current.resize(nneighbors);
        for (auto cIt = current.begin(); cIt != current.end(); ++cIt, ++iptr, ++dptr) {
            cIt->first = *iptr;
            cIt->second = *dptr;
        }
    }

    return neighbors;
}

template<typename Idx = int>
std::vector<std::vector<Idx> > unpack_neighbors(Rcpp::IntegerMatrix nnidx) {
    size_t nobs = nnidx.cols();
    size_t nneighbors = nnidx.rows();
    const int* iptr = static_cast<const int*>(nnidx.begin());

    std::vector<std::vector<Idx> > neighbors(nobs);
    for (auto& current : neighbors) {
        current.resize(nneighbors);
        for (auto cIt = current.begin(); cIt != current.end(); ++cIt, ++iptr) {
            *cIt = *iptr;
        }
    }

    return neighbors;
}

#endif

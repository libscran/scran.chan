#include "config.h"

#include "Rcpp.h"
#include "knncolle.h"
#include "knncolle/Annoy/Annoy.hpp"

// [[Rcpp::export(rng=false)]]
SEXP build_nn_index(Rcpp::NumericMatrix data) {
    size_t nr = data.nrow(), nc = data.ncol();
    auto ptr = static_cast<const double*>(data.begin());
    return KnncollePtr(new knncolle::AnnoyEuclidean<int, float>(nr, nc, ptr));
}

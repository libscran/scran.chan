#include "scran/normalization/LogNormCounts.hpp"
#include "Rcpp.h"
#include "tatamize.h"

//[[Rcpp::export(rng=false)]]
SEXP log_norm_counts(SEXP x, Rcpp::Nullable<Rcpp::NumericVector> size_factors) {
    auto mat = extract_NumericMatrix_shared(x);
    scran::LogNormCounts norm;

    std::vector<double> sf;
    if (size_factors.isNotNull()) {
        Rcpp::NumericVector sf_(size_factors);
        sf = std::vector<double>(sf_.begin(), sf_.end());
    } else {
        sf = tatami::column_sums(mat.get());
    }

    return new_MatrixChan(norm.run(std::move(mat), std::move(sf)));
}

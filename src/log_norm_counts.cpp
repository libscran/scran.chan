#include "scran/normalization/LogNormCounts.hpp"
#include "Rcpp.h"
#include "tatamize.h"
#include "ResolvedBatch.h"

//[[Rcpp::export(rng=false)]]
SEXP log_norm_counts(SEXP x, Rcpp::Nullable<Rcpp::NumericVector> size_factors, Rcpp::Nullable<Rcpp::IntegerVector> batch, std::string batch_mode) {
    auto mat = extract_NumericMatrix_shared(x);

    scran::LogNormCounts norm;
    if (batch_mode == "perblock") {
        norm.set_block_mode(scran::CenterSizeFactors::PER_BLOCK);
    }

    std::vector<double> sf;
    if (size_factors.isNotNull()) {
        Rcpp::NumericVector sf_(size_factors);
        sf = std::vector<double>(sf_.begin(), sf_.end());
    } else {
        sf = tatami::column_sums(mat.get());
    }

    auto batch_info = ResolvedBatch(batch);
    auto bptr = batch_info.ptr;

    return new_MatrixChan(norm.run_blocked(std::move(mat), std::move(sf), bptr));
}

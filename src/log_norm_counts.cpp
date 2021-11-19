#include "scran/normalization/LogNormCounts.hpp"
#include "Rcpp.h"
#include "tatamize.h"
#include "ResolvedBatch.h"

//[[Rcpp::export(rng=false)]]
SEXP log_norm_counts(SEXP x, Rcpp::Nullable<Rcpp::NumericVector> size_factors, Rcpp::Nullable<Rcpp::IntegerVector> batch, std::string batch_mode) {
    auto mat = extract_NumericMatrix_shared(x);

    std::vector<double> sf;
    if (size_factors.isNotNull()) {
        Rcpp::NumericVector sf_(size_factors);
        sf = std::vector<double>(sf_.begin(), sf_.end());
    } else {
        sf = tatami::column_sums(mat.get());
    }

    auto batch_info = ResolvedBatch(batch);
    auto bptr = batch_info.ptr;

    scran::CenterSizeFactors cenc;
    if (batch_mode == "perblock") {
        cenc.set_block_mode(scran::CenterSizeFactors::PER_BLOCK);
    }
    cenc.run_blocked(sf.size(), sf.data(), bptr);
    Rcpp::NumericVector centered(sf.begin(), sf.end());

    scran::LogNormCounts norm;
    norm.set_center(false);

    return Rcpp::List::create(
        // The blocking doesn't really have much effect at this point, as the
        // centering was already done; but we just throw it in, just in case
        // some future behavior relies on knowledge of blocking.
        Rcpp::Named("pointer") = new_MatrixChan(norm.run_blocked(std::move(mat), std::move(sf), bptr)),

        Rcpp::Named("size_factors") = centered
    );
}

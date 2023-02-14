#include "config.h"

#include "Rcpp.h"

#include "scran/aggregation/ScoreFeatureSet.hpp"

#include "tatamize.h"
#include "ResolvedFeatures.h"
#include "ResolvedBatch.h"

// [[Rcpp::export(rng=false)]]
Rcpp::List score_feature_set(SEXP x, Rcpp::LogicalVector features, Rcpp::Nullable<Rcpp::IntegerVector> batch, bool scale, int nthreads) {
    scran::ScoreFeatureSet scorer;
    scorer.set_num_threads(nthreads);
    scorer.set_scale(scale);

    ResolvedFeatures feat(features);
    const int* fptr = feat.ptr;

    auto batch_info = ResolvedBatch(batch);
    auto bptr = batch_info.ptr;

    auto mat = extract_NumericMatrix(x);
    auto res = scorer.run_blocked(mat, fptr, bptr);

    return Rcpp::List::create(
        Rcpp::Named("scores") = Rcpp::NumericVector(res.scores.begin(), res.scores.end()),
        Rcpp::Named("weights") = Rcpp::NumericVector(res.weights.begin(), res.weights.end())
    );
}

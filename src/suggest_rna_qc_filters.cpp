#include "config.h"

#include "Rcpp.h"
#include "scran/quality_control/SuggestRnaQcFilters.hpp"
#include "ResolvedBatch.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List suggest_rna_qc_filters(Rcpp::NumericVector sums, Rcpp::IntegerVector detected, Rcpp::List subsets, Rcpp::Nullable<Rcpp::IntegerVector> batch, double nmads) {
    std::vector<Rcpp::NumericVector> in_subsets;
    const size_t nsubs = subsets.size();
    for (auto sIt = subsets.begin(); sIt != subsets.end(); ++sIt) {
        in_subsets.emplace_back(*sIt);
    }

    scran::PerCellRnaQcMetrics::Buffers<double, int> buffers;
    buffers.sums = sums.begin();
    buffers.detected = detected.begin();
    for (auto& s : in_subsets) {
        buffers.subset_proportions.push_back(s.begin());
    }

    auto batch_info = ResolvedBatch(batch);
    auto bptr = batch_info.ptr;
    size_t ncells = sums.size();

    // Setting up the outputs.
    scran::SuggestRnaQcFilters qc;
    qc.set_num_mads(nmads);
    auto thresh = qc.run_blocked(ncells, bptr, buffers);

    // Converting threshold statistics.
    Rcpp::NumericVector thresh_sums(thresh.sums.begin(), thresh.sums.end());
    Rcpp::NumericVector thresh_detected(thresh.detected.begin(), thresh.detected.end());
    Rcpp::List thresh_subs(nsubs);
    for (size_t s = 0; s < nsubs; ++s) {
        const auto& current = thresh.subset_proportions[s];
        thresh_subs[s] = Rcpp::NumericVector(current.begin(), current.end());
    }

    // Computing the filters.
    Rcpp::LogicalVector by_overall(ncells);
    thresh.filter_blocked(ncells, bptr, buffers, static_cast<int*>(by_overall.begin()));

    return Rcpp::List::create(
        Rcpp::Named("filter") = by_overall,
        Rcpp::Named("thresholds") = Rcpp::List::create(
            Rcpp::Named("sums") = thresh_sums,
            Rcpp::Named("detected") = thresh_detected,
            Rcpp::Named("subsets") = thresh_subs
        )
    );
}

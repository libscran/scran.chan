#include "config.h"

#include "Rcpp.h"
#include "scran/quality_control/PerCellRnaQcFilters.hpp"
#include "ResolvedBatch.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List per_cell_qc_filters(Rcpp::NumericVector sums, Rcpp::IntegerVector detected, Rcpp::List subsets, Rcpp::Nullable<Rcpp::IntegerVector> batch, double nmads) {
    std::vector<const double*> in_sub_ptrs;
    std::vector<Rcpp::NumericVector> in_subsets;
    const size_t nsubs = subsets.size();

    for (auto sIt = subsets.begin(); sIt != subsets.end(); ++sIt) {
        in_subsets.emplace_back(*sIt);
    }
    for (const auto& s : in_subsets) {
        in_sub_ptrs.push_back(s.begin());
    }

    auto batch_info = ResolvedBatch(batch);
    auto bptr = batch_info.ptr;

    // Setting up the outputs.
    const size_t N = sums.size();
    Rcpp::LogicalVector by_sums(N), by_detected(N), by_overall(N);

    std::vector<Rcpp::LogicalVector> by_subsets;
    by_subsets.reserve(in_subsets.size());
    std::vector<int*> out_sub_ptrs;

    for (size_t i = 0; i < in_subsets.size(); ++i) {
        by_subsets.push_back(Rcpp::LogicalVector(N));
        out_sub_ptrs.push_back(by_subsets.back().begin());
    }

    scran::PerCellQCFilters qc;
    qc.set_nmads(nmads);
    auto thresh = qc.run_blocked(N, 
        bptr,
        static_cast<const double*>(sums.begin()),
        static_cast<const int*>(detected.begin()),
        std::move(in_sub_ptrs),
        static_cast<int*>(by_sums.begin()),
        static_cast<int*>(by_detected.begin()),
        std::move(out_sub_ptrs),
        static_cast<int*>(by_overall.begin())
    );

    // Converting threshold statistics.
    Rcpp::NumericVector thresh_sums(thresh.sums.begin(), thresh.sums.end());
    Rcpp::NumericVector thresh_detected(thresh.detected.begin(), thresh.detected.end());
    Rcpp::List thresh_subs(nsubs);
    for (size_t s = 0; s < nsubs; ++s) {
        const auto& current = thresh.subset_proportions[s];
        thresh_subs[s] = Rcpp::NumericVector(current.begin(), current.end());
    }

    return Rcpp::List::create(
        Rcpp::Named("filters") = Rcpp::List::create(
            Rcpp::Named("sums") = by_sums,
            Rcpp::Named("detected") = by_detected,
            Rcpp::Named("subsets") = by_subsets,
            Rcpp::Named("overall") = by_overall
        ),
        Rcpp::Named("thresholds") = Rcpp::List::create(
            Rcpp::Named("sums") = thresh_sums,
            Rcpp::Named("detected") = thresh_detected,
            Rcpp::Named("subsets") = thresh_subs
        )
    );
}

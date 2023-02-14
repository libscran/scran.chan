#include "config.h"

#include "Rcpp.h"
#include "scran/quality_control/SuggestCrisprQcFilters.hpp"
#include "ResolvedBatch.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List suggest_crispr_qc_filters(Rcpp::NumericVector sums, Rcpp::NumericVector max_prop, Rcpp::Nullable<Rcpp::IntegerVector> batch, double nmads) {
    scran::PerCellCrisprQcMetrics::Buffers<double, int> buffers;
    buffers.sums = sums.begin();
    buffers.max_proportion = max_prop.begin();

    auto batch_info = ResolvedBatch(batch);
    auto bptr = batch_info.ptr;
    size_t ncells = sums.size();

    // Setting up the outputs.
    scran::SuggestCrisprQcFilters qc;
    qc.set_num_mads(nmads);
    auto thresh = qc.run_blocked(ncells, bptr, buffers);

    // Converting threshold statistics.
    Rcpp::NumericVector thresh_max(thresh.max_count.begin(), thresh.max_count.end());

    // Computing the filters.
    Rcpp::LogicalVector by_overall(ncells);
    thresh.filter_blocked(ncells, bptr, buffers, static_cast<int*>(by_overall.begin()));

    return Rcpp::List::create(
        Rcpp::Named("filter") = by_overall,
        Rcpp::Named("thresholds") = Rcpp::List::create(
            Rcpp::Named("max.count") = thresh_max
        )
    );
}

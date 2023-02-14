#include "config.h"

#include "Rcpp.h"
#include "scran/quality_control/PerCellAdtQcMetrics.hpp"
#include "tatamize.h"

// [[Rcpp::export(rng=false)]]
Rcpp::List per_cell_adt_qc_metrics(SEXP x, Rcpp::List subsets, int nthreads) {
    std::vector<const int*> in_con_ptrs;
    std::vector<Rcpp::LogicalVector> in_subsets;

    for (auto sIt = subsets.begin(); sIt != subsets.end(); ++sIt) {
        in_subsets.emplace_back(*sIt);
    }
    for (const auto& s : in_subsets) {
        in_con_ptrs.push_back(s.begin());
    }

    auto mat = extract_NumericMatrix(x);

    // Creating output containers.
    size_t nc = mat->ncol();
    Rcpp::NumericVector sums(nc);
    Rcpp::IntegerVector detected(nc);
    std::vector<Rcpp::NumericVector> out_subsets;
    for (size_t s = 0; s < subsets.size(); ++s) {
        out_subsets.emplace_back(nc);
    }

    scran::PerCellAdtQcMetrics::Buffers<double, int> buffers;
    buffers.sums = sums.begin();
    buffers.detected = detected.begin();
    for (auto& s : out_subsets) {
        buffers.subset_totals.push_back(s.begin());
    }

    // Running QC code.
    scran::PerCellAdtQcMetrics qc;
    qc.set_num_threads(nthreads);
    qc.run(mat, std::move(in_con_ptrs), buffers);

    return Rcpp::List::create(
        Rcpp::Named("sums") = sums,
        Rcpp::Named("detected") = detected,
        Rcpp::Named("subsets") = Rcpp::List(out_subsets.begin(), out_subsets.end())
    );
}

#include "config.h"

#include "Rcpp.h"
#include "scran/quality_control/PerCellCrisprQcMetrics.hpp"
#include "tatamize.h"

// [[Rcpp::export(rng=false)]]
Rcpp::List per_cell_crispr_qc_metrics(SEXP x, int nthreads) {
    auto mat = extract_NumericMatrix(x);

    // Creating output containers.
    size_t nc = mat->ncol();
    Rcpp::NumericVector sums(nc);
    Rcpp::IntegerVector detected(nc);
    Rcpp::NumericVector max_proportion(nc);
    Rcpp::IntegerVector max_index(nc);

    scran::PerCellCrisprQcMetrics::Buffers<double, int> buffers;
    buffers.sums = sums.begin();
    buffers.detected = detected.begin();
    buffers.max_proportion = max_proportion.begin();
    buffers.max_index = max_index.begin();

    // Running QC code.
    scran::PerCellCrisprQcMetrics qc;
    qc.set_num_threads(nthreads);
    qc.run(mat, buffers);

    return Rcpp::List::create(
        Rcpp::Named("sums") = sums,
        Rcpp::Named("detected") = detected,
        Rcpp::Named("max.proportion") = max_proportion,
        Rcpp::Named("max.index") = max_index
    );
}

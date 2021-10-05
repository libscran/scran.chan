#undef _OPENMP
#include "Rcpp.h"
#include "scran/quality_control/PerCellQCMetrics.hpp"
#include "tatamize.h"
#ifdef _OPENMP
#include "omp.h"
#endif

// [[Rcpp::export(rng=false)]]
Rcpp::List per_cell_qc_metrics(SEXP x, Rcpp::List subsets, int nthreads) {
#ifdef _openmp
    omp_set_num_threads(nthreads);
#endif
    
    std::vector<const int*> in_sub_ptrs;
    std::vector<Rcpp::LogicalVector> in_subsets;

    for (auto sIt = subsets.begin(); sIt != subsets.end(); ++sIt) {
        in_subsets.emplace_back(*sIt);
    }
    for (const auto& s : in_subsets) {
        in_sub_ptrs.push_back(s.begin());
    }

    auto mat = extract_NumericMatrix(x);

    // Creating output containers.
    size_t nc = mat->ncol();
    Rcpp::NumericVector sums(nc);
    Rcpp::IntegerVector detected(nc);
    std::vector<Rcpp::NumericVector> out_subsets;
    std::vector<double*> out_sub_ptrs;
    for (size_t s = 0; s < subsets.size(); ++s) {
        out_subsets.emplace_back(nc);
    }
    for (auto& s : out_subsets) {
        out_sub_ptrs.push_back(s.begin());
    }

    // Running QC code.
    scran::PerCellQCMetrics qc;
    qc.run(mat, 
        std::move(in_sub_ptrs), 
        static_cast<double*>(sums.begin()),
        static_cast<int*>(detected.begin()),
        std::move(out_sub_ptrs)
    );

    return Rcpp::List::create(
        Rcpp::Named("sums") = sums,
        Rcpp::Named("detected") = detected,
        Rcpp::Named("subsets") = Rcpp::List(out_subsets.begin(), out_subsets.end())
    );
}

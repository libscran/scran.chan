#include "config.h"

#include "tatamize.h"
#include "scran/quality_control/FilterCells.hpp"
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
SEXP filter_cells(SEXP x, Rcpp::LogicalVector discard) {
    scran::FilterCells qc;
    auto y = qc.run(extract_NumericMatrix_shared(x), static_cast<const int*>(discard.begin()));
    return new_MatrixChan(std::move(y));
}

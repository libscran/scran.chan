#include "config.h"

#include "Rcpp.h"
#include "tatamize.h"
#include "ResolvedBatch.h"

template<int MARGIN, class T>
void subset_dimensions(T& mat, Rcpp::RObject sub) {
    if (sub == R_NilValue) {
        return;
    }

    Rcpp::IntegerVector sub_(sub);
    bool is_consecutive = true;
    for (size_t i = 1, end = sub_.size(); i < end; ++i) {
        if (sub_[i] != sub_[i] + 1) {
            is_consecutive = false;
        }
    } 

    T output;
    if (is_consecutive) {
        if (sub_.size()) {
            output = tatami::make_DelayedSubsetBlock<MARGIN>(std::move(mat), sub_[0], sub_[sub_.size() - 1] + 1);
        } else {
            output = tatami::make_DelayedSubsetBlock<MARGIN>(std::move(mat), 0, 0);
        }
    } else {
        std::vector<int> copy(sub_.begin(), sub_.end());
        for (auto& x : copy) {
            --x; // 1-based to 0-based.
        }
        output = tatami::make_DelayedSubset<MARGIN>(std::move(mat), std::move(copy));
    }

    mat = output;
}

//[[Rcpp::export(rng=false)]]
SEXP subset_matrix(SEXP x, Rcpp::RObject i, Rcpp::RObject j) {
    auto mat = extract_NumericMatrix_shared(x);
    subset_dimensions<0>(mat, i);
    subset_dimensions<1>(mat, j);
    return new_MatrixChan(std::move(mat));
}

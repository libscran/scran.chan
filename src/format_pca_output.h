#ifndef FORMAT_PCA_OUTPUT_H
#define FORMAT_PCA_OUTPUT_H

#include "Rcpp.h"
#include <algorithm>

template<class Mat, class Result>
Rcpp::List format_pca_output(Mat mat, const Result& res) {
    size_t ndim = res.pcs.rows();

    size_t ncol = mat->ncol();
    Rcpp::NumericMatrix output_pcs(ndim, ncol);
    std::copy(res.pcs.data(), res.pcs.data() + ndim * ncol, output_pcs.begin());

    Rcpp::NumericVector output_prop(res.variance_explained.begin(), res.variance_explained.end());
    for (auto& p : output_prop) {
        p /= res.total_variance;
    }

    return Rcpp::List::create(
        Rcpp::Named("components") = output_pcs, 
        Rcpp::Named("prop.variance") = output_prop
    );
}

#endif

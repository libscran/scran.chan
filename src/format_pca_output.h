#ifndef FORMAT_PCA_OUTPUT_H
#define FORMAT_PCA_OUTPUT_H

#include "Rcpp.h"
#include <algorithm>

template<class Result>
Rcpp::List format_pca_output(size_t ncol, const Result& res, bool rotation) {
    size_t ndim = res.pcs.rows();

    Rcpp::NumericMatrix output_pcs(ndim, ncol);
    std::copy(res.pcs.data(), res.pcs.data() + ndim * ncol, output_pcs.begin());

    Rcpp::RObject output_rot;
    if (rotation) {
        size_t ngenes = res.rotation.rows();
        Rcpp::NumericMatrix tmp(ngenes, ndim);
        std::copy(res.rotation.data(), res.rotation.data() + ngenes * ndim, tmp.begin());
        output_rot = tmp;
    } else {
        output_rot = R_NilValue;
    }

    Rcpp::NumericVector output_prop(res.variance_explained.begin(), res.variance_explained.end());
    for (auto& p : output_prop) {
        p /= res.total_variance;
    }

    if (rotation) {
        return Rcpp::List::create(
            Rcpp::Named("components") = output_pcs, 
            Rcpp::Named("rotation") = output_rot, 
            Rcpp::Named("prop.variance") = output_prop
        );
    } else {
        return Rcpp::List::create(
            Rcpp::Named("components") = output_pcs, 
            Rcpp::Named("prop.variance") = output_prop
        );
    }
}

#endif

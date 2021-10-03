#include "Rcpp.h"
#include "scran/dimensionality_reduction/RunPCA.hpp"
#include "tatamize.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List run_pca(SEXP x, int ndim, Rcpp::Nullable<Rcpp::LogicalVector> features) {
    scran::RunPCA pcs;
    pcs.set_rank(ndim);

    const int* fptr = NULL;
    Rcpp::LogicalVector feats_;
    if (features.isNotNull()) {
        feats_ = Rcpp::LogicalVector(features);
        fptr = static_cast<const int*>(feats_.begin());
    }

    auto mat = extract_NumericMatrix_shared(x);
    auto res = pcs.run(mat, fptr);
    res.pcs.adjointInPlace();

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

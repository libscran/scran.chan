#include "config.h"

#include "Rcpp.h"

#include "scran/dimensionality_reduction/RunPCA.hpp"

#include "tatamize.h"
#include "ResolvedFeatures.h"
#include "format_pca_output.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List run_pca(SEXP x, int ndim, Rcpp::Nullable<Rcpp::LogicalVector> features, bool rotation, int nthreads) {
    scran::RunPCA pcs;
    pcs.set_rank(ndim).set_num_threads(nthreads);

    ResolvedFeatures feat(features);
    const int* fptr = feat.ptr;

    auto mat = extract_NumericMatrix(x);
    auto res = pcs.run(mat, fptr);

    return format_pca_output(mat->ncol(), res, rotation);
}

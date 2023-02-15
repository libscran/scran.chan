#include "config.h"

#include "Rcpp.h"

#include "scran/dimensionality_reduction/MultiBatchPCA.hpp"

#include "tatamize.h"
#include "ResolvedFeatures.h"
#include "format_pca_output.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List run_multibatch_pca(SEXP x, int ndim, Rcpp::IntegerVector batch, Rcpp::Nullable<Rcpp::LogicalVector> features, bool scale, bool rotation, int nthreads) {
    scran::MultiBatchPCA pcs;
    pcs
        .set_rank(ndim)
        .set_scale(scale)
        .set_num_threads(nthreads);

    ResolvedFeatures feat(features);
    const int* fptr = feat.ptr;

    auto mat = extract_NumericMatrix(x);
    auto res = pcs.run(mat, static_cast<const int*>(batch.begin()), fptr);

    return format_pca_output(mat->ncol(), res, rotation);
}

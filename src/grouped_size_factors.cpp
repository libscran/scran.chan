#include "config.h"

#include "Rcpp.h"
#include "scran/normalization/GroupedSizeFactors.hpp"
#include "tatamize.h"

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector grouped_size_factors(SEXP x, Rcpp::IntegerVector clusters, bool center, double prior_count, int reference, int nthreads) {
    scran::GroupedSizeFactors runner;
    runner.set_center(center).set_prior_count(prior_count).set_num_threads(nthreads);

    auto mat = extract_NumericMatrix(x);
    Rcpp::NumericVector output(mat->ncol());
    double* optr = output.begin();

    const int* cptr = clusters.begin();
    if (reference >= 0) {
        runner.run(mat, cptr, reference, optr);
    } else {
        runner.run(mat, cptr, optr);
    }

    return output;
}

#include "config.h"

#include "tatamize.h"
#include "Rcpp.h"

#include "scran/feature_selection/ModelGeneVar.hpp"
#include "scran/utils/average_vectors.hpp"

#include "ResolvedBatch.h"
#include "vector_to_pointers.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List model_gene_var(SEXP x, Rcpp::Nullable<Rcpp::IntegerVector> batch, double span, bool use_fixed, double fixed_width, int min_count, int nthreads) {
    scran::ModelGeneVar mvar;
    mvar.set_num_threads(nthreads);
    mvar.set_span(span);
    mvar.set_use_fixed_width(use_fixed);
    mvar.set_fixed_width(fixed_width);
    mvar.set_minimum_window_count(min_count);

    auto mat = extract_NumericMatrix(x);
    size_t NR = mat->nrow(), NC = mat->ncol();

    auto batch_info = ResolvedBatch(batch);
    auto bptr = batch_info.ptr;
    size_t nbatches = batch_info.number();

    // Setting up the outputs. 
    std::vector<Rcpp::NumericVector> means, variances, fitted, residuals;
    means.reserve(nbatches);
    variances.reserve(nbatches);
    fitted.reserve(nbatches);
    residuals.reserve(nbatches);
    
    for (size_t b = 0; b < nbatches; ++b) {
        means.emplace_back(NR);
        variances.emplace_back(NR);
        fitted.emplace_back(NR);
        residuals.emplace_back(NR);
    }

    mvar.run_blocked(mat, 
        bptr,
        vector_to_pointers<double>(means),
        vector_to_pointers<double>(variances),
        vector_to_pointers<double>(fitted),
        vector_to_pointers<double>(residuals)
    );

    auto createDF = [](Rcpp::NumericVector m, Rcpp::NumericVector v, Rcpp::NumericVector f, Rcpp::NumericVector r) -> Rcpp::DataFrame {
        return Rcpp::DataFrame::create(
            Rcpp::Named("means") = m, 
            Rcpp::Named("variances") = v,
            Rcpp::Named("fitted") = f,
            Rcpp::Named("residuals") = r
        );
    };

    // Deciding what to return.
    if (bptr == NULL) {
        return Rcpp::List::create(
            Rcpp::Named("statistics") = createDF(means[0], variances[0], fitted[0], residuals[0])
        );

    } else {
        Rcpp::List per_batch(nbatches);
        for (size_t b = 0; b < nbatches; ++b) {
            per_batch[b] = createDF(means[b], variances[b], fitted[b], residuals[b]);
        }
        
        Rcpp::NumericVector ave_mean(NR), ave_var(NR), ave_fit(NR), ave_res(NR);
        scran::average_vectors(NR, vector_to_pointers<double>(means), static_cast<double*>(ave_mean.begin()));
        scran::average_vectors(NR, vector_to_pointers<double>(variances), static_cast<double*>(ave_var.begin()));
        scran::average_vectors(NR, vector_to_pointers<double>(fitted), static_cast<double*>(ave_fit.begin()));
        scran::average_vectors(NR, vector_to_pointers<double>(residuals), static_cast<double*>(ave_res.begin()));

        return Rcpp::List::create(
            Rcpp::Named("statistics") = createDF(ave_mean, ave_var, ave_fit, ave_res),
            Rcpp::Named("per.batch") = per_batch
        );
    }
}

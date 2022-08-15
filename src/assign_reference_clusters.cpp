#include "config.h"

#include "knncolle.h"
#include "scran/clustering/AssignReferenceClusters.hpp"
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
SEXP assign_reference_clusters(SEXP ref_index, Rcpp::IntegerVector ref_clusters, Rcpp::NumericMatrix test_data, int k, bool approximate, int nthreads) {
    KnncollePtr index(ref_index);
    size_t ref_nr = index->ndim(), ref_nc = index->nobs();

    size_t test_nr = test_data.nrow(), test_nc = test_data.ncol();
    auto test_ptr = static_cast<const double*>(test_data.begin());
    if (test_nr != ref_nr) {
        throw std::runtime_error("test and reference matrices should have the same number of rows");
    }

    scran::AssignReferenceClusters runner;
    runner
        .set_num_neighbors(k)
        .set_approximate(approximate)
        .set_num_threads(nthreads);

    Rcpp::IntegerVector output(test_nc);
    Rcpp::NumericVector best_prop(test_nc);
    Rcpp::NumericVector second_prop(test_nc);
    runner.run(
        index.get(),
        static_cast<const int*>(ref_clusters.begin()), 
        test_nc, 
        test_ptr, 
        static_cast<int*>(output.begin()),
        static_cast<double*>(best_prop.begin()),
        static_cast<double*>(second_prop.begin())
    );
   
    return Rcpp::List::create(
        Rcpp::Named("assigned") = output,
        Rcpp::Named("best.prop") = best_prop,
        Rcpp::Named("second.prop") = second_prop
    );
}

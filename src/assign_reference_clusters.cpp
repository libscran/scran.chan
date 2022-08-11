#include "config.h"

#include "scran/clustering/AssignReferenceClusters.hpp"
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
SEXP assign_reference_clusters(Rcpp::NumericMatrix ref_data, Rcpp::IntegerVector ref_clusters, Rcpp::NumericMatrix test_data, double quantile, bool approximate, int nthreads) {
    size_t ref_nr = ref_data.nrow(), ref_nc = ref_data.ncol();
    auto ref_ptr = static_cast<const double*>(ref_data.begin());

    size_t test_nr = test_data.nrow(), test_nc = test_data.ncol();
    auto test_ptr = static_cast<const double*>(test_data.begin());
    if (test_nr != ref_nr) {
        throw std::runtime_error("test and reference matrices should have the same number of rows");
    }

    scran::AssignReferenceClusters runner;
    runner
        .set_quantile(quantile)
        .set_approximate(approximate)
        .set_num_threads(nthreads);

    Rcpp::IntegerVector output(test_nc);
    runner.run(ref_nr, ref_nc, ref_ptr, static_cast<const int*>(ref_clusters.begin()), test_nc, test_ptr, static_cast<int*>(output.begin()));
   
    return output; 
}

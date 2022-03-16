#include "mnncorrect/MnnCorrect.hpp"
#include "Rcpp.h"
#ifdef _OPENMP
#include "omp.h"
#endif

//[[Rcpp::export(rng=false)]]
Rcpp::List mnn_correct(Rcpp::NumericMatrix x, Rcpp::IntegerVector batch, int k, double nmads, int nthreads) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif

    mnncorrect::MnnCorrect<> runner;
    runner.set_num_neighbors(k).set_num_mads(nmads);

    Rcpp::NumericMatrix output(x.nrow(), x.ncol());
    auto res = runner.run(x.nrow(), x.ncol(), 
        static_cast<const double*>(x.begin()), 
        static_cast<const int*>(batch.begin()), 
        static_cast<double*>(output.begin()));

    return Rcpp::List::create(
        Rcpp::Named("corrected") = output,
        Rcpp::Named("merge.order") = Rcpp::IntegerVector(res.merge_order.begin(), res.merge_order.end()),
        Rcpp::Named("num.pairs") = Rcpp::IntegerVector(res.num_pairs.begin(), res.num_pairs.end())
    );
}

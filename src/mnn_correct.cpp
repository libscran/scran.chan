#include "config.h"

#include "mnncorrect/MnnCorrect.hpp"
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List mnn_correct(Rcpp::NumericMatrix x, Rcpp::IntegerVector batch, int k, double nmads, int nthreads, Rcpp::Nullable<Rcpp::IntegerVector> order, std::string ref_policy) {
    mnncorrect::MnnCorrect<> runner;
    runner.set_approximate(true).set_num_neighbors(k).set_num_mads(nmads);
    runner.set_num_threads(nthreads);

    std::vector<int> ordering;
    const int* optr = NULL;
    if (order.isNotNull()) {
        Rcpp::IntegerVector order_(order);
        ordering.insert(ordering.end(), order_.begin(), order_.end());
        optr = ordering.data();
    }

    if (ref_policy == "input") {
        runner.set_reference_policy(mnncorrect::Input);
    } else if (ref_policy == "max-variance") {
        runner.set_reference_policy(mnncorrect::MaxVariance);
    } else if (ref_policy == "max-rss") {
        runner.set_reference_policy(mnncorrect::MaxRss);
    } else if (ref_policy != "max-size") {
        throw std::runtime_error("unknown reference policy");
    }

    Rcpp::NumericMatrix output(x.nrow(), x.ncol());
    auto res = runner.run(x.nrow(), x.ncol(), 
        static_cast<const double*>(x.begin()), 
        static_cast<const int*>(batch.begin()), 
        static_cast<double*>(output.begin()),
        optr);

    return Rcpp::List::create(
        Rcpp::Named("corrected") = output,
        Rcpp::Named("merge.order") = Rcpp::IntegerVector(res.merge_order.begin(), res.merge_order.end()),
        Rcpp::Named("num.pairs") = Rcpp::IntegerVector(res.num_pairs.begin(), res.num_pairs.end())
    );
}

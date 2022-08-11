#include "config.h"

#include "scran/aggregation/DownsampleByNeighbors.hpp"
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
SEXP downsample_by_neighbors(Rcpp::NumericMatrix data, int k, bool approximate, int nthreads) {
    size_t nr = data.nrow(), nc = data.ncol();
    auto ptr = static_cast<const double*>(data.begin());

    scran::DownsampleByNeighbors runner;
    runner
        .set_num_neighbors(k)
        .set_approximate(approximate)
        .set_num_threads(nthreads);

    auto res = runner.run(nr, nc, ptr);
    return Rcpp::IntegerVector(res.begin(), res.end());
}

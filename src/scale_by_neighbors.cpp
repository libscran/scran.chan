#include "config.h"

#include "scran/dimensionality_reduction/ScaleByNeighbors.hpp"
#include "Rcpp.h"

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector scale_by_neighbors(Rcpp::List matrices, int k, bool approximate, int nthreads) {
    scran::ScaleByNeighbors runner;
    runner
        .set_neighbors(k)
        .set_approximate(approximate)
        .set_num_threads(nthreads);

    std::vector<std::pair<double, double> > values;
    values.reserve(matrices.size());

    for (size_t x = 0, end = matrices.size(); x < end; ++x) {
        Rcpp::NumericMatrix mat(matrices[x]);
        values.push_back(runner.compute_distance(mat.nrow(), mat.ncol(), mat.begin()));
    }

    auto output = scran::ScaleByNeighbors::compute_scale(values);
    return Rcpp::NumericVector(output.begin(), output.end());
}

#include "config.h"

#include "Rcpp.h"
#include "knncolle.h"
#include "umappp/Umap.hpp"

//[[Rcpp::export(rng=false)]]
SEXP run_umap(Rcpp::IntegerMatrix nnidx, Rcpp::NumericMatrix nndist, double min_dist, int seed, int nthreads) {
    auto neighbors = unpack_neighbors<int, float>(nnidx, nndist);
    size_t nobs = neighbors.size();

    umappp::Umap<float> runner;
    runner.set_min_dist(min_dist).set_seed(seed).set_num_threads(nthreads);
    std::vector<float> embedding(2 * nobs);
    runner.run(std::move(neighbors), 2, embedding.data());

    Rcpp::NumericMatrix output(2, nobs);
    std::copy(embedding.begin(), embedding.end(), output.begin());
    return output;
}

#include "config.h"

#include "Rcpp.h"
#include "knncolle.h"
#include "qdtsne/qdtsne.hpp"

//[[Rcpp::export(rng=false)]]
SEXP run_tsne(Rcpp::IntegerMatrix nnidx, Rcpp::NumericMatrix nndist, double perplexity, int interpolate, int max_depth, int max_iter, int seed, int nthreads) {
    auto neighbors = unpack_neighbors<int, float>(nnidx, nndist);
    size_t nobs = neighbors.size();

    qdtsne::Tsne<2, float> runner;
    runner
        .set_perplexity(perplexity)
        .set_max_depth(max_depth)
        .set_num_threads(nthreads)
        .set_max_iter(max_iter);
    
    if (interpolate) {
        if (interpolate > 0) {
            runner.set_interpolation(interpolate);
        } else if (nobs > 50000) { // deduce it automatically.
            runner.set_interpolation(200);
        }
    }

    std::vector<float> embedding(2 * nobs);
    qdtsne::initialize_random(embedding.data(), nobs, seed);
    runner.run(std::move(neighbors), embedding.data());

    Rcpp::NumericMatrix output(2, nobs);
    std::copy(embedding.begin(), embedding.end(), output.begin());
    return output;
}

//[[Rcpp::export(rng=false)]]
int perplexity_to_neighbors(double p) {
    return qdtsne::perplexity_to_k(p);
}

#include "config.h"

#include "Rcpp.h"
#include "knncolle.h"
#include "scran/clustering/BuildSNNGraph.hpp"
#include "scran/clustering/ClusterSNNGraph.hpp"

#include <string>

//[[Rcpp::export(rng=false)]]
SEXP cluster_snn_graph(Rcpp::IntegerMatrix nnidx, std::string weight_scheme, std::string method, double resolution, int steps, int seed, int nthreads) {
    auto neighbors = unpack_neighbors<int>(nnidx);
    auto ncells = neighbors.size();

    scran::BuildSNNGraph::Scheme weight;
    if (weight_scheme == "jaccard") {
        weight = scran::BuildSNNGraph::JACCARD;
    } else if (weight_scheme == "number") {
        weight = scran::BuildSNNGraph::NUMBER;
    } else if (weight_scheme == "rank") {
        weight = scran::BuildSNNGraph::RANKED;
    } else {
        throw std::runtime_error("unknown weighting scheme '" + weight_scheme + "'");
    }

    scran::BuildSNNGraph builder;
    builder.set_num_threads(nthreads).set_weighting_scheme(weight);
    auto edges = builder.run(neighbors);

    Rcpp::RObject output;

    if (method == "multilevel") {
        scran::ClusterSNNGraphMultiLevel runner;
        runner.set_resolution(resolution).set_seed(seed);
        auto res = runner.run(edges);

        Rcpp::List memberships(res.membership.size());
        for (size_t c = 0; c < memberships.size(); ++c) {
            const auto& current = res.membership[c];
            memberships[c] = Rcpp::IntegerVector(current.begin(), current.end());
        }
    
        output = Rcpp::List::create(
            Rcpp::Named("membership") = Rcpp::IntegerVector(memberships[res.max]),
            Rcpp::Named("best") = Rcpp::IntegerVector::create(res.max),
            Rcpp::Named("levels") = memberships,
            Rcpp::Named("modularity") = Rcpp::NumericVector(res.modularity.begin(), res.modularity.end())
        );

    } else if (method == "leiden") {
        scran::ClusterSNNGraphLeiden runner;
        runner.set_resolution(resolution).set_seed(seed);
        auto res = runner.run(edges);

        output = Rcpp::List::create(
            Rcpp::Named("membership") = Rcpp::IntegerVector(res.membership.begin(), res.membership.end()),
            Rcpp::Named("quality") = Rcpp::NumericVector::create(res.quality)
        );

    } else if (method == "walktrap") {
        scran::ClusterSNNGraphWalktrap runner;
        auto res = runner.run(edges);

        Rcpp::IntegerMatrix merges(res.merges.size(), 2);
        for (size_t m = 0; m < res.merges.size(); ++m) {
            const auto& step = res.merges[m];
            merges(m, 0) = step.first;
            merges(m, 1) = step.first;
        }

        output = Rcpp::List::create(
            Rcpp::Named("membership") = Rcpp::IntegerVector(res.membership.begin(), res.membership.end()),
            Rcpp::Named("merges") = merges,
            Rcpp::Named("modularity") = Rcpp::NumericVector(res.modularity.begin(), res.modularity.end())
        );

    } else {
        throw std::runtime_error("unknown community detection method '" + method + "'");
    }

    return output;
}

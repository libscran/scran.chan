#ifndef CLUSTER_SNN_GRAPH_H
#define CLUSTER_SNN_GRAPH_H
#include "scran/clustering/ClusterSNNGraph.hpp"
#include <deque>
#include "Rcpp.h"

struct SNNGraph {
    SNNGraph (size_t n, std::deque<scran::BuildSNNGraph::WeightedEdge> s, double r) : ncells(n), store(std::move(s)), resolution(r) {}

    auto run() {
        scran::ClusterSNNGraphMultiLevel runner;
        runner.set_resolution(resolution);
        return runner.run(ncells, store);
    }

    template<class V>
    Rcpp::List yield(const V& clusters) {
        Rcpp::List membership(clusters.membership.size());
        for (size_t c = 0; c < membership.size(); ++c) {
            const auto& current = clusters.membership[c];
            membership[c] = Rcpp::IntegerVector(current.begin(), current.end());
        }
    
        return Rcpp::List::create(
            Rcpp::Named("best") = Rcpp::IntegerVector::create(clusters.max),
            Rcpp::Named("membership") = membership,
            Rcpp::Named("modularity") = Rcpp::NumericVector(clusters.modularity.begin(), clusters.modularity.end())
        );
    }

    size_t ncells;
    std::deque<scran::BuildSNNGraph::WeightedEdge> store;
    double resolution;
};

#endif

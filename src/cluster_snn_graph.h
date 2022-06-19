#ifndef CLUSTER_SNN_GRAPH_H
#define CLUSTER_SNN_GRAPH_H

#include "scran/clustering/ClusterSNNGraph.hpp"
#include <deque>
#include "Rcpp.h"
#include <string>

struct SNNGraph {
    SNNGraph (size_t n, std::deque<scran::BuildSNNGraph::WeightedEdge> s, std::string m, double r) : 
        ncells(n), store(std::move(s)), method(m), resolution(r) {}

    void run() {
        if (method == "multilevel") {
            scran::ClusterSNNGraphMultiLevel runner;
            runner.set_resolution(resolution);
            res_multilevel = runner.run(ncells, store);
        } else if (method == "leiden") {
            scran::ClusterSNNGraphLeiden runner;
            runner.set_resolution(resolution);
            res_leiden = runner.run(ncells, store);
        } else {
            scran::ClusterSNNGraphWalktrap runner;
            res_walktrap = runner.run(ncells, store);
        }
        return;
    }

    Rcpp::List yield() {
        if (method == "multilevel") {
            const auto& clusters = res_multilevel;
            Rcpp::List memberships(clusters.membership.size());
            for (size_t c = 0; c < memberships.size(); ++c) {
                const auto& current = clusters.membership[c];
                memberships[c] = Rcpp::IntegerVector(current.begin(), current.end());
            }
        
            return Rcpp::List::create(
                Rcpp::Named("membership") = Rcpp::IntegerVector(memberships[clusters.max]),
                Rcpp::Named("best") = Rcpp::IntegerVector::create(clusters.max),
                Rcpp::Named("levels") = memberships,
                Rcpp::Named("modularity") = Rcpp::NumericVector(clusters.modularity.begin(), clusters.modularity.end())
            );
        } else if (method == "leiden") {
            const auto& clusters = res_leiden;
            return Rcpp::List::create(
                Rcpp::Named("membership") = Rcpp::IntegerVector(clusters.membership.begin(), clusters.membership.end()),
                Rcpp::Named("quality") = Rcpp::NumericVector::create(clusters.quality)
            );
        } else {
            const auto& clusters = res_walktrap;
            Rcpp::IntegerMatrix merges(clusters.merges.size(), 2);
            for (size_t m = 0; m < clusters.merges.size(); ++m) {
                const auto& step = clusters.merges[m];
                merges(m, 0) = step.first;
                merges(m, 1) = step.first;
            }

            return Rcpp::List::create(
                Rcpp::Named("membership") = Rcpp::IntegerVector(clusters.membership.begin(), clusters.membership.end()),
                Rcpp::Named("merges") = Rcpp::IntegerVector(clusters.membership.begin(), clusters.membership.end()),
                Rcpp::Named("modularity") = Rcpp::NumericVector(clusters.modularity.begin(), clusters.modularity.end())
            );
        }
    }

private:
    size_t ncells;
    std::deque<scran::BuildSNNGraph::WeightedEdge> store;
    std::string method;
    double resolution;

    scran::ClusterSNNGraphMultiLevel::Results res_multilevel;
    scran::ClusterSNNGraphLeiden::Results res_leiden;
    scran::ClusterSNNGraphWalktrap::Results res_walktrap;
};

#endif

#include "Rcpp.h"
#include "scran/clustering/ClusterSNNGraph.hpp"
#include "knncolle.h"
#ifdef _OPENMP
#include "omp.h"
#endif

struct SNNGraph {
    SNNGraph (size_t n, std::deque<scran::BuildSNNGraph::WeightedEdge> s) : ncells(n), store(std::move(s)) {}
    size_t ncells;
    std::deque<scran::BuildSNNGraph::WeightedEdge> store;
};

//[[Rcpp::export(rng=false)]]
SEXP build_graph(SEXP nnptr, int k, int nthreads) {
    KnncollePtr nns(nnptr);
    scran::BuildSNNGraph builder;
    builder.set_neighbors(k);

#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    return Rcpp::XPtr<SNNGraph>(new SNNGraph(nns->nobs(), builder.run(nns.get())));
}

//[[Rcpp::export(rng=false)]]
SEXP cluster_multilevel(SEXP ptr, double res) {
    Rcpp::XPtr<SNNGraph> x(ptr);
    scran::ClusterSNNGraph runner;
    auto clusters = runner.run_multilevel(x->ncells, x->store, res);

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


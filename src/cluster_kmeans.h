#ifndef CLUSTER_KMEANS_H
#define CLUSTER_KMEANS_H
#include "scran/clustering/ClusterKmeans.hpp"
#include <deque>
#include "Rcpp.h"
#include <string>

struct Kmeans {
    Kmeans(Rcpp::NumericMatrix data, int nc) : 
        ndim(data.nrow()), nobs(data.ncol()), ptr(data.begin()),
        nclusters(nc), centers(ndim, nclusters), clusters(data.ncol()),
        center_ptr(centers.begin()), cluster_ptr(clusters.begin()),
        withinss(nclusters), ss_ptr(withinss.begin())
    {}

    void run() {
        scran::ClusterKmeans runner;
        auto out = runner.run(ndim, nobs, ptr, nclusters, center_ptr, cluster_ptr);
        std::copy(out.withinss.begin(), out.withinss.end(), ss_ptr);
        iterations = out.iterations;
        return;
    }

    Rcpp::List yield() {
        return Rcpp::List::create(
            Rcpp::Named("clusters") = clusters, 
            Rcpp::Named("centers") = centers,
            Rcpp::Named("iterations") = Rcpp::IntegerVector::create(iterations),
            Rcpp::Named("withinss") = Rcpp::NumericVector(withinss.begin(), withinss.end())
        );
    }

private:
    int ndim;
    size_t nobs;
    const double * ptr;

    int nclusters;
    Rcpp::NumericMatrix centers;
    Rcpp::IntegerVector clusters;
    double* center_ptr;
    int* cluster_ptr;

    int iterations;
    Rcpp::NumericVector withinss;
    double* ss_ptr;
};

#endif

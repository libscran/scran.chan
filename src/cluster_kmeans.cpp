#include "Rcpp.h"
#include "scran/scran.hpp"
#include "kmeans/InitializePCAPartition.hpp"
#ifdef _OPENMP
#include "omp.h"
#endif

#include <memory>
#include <string>
#include <stdexcept>

//[[Rcpp::export(rng=false)]]
SEXP cluster_kmeans(Rcpp::NumericMatrix data, int nclusters, std::string init_method, int seed, int nthreads) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif

    int ndim = data.nrow();
    size_t nobs = data.ncol();
    auto ptr = static_cast<const double*>(data.begin());

    Rcpp::NumericMatrix centers(ndim, nclusters);
    Rcpp::IntegerVector clusters(nobs);
    auto center_ptr = static_cast<double*>(centers.begin());
    auto cluster_ptr = static_cast<int*>(clusters.begin());

    Rcpp::NumericVector withinss(nclusters);
    auto ss_ptr = static_cast<double*>(withinss.begin());

    std::shared_ptr<kmeans::Initialize<> > iptr;
    if (init_method == "random") {
        iptr.reset(new kmeans::InitializeRandom);
    } else if (init_method == "kmeans++") {
        iptr.reset(new kmeans::InitializeKmeansPP);
    } else if (init_method == "pca-part") {
        iptr.reset(new kmeans::InitializePCAPartition);
    } else {
        throw std::runtime_error("unknown init_method '" + init_method + "'");
    }

    scran::ClusterKmeans runner;
    auto out = runner.run(ndim, nobs, ptr, nclusters, center_ptr, cluster_ptr, iptr.get());
    std::copy(out.withinss.begin(), out.withinss.end(), ss_ptr);

    return Rcpp::List::create(
        Rcpp::Named("clusters") = clusters, 
        Rcpp::Named("centers") = centers,
        Rcpp::Named("iterations") = Rcpp::IntegerVector::create(out.iterations),
        Rcpp::Named("withinss") = withinss
    );
}

# Basic testing of the k-means clustering.
# library(testthat); library(scran.chan); source("test-clusterKmeans.R")

x <- t(as.matrix(iris[,1:4]))

test_that("clusterKmeans works in basic mode", {
    clustering <- clusterKmeans.chan(x)
    again <- clusterKmeans.chan(x)
    expect_identical(clustering, again)

    expect_identical(length(clustering$clusters), ncol(x))
    expect_identical(length(unique(clustering$clusters)), 10L)
    expect_identical(ncol(clustering$centers), 10L)
    expect_identical(length(clustering$withinss), 10L)

    alt <- clusterKmeans.chan(x, k=5)
    expect_identical(length(unique(alt$clusters)), 5L)
    expect_identical(ncol(alt$centers), 5L)
})

test_that("clusterKmeans works in sweep mode", {
    clustering <- clusterKmeans.chan(x, k=c(3,5), init.method=c("pca-part", "kmeans++"), seed=c(2,6))
    expect_true(nrow(clustering$parameters) > 1L)

    all.valid <- vapply(clustering$results, function(x) !is.null(x$clusters), TRUE)
    expect_true(all(all.valid))

    # Works in parallel.
    parallel <- clusterKmeans.chan(x, k=c(3,5), init.method=c("pca-part", "kmeans++"), seed=c(2,6), num.threads=2)
    expect_identical(clustering, parallel)
})

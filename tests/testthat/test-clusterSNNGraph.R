# Basic testing of the SNN clustering.
# library(testthat); library(scran.chan); source("test-clusterSNNGraph.R")

x <- t(as.matrix(iris[,1:4]))

test_that("clusterSNNGraph works in basic mode", {
    clustering <- clusterSNNGraph.chan(x)
    again <- clusterSNNGraph.chan(x)
    expect_identical(clustering, again)
    expect_identical(length(clustering$membership), ncol(x))

    # Trying some other parameters.
    alt <- clusterSNNGraph.chan(x, weight.scheme="number", method="walktrap")
    expect_identical(length(alt$membership), ncol(x))
    expect_false(identical(alt$membership, clustering$membership))

    alt2 <- clusterSNNGraph.chan(x, weight.scheme="jaccard", method="leiden")
    expect_identical(length(alt$membership), ncol(x))
})

test_that("clusterSNNGraph works in sweep mode", {
    clustering <- clusterSNNGraph.chan(x, num.neighbors=c(10, 20), resolution=c(0.5, 1.5), seed=c(10, 20))
    expect_true(nrow(clustering$parameters) > 1L)

    all.valid <- vapply(clustering$results, function(x) !is.null(x$membership), TRUE)
    expect_true(all(all.valid))

    # Works in parallel.
    parallel <- clusterSNNGraph.chan(x, num.neighbors=c(10, 20), resolution=c(0.5, 1.5), seed=c(10, 20), num.threads=2)
    expect_identical(clustering, parallel)
})

test_that("clusterSNNGraph works after downsampling", {
    clustering <- clusterSNNGraph.chan(x, downsample=5)
    expect_identical(length(clustering$membership), ncol(x))
    expect_true(!anyNA(clustering$membership))
})

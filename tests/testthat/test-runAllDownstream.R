# Tests the runAllDownstream function.
# library(testthat); library(scran.chan); source("test-runAllDownstream.R")

x <- t(as.matrix(iris[,1:4]))

test_that("runAllDownstream works more or less as expected", {
    res <- runAllDownstream(x)
    expect_identical(nrow(res$tsne), ncol(x))
    expect_identical(nrow(res$umap), ncol(x))
    expect_identical(length(res$cluster.snn$membership), ncol(x))

    # Same results in parallel
    alt <- runAllDownstream(x, num.threads=3)
    expect_identical(res, alt)
})

test_that("runAllDownstream works correctly with sweeps", {
    res <- runAllDownstream(x, 
        do.cluster.kmeans=TRUE, 
        cluster.kmeans.args=list(k=c(5,7,10)), 
        cluster.snn.args=list(method=c("leiden", "multilevel"))
    )

    expect_identical(length(res$cluster.kmeans$results[[1]]$clusters), ncol(x))

    # Still dropped for non-sweep parameters.
    expect_identical(nrow(res$tsne), ncol(x))
    expect_identical(nrow(res$umap), ncol(x))

    expect_true(nrow(res$cluster.snn$parameters) > 1L)
    expect_identical(length(res$cluster.snn$results[[1]]$membership), ncol(x))

    expect_true(nrow(res$cluster.kmeans$parameters) > 1L)
    expect_identical(length(res$cluster.kmeans$results[[1]]$clusters), ncol(x))

    # Same results in parallel
    alt <- runAllDownstream(x, 
        do.cluster.kmeans=TRUE, 
        cluster.kmeans.args=list(k=c(5,7,10)), 
        cluster.snn.args=list(method=c("leiden", "multilevel")),
        num.threads=3
    )
    expect_identical(res, alt)
})

test_that("runAllDownstream works correctly with downsampling", {
    res <- runAllDownstream(x, downsample=5)
    expect_identical(nrow(res$tsne), ncol(x))
    expect_identical(nrow(res$umap), ncol(x))
    expect_identical(length(res$cluster.snn$membership), ncol(x))

    res <- runAllDownstream(x, do.cluster.kmeans=TRUE, do.umap=FALSE, do.tsne=FALSE, do.cluster.snn=FALSE, downsample=5)
    expect_identical(length(res$cluster.kmeans$clusters), ncol(x))
})

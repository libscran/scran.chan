# Basic testing of the UMAP.
# library(testthat); library(scran.chan); source("test-runUMAP.R")

x <- t(as.matrix(iris[,1:4]))

test_that("runUMAP works in basic mode", {
    umap <- runUMAP.chan(x)
    again <- runUMAP.chan(x)
    expect_identical(umap, again)

    expect_identical(nrow(umap), ncol(x))
    expect_identical(ncol(umap), 2L)

    alt <- runUMAP.chan(x, num.neighbors=20)
    expect_identical(dim(alt), dim(umap))
    expect_false(identical(alt, umap))
})

test_that("runUMAP works in sweep mode", {
    umap <- runUMAP.chan(x, num.neighbors=c(10, 15), min.dist=c(0.01, 0.05), seed=c(2,6))
    expect_true(nrow(umap$parameters) > 1L)

    all.valid <- vapply(umap$results, function(y) nrow(y) == ncol(x), TRUE)
    expect_true(all(all.valid))

    # Works in parallel.
    parallel <- runUMAP.chan(x, num.neighbors=c(10, 15), min.dist=c(0.01, 0.05), seed=c(2,6), num.threads=2)
    expect_identical(umap, parallel)
})

test_that("runUMAP works with downsampling", {
    umap <- runUMAP.chan(x, downsample=5)
    expect_identical(nrow(umap), ncol(x))
})

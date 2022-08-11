# library(testthat); library(scran.chan); source("test-downsampleByNeighbors.R")

x <- t(as.matrix(iris[,1:4]))

test_that("downsampleByNeighbors works as expected", {
    keep <- downsampleByNeighbors.chan(x, k=5)
    expect_true(length(keep) < ncol(x))
    expect_true(min(keep) > 0)
    expect_true(max(keep) <= ncol(x))

    keep2 <- downsampleByNeighbors.chan(x, k=10)
    expect_true(length(keep2) < length(keep))
    expect_true(min(keep2) > 0)
    expect_true(max(keep2) <= ncol(x))
})


# library(testthat); library(scran.chan); source("test-downsampleByNeighbors.R")

x <- t(as.matrix(iris[,1:4]))

test_that("downsampleByNeighbors works as expected", {
    chosen <- downsampleByNeighbors.chan(x, k=5)
    keep <- chosen$chosen
    expect_true(length(keep) < ncol(x))
    expect_true(min(keep) > 0)
    expect_true(max(keep) <= ncol(x))
    expect_true(all(keep %in% chosen$assigned))
    expect_true(all(chosen$assigned %in% keep))

    chosen2 <- downsampleByNeighbors.chan(x, k=10)
    keep2 <- chosen2$chosen
    expect_true(length(keep2) < length(keep))
    expect_true(min(keep2) > 0)
    expect_true(max(keep2) <= ncol(x))
    expect_true(all(keep2 %in% chosen2$assigned))
    expect_true(all(chosen2$assigned %in% keep2))
})


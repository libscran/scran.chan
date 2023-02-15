# Tests the groupedSizeFactors functionality.
# library(testthat); library(scran.chan); source("test-groupedSizeFactors.R")

library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
y <- initializeSparseMatrix(x)

test_that("groupedSizeFactors works as expected", {
    sf <- groupedSizeFactors(y, group=NULL)
    expect_identical(length(sf), ncol(x))
    expect_true(all(sf > 0))

    # Providing my own clustering.
    clust <- LETTERS[sample(5, ncol(x), replace=TRUE)]
    sf <- groupedSizeFactors(y, clust, reference="A")
    expect_identical(length(sf), ncol(x))
    expect_true(all(sf > 0))
})

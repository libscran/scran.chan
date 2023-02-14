# This tests the scoreFeatureSet machinery.
# library(testthat); library(scran.chan); source("test-scoreFeatureSet.R")

library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.5) * 100))
y <- initializeSparseMatrix(x)
normed <- logNormCounts.chan(y)

test_that("scoreFeatureSet works more or less as expected", {
    chosen <- c(1,3,5,10,20,100)
    scores <- scoreFeatureSet.chan(normed, chosen)
    expect_identical(length(chosen), length(scores$weights))
    expect_identical(length(scores$scores), ncol(x))
})

test_that("scoreFeatureSet works with blocking", {
    chosen <- c(2, 4, 6, 8, 10, 12, 14, 16, 18)
    block <- sample(3, ncol(x), replace=TRUE)
    scores <- scoreFeatureSet.chan(normed, chosen, batch=block)
    expect_identical(length(chosen), length(scores$weights))
    expect_identical(length(scores$scores), ncol(x))
})


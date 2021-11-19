# Tests for correct centering of the size factors.
# library(testthat); library(scran.chan); source("test-logNormCounts.R")

set.seed(100)
test_that("logNormCounts works as expected in basic usage", {
    library(Matrix)
    x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
    y <- initializeSparseMatrix(x)
    normed <- logNormCounts.chan(y)

    expected <- colSums(x)
    expected <- expected / mean(expected)
    expect_equal(expected, normed$size.factors)

    for (i in seq_len(ncol(x))) {
        ref <- log2(x[,i] / expected[i] + 1)
        expect_equal(ref, scran.chan:::tatami_columns(normed$pointer, i-1, 0, nrow(x))[,1])
    }

    # Listens to the size factors.
    sf <- runif(ncol(x))
    normed <- logNormCounts.chan(y, size.factors=sf)
    expected <- sf / mean(sf)
    expect_equal(expected, normed$size.factors)
})

test_that("logNormCounts works as expected when dealing with blocks", {
    x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
    y <- initializeSparseMatrix(x)

    block <- sample(3, ncol(x), replace=TRUE)
    normed <- logNormCounts.chan(y, batch=block)

    expected <- colSums(x)
    scale <- min(vapply(split(expected, block), mean, 0))
    expected <- expected / scale
    expect_equal(expected, normed$size.factors)

    for (i in seq_len(ncol(x))) {
        ref <- log2(x[,i] / expected[i] + 1)
        expect_equal(ref, scran.chan:::tatami_columns(normed$pointer, i-1, 0, nrow(x))[,1])
    }

    # Trying again with the other batch mode.
    normed <- logNormCounts.chan(y, batch=block, batch.mode="perblock")
    expected <- colSums(x)
    scale <- vapply(split(expected, factor(block, seq_len(3))), mean, 0)
    expected <- expected / unname(scale)[block]
    expect_equal(expected, normed$size.factors)
})

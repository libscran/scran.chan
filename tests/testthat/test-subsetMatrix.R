# This contains some brief tests for subsetMatrix
# library(testthat); library(scran.chan); source("test-subsetMatrix.R")

library(Matrix)
x <- round(abs(rsparsematrix(100, 10, 0.1) * 100))
rownames(x) <- sprintf("GENE_%i", seq_len(100))
colnames(x) <- LETTERS[1:10]

compare <- function(ref, ptr) {
    expect_identical(dim(ref), scran.chan:::tatami_dim(ptr$pointer))
    for (i in seq_len(ncol(ref))) {
        expect_identical(unname(ref[,i]), scran.chan:::tatami_columns(ptr$pointer, i-1, 0, nrow(ref))[,1])
    }
}

test_that("subsetMatrix works as expected", {
    y <- initializeSparseMatrix(x)

    sub <- subsetMatrix(y, NULL, NULL) # no-op.
    compare(x, sub)

    ref <- x[1:5, c("A", "C", "E")]
    sub <- subsetMatrix(y, 1:5, c("A", "C", "E"))
    compare(ref, sub)

    ref <- x[10:5,]
    sub <- subsetMatrix(y, 10:5, NULL)
    compare(ref, sub)

    keep <- rep(c(TRUE, FALSE), length.out=10)
    ref <- x[,keep]
    sub <- subsetMatrix(y, NULL, keep)
    compare(ref, sub)

    keep <- sample(rownames(x), 100)
    ref <- x[keep, 5:1]
    sub <- subsetMatrix(y, keep, 5:1)
    compare(ref, sub)
})

test_that("subsetMatrix slices the dimnames", {
    y <- initializeSparseMatrix(x)

    ref <- x[1:5, c("A", "C", "E")]
    sub <- subsetMatrix(y, 1:5, c("A", "C", "E"))
    compare(ref, sub)

    expect_identical(sub$rownames, rownames(ref))
    expect_identical(sub$colnames, colnames(ref))
})

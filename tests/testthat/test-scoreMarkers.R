# This contains some brief tests for scoreMarkers
# library(testthat); library(scran.chan); source("test-scoreMarkers.R")

test_that("scoreMarkers works as expected for simple cases", {
    library(Matrix)
    x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
    y <- initializeSparseMatrix(x)

    g <- sample(4, ncol(x), replace=TRUE)
    markers <- scoreMarkers.chan(y, g)

    df1 <- markers$statistics[["1"]]
    expect_equal(df1$mean, rowMeans(x[as.integer(rownames(df1)),g==1]))

    df4 <- markers$statistics[["4"]]
    expect_equal(df4$detected, rowMeans(x[as.integer(rownames(df4)),g==4] > 0))
})

test_that("scoreMarkers works as expected for the full set of statistics", {
    library(Matrix)
    x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
    y <- initializeSparseMatrix(x)

    g <- sample(4, ncol(x), replace=TRUE)
    markers <- scoreMarkers.chan(y, g)
    full <- scoreMarkers.chan(y, g, simple.means.only=FALSE)

    expect_identical(markers$statistics[["2"]]$logFC, full$statistics[["2"]]$logFC.mean) 
    expect_identical(markers$statistics[["3"]]$logFC, full$statistics[["3"]]$logFC.mean) 
})

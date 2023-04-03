# Tests the calculation and filtering of RNA-based QC metrics
# library(testthat); library(scran.chan); source("test-perCellCrisprQcMetrics.R")

library(Matrix)
x <- round(abs(rsparsematrix(25, 100, 0.1) * 100))
y <- initializeSparseMatrix(x)
z <- as.matrix(x)

test_that("perCellCrisprQcMetrics works as expected", { 
    qc <- perCellCrisprQcMetrics.chan(y)
    expect_identical(qc$sum, Matrix::colSums(x))
    expect_identical(qc$detected, Matrix::colSums(x > 0))

    top <- max.col(t(z))
    expect_equal(qc$max.proportion, z[cbind(top, seq_along(top))] / qc$sum)

    # Can't compare top and top2 directly in case of ties.
    top2 <- qc$max.index
    expect_equal(qc$max.proportion, z[cbind(top2, seq_along(top2))] / qc$sum)
})

test_that("suggestCrisprQcFilters works as expected", { 
    qc <- perCellCrisprQcMetrics.chan(y)
    filtered <- suggestCrisprQcFilters.chan(qc$sums, qc$max.proportion)

    expected <- qc$max.proportion * qc$sums >= filtered$thresholds$max.count
    expected[is.na(expected)] <- TRUE
    expect_identical(!expected, filtered$filter)
})

test_that("suggestCrisprQcFilters works as expected with blocking", { 
    qc <- perCellCrisprQcMetrics.chan(y)
    block <- sample(3, ncol(x), replace=TRUE)
    filtered <- suggestCrisprQcFilters.chan(qc$sums, qc$max.proportion, batch=block)

    expected <- qc$max.proportion * qc$sums >= unname(filtered$thresholds$max.count[block])
    expected[is.na(expected)] <- TRUE
    expect_identical(!expected, filtered$filter)
})


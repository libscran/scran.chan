# Tests the calculation and filtering of RNA-based QC metrics
# library(testthat); library(scran.chan); source("test-perCellAdtQcMetrics.R")

library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
y <- initializeSparseMatrix(x)

test_that("perCellAdtQcMetrics works as expected", { 
    sub <- list(IgG=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- perCellAdtQcMetrics.chan(y, sub)
    expect_identical(qc$sum, Matrix::colSums(x))
    expect_identical(qc$detected, Matrix::colSums(x > 0))
    expect_identical(qc$subsets$IgG, Matrix::colSums(x[sub$IgG,]))
})

test_that("suggestAdtQcFilters works as expected", { 
    sub <- list(IgG=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- perCellAdtQcMetrics.chan(y, sub)
    filtered <- suggestAdtQcFilters.chan(qc$detected, qc$subsets)

    # Check the thresholds.
    ldet <- log(qc$detected)
    expect_equal(filtered$thresholds$detected, exp(median(ldet) - 3 * mad(ldet)))
    lsub <- log(qc$subsets[[1]])
    expect_equal(filtered$thresholds$subsets[[1]], exp(median(lsub) + 3 * mad(lsub)))

    # Check the filter.
    expected <- qc$detected >= filtered$thresholds$detected & qc$subsets[[1]] <= filtered$thresholds$subsets[[1]]
    expect_identical(!expected, filtered$filter)
})

test_that("suggestAdtQcFilters works as expected with blocking", { 
    sub <- list(IgG=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- perCellAdtQcMetrics.chan(y, sub)
    block <- sample(3, ncol(x), replace=TRUE)
    filtered <- suggestAdtQcFilters.chan(qc$detected, qc$subsets, batch=block)

    # Check the thresholds.
    for (b in 1:3) {
        keep <- block == b
        ldet <- log(qc$detected[keep])
        expect_equal(filtered$thresholds$detected[[b]], exp(median(ldet) - 3 * mad(ldet)))
        lsub <- log(qc$subsets[[1]][keep])
        expect_equal(filtered$thresholds$subsets[[1]][[b]], exp(median(lsub) + 3 * mad(lsub)))
    }

    # Check the filter.
    expected <- qc$detected >= filtered$thresholds$detected[block] & qc$subsets[[1]] <= filtered$thresholds$subsets[[1]][block]
    expect_identical(!unname(expected), filtered$filter)
})


# Tests the calculation and filtering of RNA-based QC metrics
# library(testthat); library(scran.chan); source("test-perCellRnaQcMetrics.R")

library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
y <- initializeSparseMatrix(x)

test_that("perCellRnaQcMetrics works as expected", { 
    sub <- list(Mito=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- perCellRnaQcMetrics.chan(y, sub)
    expect_identical(qc$sum, Matrix::colSums(x))
    expect_identical(qc$detected, Matrix::colSums(x > 0))
    expect_identical(qc$subsets$Mito, Matrix::colSums(x[sub$Mito,]) / qc$sum)
})

test_that("suggestRnaQcFilters works as expected", { 
    sub <- list(Mito=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- perCellRnaQcMetrics.chan(y, sub)
    filtered <- suggestRnaQcFilters.chan(qc$sums, qc$detected, qc$subsets)

    # Check the thresholds.
    lsums <- log(qc$sums)
    expect_equal(filtered$thresholds$sums, exp(median(lsums) - 3 * mad(lsums)))
    ldet <- log(qc$detected)
    expect_equal(filtered$thresholds$detected, exp(median(ldet) - 3 * mad(ldet)))

    sub <- qc$subsets[[1]]
    expect_equal(filtered$thresholds$subsets[[1]], median(sub) + 3 * mad(sub))
    expect_equal(names(filtered$thresholds$subsets), "Mito")

    # Check the filter.
    expected <- qc$sums >= filtered$thresholds$sums & qc$detected >= filtered$thresholds$detected & qc$subsets[[1]] <= filtered$thresholds$subsets[[1]]
    expect_identical(!expected, filtered$filter)

    # Also check back-compatibility.
    old <- perCellQCFilters.chan(qc$sums, qc$detected, qc$subsets)
    expect_identical(filtered$filter, old$filters$sums | old$filters$detected | old$filters$subsets[[1]])
})

test_that("suggestRnaQcFilters works as expected with blocking", { 
    sub <- list(Mito=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- perCellRnaQcMetrics.chan(y, sub)
    block <- sample(3, ncol(x), replace=TRUE)
    filtered <- suggestRnaQcFilters.chan(qc$sums, qc$detected, qc$subsets, batch=block)

    # Check the thresholds.
    for (b in 1:3) {
        keep <- block == b
        lsums <- log(qc$sums[keep])
        expect_equal(filtered$thresholds$sums[[b]], exp(median(lsums) - 3 * mad(lsums)))
        ldet <- log(qc$detected[keep])
        expect_equal(filtered$thresholds$detected[[b]], exp(median(ldet) - 3 * mad(ldet)))
        sub <- qc$subsets[[1]][keep]
        expect_equal(filtered$thresholds$subsets[[1]][[b]], median(sub) + 3 * mad(sub))
    }

    # Check the filter.
    expected <- qc$sums >= filtered$thresholds$sums[block] & qc$detected >= filtered$thresholds$detected[block] & qc$subsets[[1]] <= filtered$thresholds$subsets[[1]][block]
    expect_identical(!unname(expected), filtered$filter)

    # Also check back-compatibility.
    old <- perCellQCFilters.chan(qc$sums, qc$detected, qc$subsets, batch=block)
    expect_identical(filtered$filter, old$filters$sums | old$filters$detected | old$filters$subsets[[1]])
})


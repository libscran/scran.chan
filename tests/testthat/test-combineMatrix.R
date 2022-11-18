# This contains some brief tests for comsetMatrix
# library(testthat); library(scran.chan); source("setup.R"); source("test-combineMatrix.R")

library(Matrix)
x1 <- round(abs(rsparsematrix(100, 10, 0.1) * 100))
rownames(x1) <- sprintf("GENE_%i", seq_len(100))
colnames(x1) <- LETTERS[1:10]

test_that("rbindMatrix works as expected", {
    y1 <- initializeSparseMatrix(x1)
    com <- rbindMatrix(y1)
    compare(x1, com)
    expect_identical(com$rownames, rownames(x1))
    expect_identical(com$colnames, colnames(x1))

    # Adding a matrix.
    x2 <- round(abs(rsparsematrix(50, 10, 0.1) * 100))
    rownames(x2) <- sprintf("gene_%i", seq_len(50))
    colnames(x2) <- letters[1:10]
    y2 <- initializeSparseMatrix(x2)

    com <- rbindMatrix(y1, y2)
    ref <- rbind(x1, x2)
    compare(ref, com)
    expect_identical(com$rownames, c(rownames(x1), rownames(x2)))
    expect_identical(com$colnames, colnames(x1))

    # Adding another matrix.
    x3 <- round(abs(rsparsematrix(200, 10, 0.1) * 100))
    rownames(x3) <- sprintf("feat_%i", seq_len(200))
    colnames(x3) <- letters[26:17]
    y3 <- initializeSparseMatrix(x3)

    com <- rbindMatrix(y1, y2, y3)
    ref <- rbind(x1, x2, x3)
    compare(ref, com)
    expect_identical(com$rownames, c(rownames(x1), rownames(x2), rownames(x3)))
    expect_identical(com$colnames, colnames(x1))

    # Removing some of the row names.
    y1$rownames <- NULL
    com <- rbindMatrix(y1, y2, y3)
    expect_identical(com$rownames, c(character(nrow(x1)), rownames(x2), rownames(x3)))

    # Removing all of the row names.
    y2$rownames <- NULL
    y3$rownames <- NULL
    com <- rbindMatrix(y1, y2, y3)
    expect_null(com$rownames)
})

test_that("cbindMatrix works as expected", {
    y1 <- initializeSparseMatrix(x1)
    com <- cbindMatrix(y1)
    compare(x1, com)
    expect_identical(com$rownames, rownames(x1))
    expect_identical(com$colnames, colnames(x1))

    # Adding a matrix.
    x2 <- round(abs(rsparsematrix(100, 20, 0.1) * 100))
    rownames(x2) <- sprintf("gene_%i", seq_len(100))
    colnames(x2) <- letters[1:20]
    y2 <- initializeSparseMatrix(x2)

    com <- cbindMatrix(y1, y2)
    ref <- cbind(x1, x2)
    compare(ref, com)
    expect_identical(com$rownames, rownames(x1))
    expect_identical(com$colnames, c(colnames(x1), colnames(x2)))

    # Adding another matrix.
    x3 <- round(abs(rsparsematrix(100, 25, 0.1) * 100))
    rownames(x3) <- sprintf("feat_%i", seq_len(100))
    colnames(x3) <- letters[26:2]
    y3 <- initializeSparseMatrix(x3)

    com <- cbindMatrix(y1, y2, y3)
    ref <- cbind(x1, x2, x3)
    compare(ref, com)
    expect_identical(com$rownames, rownames(x1))
    expect_identical(com$colnames, c(colnames(x1), colnames(x2), colnames(x3)))

    # Removing some of the row names.
    y1$colnames <- NULL
    com <- cbindMatrix(y1, y2, y3)
    expect_identical(com$colnames, c(character(ncol(x1)), colnames(x2), colnames(x3)))

    # Removing all of the col names.
    y2$colnames <- NULL
    y3$colnames <- NULL
    com <- cbindMatrix(y1, y2, y3)
    expect_null(com$colnames)
})

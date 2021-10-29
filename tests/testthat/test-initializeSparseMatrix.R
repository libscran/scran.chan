# This checks the initialization procedure.
# library(testthat); library(scran.chan); source("test-initializeSparseMatrix.R")

set.seed(1000)
x <- Matrix::rsparsematrix(1000, 100, 0.1)
y <- round(abs(x)*10)

am_i_ok <- function(ref, ptr) {
    expect_identical(dim(ref), scran.chan:::tatami_dim(ptr$pointer))
    for (i in seq_len(ncol(ref))) {
        expect_identical(ref[,i], scran.chan:::tatami_columns(ptr$pointer, i-1, 0, nrow(ref))[,1])
    }
}

test_that("initialization works correctly with a dgCMatrix", {
    ptr <- initializeSparseMatrix(y)
    am_i_ok(y, ptr)

    # Force a copy to integer.
    ptr2 <- initializeSparseMatrix(y, no.sparse.copy=FALSE)
    am_i_ok(y, ptr2)
})

test_that("initialization works correctly with a dgRMatrix", {
    z <- new("dgRMatrix", x=y@x, j=y@i, p=y@p, Dim=rev(y@Dim))
    ptr <- initializeSparseMatrix(z)
    am_i_ok(z, ptr)
})

test_that("initialization works correctly with a DelayedMatrix", {
    z <- DelayedArray::DelayedArray(y)

    # Column-based.
    stuff2 <- initializeSparseMatrix(z)
    am_i_ok(z, stuff2)

    # Row-major initialization:
    stuff2 <- initializeSparseMatrix(z, by.column=FALSE) 
    am_i_ok(z, stuff2)
})

library(rhdf5)
tmp <- tempfile(fileext=".h5")
h5createFile(tmp)
h5createGroup(tmp, "YAY")
h5write(y@i, tmp, "YAY/indices")
h5write(y@x, tmp, "YAY/data")
h5write(y@p, tmp, "YAY/indptr")
h5write(dim(y), tmp, "YAY/shape")

library(HDF5Array)
A <- TENxMatrix(filepath=tmp, group="YAY")

test_that("initialization works correctly with H5SparseMatrices", {
    z <- initializeSparseMatrix(A)
    am_i_ok(A, z)
})

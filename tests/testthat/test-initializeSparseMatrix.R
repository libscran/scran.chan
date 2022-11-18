# This checks the initialization procedure.
# library(testthat); library(scran.chan); source("setup.R"); source("test-initializeSparseMatrix.R")

set.seed(1000)
x <- Matrix::rsparsematrix(1000, 100, 0.1)
y <- round(abs(x)*10)

test_that("initialization works correctly with a dgCMatrix", {
    ptr <- initializeSparseMatrix(y)
    am_i_ok(y, ptr)

    # Force a copy to integer.
    ptr2a <- initializeSparseMatrix(y, no.sparse.copy=FALSE)
    am_i_ok(y, ptr2a)

    # Force copy to a larger integer type for 'x'.
    y2 <- y * 1e5
    ptr2b <- initializeSparseMatrix(y2, no.sparse.copy=FALSE)
    am_i_ok(y2, ptr2b)

    # Check for correct handling of doubles.
    y3 <- y * 1.5
    ptr3a <- initializeSparseMatrix(y3, no.sparse.copy=FALSE)
    am_i_ok(floor(y3), ptr3a)

    ptr3b <- initializeSparseMatrix(y3, no.sparse.copy=TRUE)
    am_i_ok(y3, ptr3b)

    ptr3c <- initializeSparseMatrix(y3, no.sparse.copy=FALSE, force.integer=FALSE)
    am_i_ok(y3, ptr3c)

    # Force the use of a larger integer type for 'i'.
    y4 <- round(abs(Matrix::rsparsematrix(1e5, 10, 0.01) * 10))
    ptr4 <- initializeSparseMatrix(y4, no.sparse.copy=FALSE)
    am_i_ok(y4, ptr4)
})

test_that("initialization works correctly with a dgRMatrix", {
    z <- new("dgRMatrix", x=y@x, j=y@i, p=y@p, Dim=rev(y@Dim))
    ptr <- initializeSparseMatrix(z)
    am_i_ok(z, ptr)

    # Force a copy.
    ptr2 <- initializeSparseMatrix(z, no.sparse.copy=FALSE)
    am_i_ok(z, ptr2)
})

test_that("initialization works correctly with a DelayedMatrix", {
    # Multiplication and reduce block size to force non-trivial block processing. 
    DelayedArray::setAutoBlockSize(200)
    z <- DelayedArray::DelayedArray(y) * 2 

    # Column-based.
    stuff1a <- initializeSparseMatrix(z)
    am_i_ok(z, stuff1a)

    # Row-major initialization:
    stuff1b <- initializeSparseMatrix(z, by.column=FALSE) 
    am_i_ok(z, stuff1b)

    # Respects the double type.
    z2 <- DelayedArray::DelayedArray(y) * 1.5
    stuff2a <- initializeSparseMatrix(z2)
    am_i_ok(floor(z2), stuff2a)

    stuff2b <- initializeSparseMatrix(z2, force.integer=FALSE)
    am_i_ok(z2, stuff2b)

    # Force the use of a larger type for 'i'.
    y3 <- DelayedArray::DelayedArray(round(abs(Matrix::rsparsematrix(7e4, 5, 0.01) * 10))) * 2
    ptr3 <- initializeSparseMatrix(y3)
    am_i_ok(y3, ptr3)
})

library(HDF5Array)
library(rhdf5)
dump_hdf5 <- function(y, xtype="H5T_NATIVE_DOUBLE", itype="H5T_NATIVE_INT") {
    tmp <- tempfile(fileext=".h5")
    h5createFile(tmp)
    h5createGroup(tmp, "YAY")
    h5write(y@p, tmp, "YAY/indptr")
    h5write(dim(y), tmp, "YAY/shape")

    h5createDataset(tmp, "YAY/indices", dims=length(y@i), H5type=itype)
    h5write(y@i, tmp, "YAY/indices")

    h5createDataset(tmp, "YAY/data", dims=length(y@x), H5type=xtype)
    h5write(y@x, tmp, "YAY/data")

    TENxMatrix(filepath=tmp, group="YAY")
}

test_that("initialization works correctly with H5SparseMatrices", {
    A <- dump_hdf5(y)
    z <- initializeSparseMatrix(A)
    am_i_ok(A, z)

    # Auto-upgrates to a larger integer size.
    A <- dump_hdf5(y * 1e5)
    z <- initializeSparseMatrix(A)
    am_i_ok(A, z)

    # Respects doubles.
    A <- dump_hdf5(y * 1.5)
    z <- initializeSparseMatrix(A)
    am_i_ok(floor(A), z)

    z <- initializeSparseMatrix(A, force.integer=FALSE)
    am_i_ok(A, z)

    # Respects integer data types.
    A <- dump_hdf5(y, xtype="H5T_NATIVE_UINT16")
    z <- initializeSparseMatrix(A)
    am_i_ok(A * 1.0, z)

    # Force the use of a larger integer type for 'i'.
    A <- dump_hdf5(round(abs(Matrix::rsparsematrix(1e5, 10, 0.01) * 10)))
    z <- initializeSparseMatrix(A, no.sparse.copy=FALSE)
    am_i_ok(A, z)
})

test_that("initialization preserves row and column names", {
    rownames(y) <- seq_len(nrow(y))
    colnames(y) <- paste0("G", seq_len(ncol(y)))
    ptr <- initializeSparseMatrix(y)
    expect_identical(rownames(y), ptr$rownames)
    expect_identical(colnames(y), ptr$colnames)

    A <- dump_hdf5(y)
    z <- initializeSparseMatrix(A)
    expect_identical(rownames(A), z$rownames)
    expect_identical(colnames(A), z$colnames)
})



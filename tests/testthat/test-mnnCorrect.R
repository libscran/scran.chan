# This tests that the MNN correction runs correctly.
# library(testthat); library(scran.chan); source("test-mnnCorrect.R")

test_that("mnnCorrect works correctly in simple cases", {
    x <- t(as.matrix(iris[,1:4]))
    b <- sample(1:2, ncol(x), replace=TRUE)
    corrected <- mnnCorrect.chan(x, b)

    expect_identical(dim(corrected$corrected), dim(x))
    expect_identical(length(corrected$num.pairs), 1L)
    expect_identical(sort(corrected$merge.order), as.character(1:2))
})

test_that("mnnCorrect works with a specified order", {
    x <- t(as.matrix(iris[,1:4]))
    b <- sample(1:2, ncol(x), replace=TRUE)
    corrected <- mnnCorrect.chan(x, b, order=c("2", "1"))

    expect_identical(dim(corrected$corrected), dim(x))
    expect_identical(length(corrected$num.pairs), 1L)
    expect_identical(corrected$merge.order, as.character(2:1))

    # Actually has an effect.
    corrected2 <- mnnCorrect.chan(x, b, order=c("1", "2"))
    expect_false(identical(corrected$corrected, corrected2$corrected))
})

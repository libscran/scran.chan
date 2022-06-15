# This tests the runPCA machinery.
# library(testthat); library(scran.chan); source("test-runPCA.R")

library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
y <- initializeSparseMatrix(x)
normed <- logNormCounts.chan(y)

test_that("runPCA works more or less as expected", {
    pcs <- runPCA.chan(normed)
    expect_false(is.unsorted(-pcs$prop.variance))

    rm <- rowMeans(pcs$components)
    expect_true(all(abs(rm) < 1e-8))
    expect_identical(length(rm), 50L)

    rv <- apply(pcs$components, 1, var)
    div <- rv/pcs$prop.variance
    expect_true(sd(div) < 1e-8)

    # Reports the rotation matrix properly.
    pcs_with_rot <- runPCA.chan(normed, rotation=TRUE)
    expect_identical(nrow(pcs_with_rot$rotation), 1000L)
    expect_identical(ncol(pcs_with_rot$rotation), 50L)

    # Only using features of interest.
    subset <- 1:80
    subbed <- runPCA.chan(normed, subset=subset, rotation=TRUE)
    expect_identical(nrow(subbed$rotation), length(subset))
    expect_false(isTRUE(all.equal(subbed$components, pcs$components)))
})

test_that("runPCA works with blocking", {
    block <- sample(3, ncol(x), replace=TRUE)
    pcs <- runPCA.chan(normed, batch=block)
    expect_false(is.unsorted(-pcs$prop.variance))

    # Same logic applies, but this time the PCs are different.
    rm <- rowMeans(pcs$components)
    expect_true(all(abs(rm) < 1e-8))
    expect_identical(length(rm), 50L)

    rv <- apply(pcs$components, 1, var)
    div <- rv/pcs$prop.variance
    expect_true(sd(div) < 1e-8)

    ref <- runPCA.chan(normed)
    expect_false(isTRUE(all.equal(ref$components, pcs$components)))

    # Reports the rotation matrix properly.
    pcs_with_rot <- runPCA.chan(normed, batch=block, rotation=TRUE)
    expect_identical(nrow(pcs_with_rot$rotation), 1000L)
    expect_identical(ncol(pcs_with_rot$rotation), 50L)

    # Only using features of interest.
    subset <- rbinom(1000, 1, 0.2) == 1
    subbed <- runPCA.chan(normed, subset=subset, batch=block, rotation=TRUE)
    expect_identical(nrow(subbed$rotation), sum(subset))
    expect_false(isTRUE(all.equal(subbed$components, pcs$components)))
})

test_that("runPCA works with weighting", {
    block <- sample(3, ncol(x), replace=TRUE)
    pcs <- runPCA.chan(normed, batch=block, batch.method="weight")
    expect_false(is.unsorted(-pcs$prop.variance))

    # Make sure it's different from the other options.
    ref1 <- runPCA.chan(normed)
    expect_false(isTRUE(all.equal(ref1$components, pcs$components)))
    ref2 <- runPCA.chan(normed, batch=block)
    expect_false(isTRUE(all.equal(ref2$components, pcs$components)))

    # Reports the rotation matrix properly.
    pcs_with_rot <- runPCA.chan(normed, batch=block, rotation=TRUE)
    expect_identical(nrow(pcs_with_rot$rotation), 1000L)
    expect_identical(ncol(pcs_with_rot$rotation), 50L)

    # Only using features of interest.
    normed$rownames <- sprintf("GENE_%s", seq_len(1000))
    subset <- sample(normed$rownames, 200)
    subbed <- runPCA.chan(normed, subset=subset, batch=block, rotation=TRUE)
    expect_identical(rownames(subbed$rotation), normed$rownames[normed$rownames %in% subset])
    expect_false(isTRUE(all.equal(subbed$components, pcs$components)))
})

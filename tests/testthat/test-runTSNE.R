# Basic testing of the TSNE.
# library(testthat); library(scran.chan); source("test-runTSNE.R")

x <- t(as.matrix(iris[,1:4]))

test_that("runTSNE works in basic mode", {
    embed <- runTSNE.chan(x)
    again <- runTSNE.chan(x)
    expect_identical(embed, again)

    expect_identical(nrow(embed), ncol(x))
    expect_identical(ncol(embed), 2L)

    alt <- runTSNE.chan(x, perplexity=20)
    expect_identical(dim(alt), dim(embed))
    expect_false(identical(alt, embed))
})

test_that("runTSNE works in sweep mode", {
    embed <- runTSNE.chan(x, perplexity=c(10, 20, 30), seed=c(2,6))
    expect_true(nrow(embed$parameters) > 1L)

    all.valid <- vapply(embed$results, function(y) nrow(y) == ncol(x), TRUE)
    expect_true(all(all.valid))

    # Works in parallel.
    parallel <- runTSNE.chan(x, perplexity=c(10, 20, 30), seed=c(2,6), num.threads=2)
    expect_identical(embed, parallel)
})

test_that("runTSNE works with downsampling", {
    umap <- runTSNE.chan(x, downsample=5)
    expect_identical(nrow(umap), ncol(x))
})

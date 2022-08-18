# library(scran.chan); library(testthat); source("test-projectNeighborEmbedding.R")

set.seed(1001)
ndim <- 8
clust.size <- c(100, 300, 200)
nobs <- sum(clust.size)
x <- matrix(rnorm(ndim * nobs, mean=rep(1:3, clust.size)), nrow=ndim, ncol=nobs, byrow=TRUE)
embed <- matrix(rnorm(nobs * 2, mean=rep(1:3, clust.size), 0.1), nrow=nobs, ncol=2)

test_that("projectNeighborEmbedding works as expected", {
    y <- matrix(3:1, nrow=ndim, ncol=3, byrow=TRUE)
    full <- projectNeighborEmbedding.chan(x, embed, y)

    expect_true(all(full[1,] > 2.8 & full[1,] < 3.2))
    expect_true(all(full[2,] > 1.8 & full[2,] < 2.2))
    expect_true(all(full[3,] > 0.8 & full[3,] < 1.2))

    # Adding names and stuff.
    colnames(embed) <- c("A", "B")
    colnames(y) <- seq_len(ncol(y))
    full <- projectNeighborEmbedding.chan(x, embed, y)
    expect_identical(rownames(full), colnames(y))
    expect_identical(colnames(full), colnames(embed))
})

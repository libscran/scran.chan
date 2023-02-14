# This tests the scaleByNeighbors function.
# library(testthat); library(scran.chan); source("test-scaleByNeighbors.R")

pcs <- list(
    gene = matrix(rnorm(10000), ncol=200),
    protein = matrix(rnorm(1000, sd=3), ncol=200),
    guide = matrix(rnorm(2000, sd=5), ncol=200)
)

manual_knn_distance <- function(x, k = 20) {
    value <- numeric(ncol(x))
    for (i in seq_len(ncol(x))) {
        chosen <- x[,i]
        distances <- sqrt(colSums((x[,-i] - chosen)^2))
        value[i] <- sort(distances)[k]
    } 
    median(value)
}

test_that("scaleByNeighbors works as expected", {
    out <- scaleByNeighbors.chan(pcs, approximate=FALSE)
    expect_identical(names(out$scaling), names(pcs))

    references <- lapply(pcs, manual_knn_distance)
    expect_equal(out$scaling, references[[1]]/unlist(references))

    expect_identical(ncol(out$combined), ncol(pcs[[1]]))
    expect_identical(nrow(out$combined), sum(vapply(pcs, nrow, 0L)))
})

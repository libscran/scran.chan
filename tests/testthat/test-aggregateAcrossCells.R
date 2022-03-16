# This contains some brief tests for aggregateAcrossCells.
# library(testthat); library(scran.chan); source("test-aggregateAcrossCells.R")

library(Matrix)

test_that("aggregateAcrossCells works as expected for simple cases", {
    x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
    y <- initializeSparseMatrix(x)
   
    clusters <- sample(LETTERS, 100, replace=TRUE)
    agg <- aggregateAcrossCells.chan(y, list(cluster=clusters))

    for (u in unique(clusters)) {
        chosen <- clusters == u
        m <- match(u, agg$combinations[,1])

        submat <- x[,chosen,drop=FALSE]
        sum.expected <- Matrix::rowSums(submat)
        expect_equal(sum.expected, agg$sums[,m])

        detected.expected <- Matrix::rowSums(submat > 0)
        expect_equal(detected.expected, agg$detected[,m])
    }

    expect_identical(agg$counts, as.integer(table(clusters)))
})

test_that("aggregateAcrossCells works as expected for simple cases", {
    x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
    y <- initializeSparseMatrix(x)
   
    clusters <- sample(LETTERS, 100, replace=TRUE)
    samples <- sample(1:3, 100, replace=TRUE)
    agg <- aggregateAcrossCells.chan(y, list(sample=samples, cluster=clusters))

    full <- data.frame(samples, clusters)
    df <- unique(full)
    df <- df[order(df[,1], df[,2]),]

    expected.counts <- integer(nrow(df))
    for (u in seq_len(nrow(df))) {
        chosen <- clusters == df$clusters[u] & samples == df$samples[u]
        expected.counts[u] <- sum(chosen)
        m <- which(agg$combinations$cluster  == df$clusters[u] & agg$combinations$sample == df$samples[u])

        submat <- x[,chosen,drop=FALSE]
        sum.expected <- Matrix::rowSums(submat)
        expect_equal(sum.expected, agg$sums[,m])

        detected.expected <- Matrix::rowSums(submat > 0)
        expect_equal(detected.expected, agg$detected[,m])
    }

    expect_identical(agg$counts, expected.counts)
})



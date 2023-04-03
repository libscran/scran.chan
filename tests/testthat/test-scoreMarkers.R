# This contains some brief tests for scoreMarkers
# library(testthat); library(scran.chan); source("test-scoreMarkers.R")

library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
y <- initializeSparseMatrix(x)

test_that("scoreMarkers works as expected for simple cases", {
    g <- sample(4, ncol(x), replace=TRUE)
    markers <- scoreMarkers.chan(y, g)

    for (current in markers$statistics) {
        expect_false(is.unsorted(current$cohen.rank))
    }

    df1 <- markers$statistics[["1"]]
    expect_equal(df1$mean, Matrix::rowMeans(x[as.integer(rownames(df1)),g==1]))

    df4 <- markers$statistics[["4"]]
    expect_equal(df4$detected, Matrix::rowMeans(x[as.integer(rownames(df4)),g==4] > 0))
})

test_that("scoreMarkers works as expected for different sorting requests", {
    g <- sample(2, ncol(x), replace=TRUE)

    # Preserves the input order.
    markers <- scoreMarkers.chan(y, g, sort.by=NULL)
    for (current in markers$statistics) {
        expect_identical(as.integer(rownames(current)), seq_len(nrow(current)))
    }

    markers <- scoreMarkers.chan(y, g, sort.by="cohen.mean")
    for (current in markers$statistics) {
        expect_false(is.unsorted(rev(current$cohen.mean))) # reverse as it sorts in descending order.
    }
})

test_that("scoreMarkers works as expected for the other summaries", {
    g <- sample(4, ncol(x), replace=TRUE)
    markers <- scoreMarkers.chan(y, g)
    full <- scoreMarkers.chan(y, g, simple.means.only=FALSE)

    expect_identical(markers$statistics[["2"]]$logFC, full$statistics[["2"]]$logFC.mean) 
    expect_identical(markers$statistics[["3"]]$logFC, full$statistics[["3"]]$logFC.mean) 
})

test_that("scoreMarkers works as expected for the full pairwise statistics", {
    g <- sample(4, ncol(x), replace=TRUE)
    markers <- scoreMarkers.chan(y, g)
    full <- scoreMarkers.chan(y, g, all.pairwise=TRUE)
    expect_equal(markers$statistics, full$statistics)

    # Checking that the effects got ETL'd correctly.
    full <- scoreMarkers.chan(y, g, simple.means.only=FALSE, sort.by=NULL, all.pairwise=TRUE)
    for (effect in c("logFC", "delta.detected", "cohen", "auc")) {
        for (group in as.character(unique(g))) {
            all.effects <- lapply(full$pairwise[[group]], function(x) x[,effect])
            mean.effect <- rowMeans(do.call(cbind, all.effects))
            ref.effect <- full$statistics[[group]][,paste0(effect, ".mean")]
            expect_equal(mean.effect, ref.effect)
        }
    }
})


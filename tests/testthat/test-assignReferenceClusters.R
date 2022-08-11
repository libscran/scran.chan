# library(testthat); library(scran.chan); source("test-assignReferenceClusters.R")

x <- t(as.matrix(iris[,1:4]))
keep <- downsampleByNeighbors.chan(x, 10)
sub <- x[,keep]
                                                     
test_that("assignReferenceClusters works more or less as expected", {
    labs <- as.character(iris[keep,"Species"])
    full <- assignReferenceClusters.chan(sub, labs, x)
    expect_identical(length(full), ncol(x))

    # Compare to the reference; should be mostly diagonal.
    mat <- table(full, iris[,"Species"]) 
    expect_true(sum(diag(mat))/sum(mat) > 0.95)

    # Preserves factorhood.
    labs <- factor(iris[keep,"Species"])
    levels(labs) <- c(levels(labs), "WHEE")
    labs <- relevel(labs, "WHEE")
    full2 <- assignReferenceClusters.chan(sub, labs, x)
    expect_identical(full, as.character(full2))
    expect_identical(levels(full2)[1], "WHEE")
})

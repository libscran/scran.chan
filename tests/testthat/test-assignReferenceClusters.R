# library(testthat); library(scran.chan); source("test-assignReferenceClusters.R")

x <- t(as.matrix(iris[,1:4]))
keep <- downsampleByNeighbors.chan(x, 5)
sub <- x[,keep]
                                                     
test_that("assignReferenceClusters works more or less as expected", {
    labs <- as.character(iris[keep,"Species"])
    full <- assignReferenceClusters.chan(sub, labs, x)
    expect_identical(nrow(full), ncol(x))
    expect_true(all(full$best.prop > 0));
    expect_true(all(full$best.prop >= full$second.prop))

    # Compare to the reference; should be mostly diagonal.
    mat <- table(full$assigned, iris[,"Species"]) 
    expect_true(sum(diag(mat))/sum(mat) > 0.9)

    # Preserves factorhood.
    labs <- factor(iris[keep,"Species"])
    levels(labs) <- c(levels(labs), "WHEE")
    labs <- relevel(labs, "WHEE")
    full2 <- assignReferenceClusters.chan(sub, labs, x)
    expect_identical(full$assigned, as.character(full2$assigned))
    expect_identical(levels(full2$assigned)[1], "WHEE")
})

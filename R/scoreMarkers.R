#' Score marker genes
#'
#' Score marker genes using a variety of effect sizes from pairwise comparisons.
#'
#' @param x A list of matrix data like that produced by \code{\link{logNormCounts.chan}}.
#' @param groups A vector specifying the group assignment for each cell in \code{x}.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' Alternatively \code{NULL} if all cells belong to the same batch.
#' @param lfc Non-negative numeric scalar specifying the log-fold change threshold to use.
#'
#' @return A list containing \code{statistics}, a list of data frame of marker statistics. 
#' Each data frame corresponds to a group in \code{groups} and contains:
#' \itemize{
#' \item \code{mean}, the mean expression across all cells in the current group.
#' \item \code{detected}, proportion of cells with detectable expression in the current group.
#' \item \code{cohen.min}, the smallest Cohen's d across all pairwise comparisons involving the current group.
#' \item \code{cohen.mean}, the mean Cohen's d across all pairwise comparisons involving the current group.
#' \item \code{cohen.rank}, the minimum rank of the Cohen's d across all pairwise comparisons.
#' \item \code{auc.min}, the smallest AUC across all pairwise comparisons involving the current group.
#' \item \code{auc.mean}, the mean AUC across all pairwise comparisons involving the current group.
#' \item \code{auc.rank}, the minimum rank of the AUC across all pairwise comparisons.
#' }
#' Data frames are sorted by \code{cohen.rank} by default,
#'
#' If \code{batch} is supplied, this list will also contain \code{per.batch}.
#' This is a list containing \code{mean} and \code{detected}, each of which are lists of data frames containing the batch-specific statistics for each group.
#' (In this case, \code{statistics} contains the averaged statistics for each gene across batches.)
#' 
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' y <- initializeSparseMatrix(x)
#' normed <- logNormCounts.chan(y)
#'
#' # Making up some groups.
#' g <- sample(4, ncol(x), replace=TRUE)
#' markers <- scoreMarkers.chan(normed, g)
#' names(markers)
#' head(markers[[1]])
#'
#' @export
scoreMarkers.chan <- function(x, groups, batch=NULL, lfc=0) {
    groups <- transform_factor(groups, n = tatami_ncol(x))
    batch <- transform_factor(batch, n = tatami_ncol(x))
    output <- score_markers(x$pointer, groups$index, batch$index, lfc=lfc)

    formatted <- vector("list", length(groups$names))
    for (i in seq_along(formatted)) {
        df <- output$statistics[[i]]
        rownames(df) <- x$rownames
        formatted[[i]] <- df[order(df$cohen.rank),,drop=FALSE]
    }

    names(formatted) <- groups$names
    output$statistics <- formatted

    if (!is.null(batch$index)) {
        names(output$per.batch$mean) <- groups$names
        names(output$per.batch$detected) <- groups$names
        for (i in seq_along(groups$names)) {
            output$per.batch$mean[[i]] <- data.frame(output$per.batch$mean[[i]])
            dimnames(output$per.batch$mean[[i]]) <- list(x$rownames, batch$names)
            output$per.batch$detected[[i]] <- data.frame(output$per.batch$detected[[i]])
            dimnames(output$per.batch$detected[[i]]) <- list(x$rownames, batch$names)
        }
    }

    output
}

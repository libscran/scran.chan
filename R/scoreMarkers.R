#' Score marker genes
#'
#' Score marker genes using a variety of effect sizes from pairwise comparisons.
#'
#' @param x A list of matrix data like that produced by \code{\link{logNormCounts.chan}}.
#' @param groups A vector specifying the group assignment for each cell in \code{x}.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' Alternatively \code{NULL} if all cells belong to the same batch.
#' @param lfc Non-negative numeric scalar specifying the log-fold change threshold to use.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param simple.means.only Logical scalar indicating whether to only report the means for the simple effect sizes, i.e., log-fold change and delta-detected.
#' @param sort.by String specifying the column to use for sorting genes in descending order
#' (except if it ends with \code{.rank}, in which case it is sorted in ascending order).
#' If \code{NULL}, no sorting is performed.
#'
#' @return A list containing \code{statistics}, a list of data frame of marker statistics. 
#' Each data frame corresponds to a group in \code{groups} and contains:
#' \itemize{
#' \item \code{mean}, the mean expression across all cells in the current group.
#' \item \code{detected}, proportion of cells with detectable expression in the current group.
#' \item \code{logFC}, the mean of the log-fold changes in expression compared to other groups.
#' \item \code{delta.detected}, the mean of the difference in the detected proportions compared to other groups.
#' \item \code{cohen.min}, the smallest Cohen's d across all pairwise comparisons involving the current group.
#' \item \code{cohen.mean}, the mean Cohen's d across all pairwise comparisons involving the current group.
#' \item \code{cohen.rank}, the minimum rank of the Cohen's d across all pairwise comparisons.
#' \item \code{auc.min}, the smallest AUC across all pairwise comparisons involving the current group.
#' \item \code{auc.mean}, the mean AUC across all pairwise comparisons involving the current group.
#' \item \code{auc.rank}, the minimum rank of the AUC across all pairwise comparisons.
#' }
#' Rows are sorted by the specified column in \code{sort.by}.
#'
#' If \code{simple.means.only=FALSE}, \code{logFC} and \code{delta.detected} are renamed to \code{logFC.mean} and \code{delta.detected.mean}, respectively.
#' In addition, the corresponding \code{*.min} and \code{*.rank} columns are also reported.
#' These can be interpreted in a similar manner as \code{cohen.min}, \code{cohen.rank}, etc. for these effect sizes.
#'
#' If \code{batch} is supplied, this list will also contain \code{per.batch}.
#' This is a list containing \code{mean} and \code{detected}, each of which are lists of data frames containing the batch-specific statistics for each group.
#' (In this case, \code{statistics} contains the averaged statistics for each gene across batches.)
#'
#' @details
#' \code{min} is the most stringent summary statistic for identifying upregulated genes, 
#' as a gene must be strongly upregulated in every pairwise comparison to achieve a large \code{min} effect.
#' \code{rank} is the most generous as a gene may be highly ranked in any one pairwise comparison to achieve a high \code{rank},
#' and is useful for identifying the combination of best genes that distinguish a cluster from the others.
#' \code{mean} lies in between these two extremes.
#' 
#' Cohen's d and the AUC are generally give quite similar rankings,
#' but both are reported for some verisimilitude.
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
#' head(markers$statistics[[1]])
#'
#' @references
#' Lun A (2021).
#' scran::ScoreMarkers Class Reference.
#' \url{https://ltla.github.io/libscran/classscran_1_1ScoreMarkers.html}
#'
#' @export
scoreMarkers.chan <- function(x, groups, batch=NULL, lfc=0, num.threads=1, simple.means.only=TRUE, sort.by="cohen.rank") {
    groups <- transform_factor(groups, n = tatami_ncol(x))
    batch <- transform_factor(batch, n = tatami_ncol(x))
    output <- score_markers(x$pointer, groups$index, batch$index, lfc=lfc, nthreads=num.threads, simple_means_only=simple.means.only)

    formatted <- vector("list", length(groups$names))
    for (i in seq_along(formatted)) {
        df <- output$statistics[[i]]
        rownames(df) <- x$rownames

        if (!is.null(sort.by)) {
            df <- df[order(df[,sort.by], decreasing=grepl(".rank$", sort.by)),,drop=FALSE]
        }
        formatted[[i]] <- df
    }

    names(formatted) <- groups$names
    output$statistics <- formatted

    if (!is.null(batch$index)) {
        names(output$per.batch$mean) <- groups$names
        names(output$per.batch$detected) <- groups$names

        for (i in seq_along(groups$names)) {
            mean.df <- data.frame(output$per.batch$mean[[i]])
            rownames(mean.df) <- x$rownames # don't use dimnames<- as this can't handle NULL rownames.
            colnames(mean.df) <- batch$names
            output$per.batch$mean[[i]] <- mean.df

            detect.df <- data.frame(output$per.batch$detected[[i]])
            rownames(detect.df) <- x$rownames
            colnames(detect.df) <- batch$names
            output$per.batch$detected[[i]] <- detect.df
        }
    }

    output
}

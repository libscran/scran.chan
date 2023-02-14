#' Score feature set activity for each cell
#'
#' Compute per-cell scores for a given feature set, defined as the column sums of a rank-1 approximation to the submatrix for the feature set.
#' This uses the same approach implemented in the \pkg{GSDecon} package from Jason Hackney.
#'
#' @param x A list containing a log-expression matrix like that produced by \code{\link{logNormCounts.chan}}.
#' @param features Integer, logical or character vector specifying the features that belong to the set.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' Alternatively \code{NULL} if all cells belong to the same batch.
#' @param scale Logical scalar indicating whether to scale rows to unit variance within each batch before computing the weighting vectors.
#' @param num.threads Number of threads to use.
#'
#' @return List containing \code{scores}, a numeric vector of per-cell scores for each column in \code{x};
#' and \code{weights}, a numeric vector of per-feature weights for each feature in \code{features} (ordered by their row indices in \code{x}).
#'
#' @author Aaron Lun
#' @examples
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' y <- initializeSparseMatrix(x)
#' normed <- logNormCounts.chan(y)
#' 
#' scoreFeatureSet.chan(normed, c(1,3,5,10,20,100))
#' 
#' @export
scoreFeatureSet.chan <- function(x, features, batch=NULL, scale=FALSE, num.threads=1) {
    features <- to_logical(features, n=tatami_dim(x$pointer)[1], names=x$rownames)
    batch <- transform_factor(batch, n = tatami_ncol(x))
    score_feature_set(x$pointer, features, batch$index, scale, num.threads)
}

#' Log-normalize the counts 
#'
#' Compute log-transformed normalized expression values by applying scaling normalization to count data.
#'
#' @param x A list of matrix data like that produced by \code{\link{initializeSparseMatrix}}.
#' @param size.factors A numeric vector of length equal to the number of cells in \code{x},
#' containing positive size factors for all cells.
#'
#' @return A list like \code{x} where the count matrix is replaced with a log-expression matrix.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' y <- initializeSparseMatrix(x)
#' normed <- logNormCounts.chan(y)
#'
#' @export
logNormCounts.chan <- function(x, size.factors=NULL) {
    x$pointer <- log_norm_counts(x$pointer, size.factors)
    x
}

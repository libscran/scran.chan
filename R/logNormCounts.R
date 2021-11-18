#' Log-normalize the counts 
#'
#' Compute log-transformed normalized expression values by applying scaling normalization to count data.
#'
#' @param x A list of matrix data like that produced by \code{\link{initializeSparseMatrix}}.
#' @param size.factors A numeric vector of length equal to the number of cells in \code{x},
#' containing positive size factors for all cells.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' Alternatively \code{NULL} if all cells belong to the same batch.
#' @param batch.mode String indicating how \code{batch} should be handled.
#' If \code{"lowest"}, we downscale to the coverage of the lowest batch.
#' If \code{"perblock"}, we scale each batch to a mean of 1.
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
logNormCounts.chan <- function(x, size.factors=NULL, batch=NULL, batch.mode=c("lowest", "perbatch")) {
    batch <- transform_factor(batch)
    x$pointer <- log_norm_counts(x$pointer, size.factors, batch=batch$index, batch_mode=match.arg(batch.mode))
    x
}

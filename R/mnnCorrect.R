#' MNN correction
#'
#' Apply mutual nearest neighbor (MNN) correction to remove batch effects from a low-dimensional matrix.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells,
#' typically generated from \code{\link{runPCA.chan}}.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' @param k Integer scalar specifying the number of neighbors to use when identifying MNN pairs.
#' @param nmads Numeric scalar specifying the number of MADs to use when computing correction vectors.
#'
#' @return List containing:
#' \itemize{
#' \item \code{corrected}, a numeric matrix of the same dimensions as \code{x}, containing the corrected values.
#' \item \code{merge.order}, character vector containing the unique levels of \code{batch} in the automatically determined merge order.
#' The first batch is used as the reference and all other batches are iteratively merged.
#' \item \code{num.pairs}, integer vector of length equal to the number of batches minus 1.
#' This contains the number of MNN pairs at each merge.
#' }
#'
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' b <- sample(1:2, ncol(x), replace=TRUE)
#' corrected <- mnnCorrect.chan(x, b)
#' str(corrected)
#'
#' @export
mnnCorrect.chan <- function(x, batch, k=15, nmads=3) {
    batch <- transform_factor(batch, n = ncol(x))
    output <- mnn_correct(x, batch$index, k=k, nmads=nmads)
    output$merge.order <- batch$names[output$merge.order + 1L]
    output
}

#' MNN correction
#'
#' Apply mutual nearest neighbor (MNN) correction to remove batch effects from a low-dimensional matrix.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells,
#' typically generated from \code{\link{runPCA.chan}}.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' @param k Integer scalar specifying the number of neighbors to use when identifying MNN pairs.
#' @param nmads Numeric scalar specifying the number of MADs to use when computing correction vectors.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param order Vector containing levels of \code{batch} in the desired merge order.
#' If \code{NULL}, a suitable merge order is automatically determined.
#' @param reference.policy String specifying the policy to use to choose the first reference batch.
#' This can be based on the largest batch (\code{"max-size"}), the most variable batch (\code{"max-variance"}), some combination of those two (\code{"max-rss"}) or the first specified input (\code{"input"}).
#' Only used for automatic merges, i.e., when \code{order=NULL}. 
#'
#' @return List containing:
#' \itemize{
#' \item \code{corrected}, a numeric matrix of the same dimensions as \code{x}, containing the corrected values.
#' \item \code{merge.order}, character vector containing the unique levels of \code{batch} in the automatically determined merge order.
#' The first level in this vector is used as the reference batch; all other batches are iteratively merged to it.
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
mnnCorrect.chan <- function(x, batch, k=15, nmads=3, num.threads=1, order=NULL, reference.policy=c("max-size", "max-variance", "max-rss", "input")) {
    batch <- transform_factor(batch, n = ncol(x))

    if (!is.null(order)) {
        order <- match(order, batch$names)
        if (!identical(sort(order), seq_along(order))) {
            stop("'order' should contain unique values in 'batch'"); 
        }
        order <- order - 1L
    }

    output <- mnn_correct(x, batch$index, k=k, nmads=nmads, nthreads=num.threads, order=order, ref_policy = match.arg(reference.policy))
    output$merge.order <- batch$names[output$merge.order + 1L]
    output
}

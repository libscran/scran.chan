#' MNN correction
#'
#' Apply mutual nearest neighbor (MNN) correction to remove batch effects from a low-dimensional matrix.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells,
#' typically generated from \code{\link{runPCA.chan}}.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' @param k Integer scalar specifying the number of neighbors to use when identifying MNN pairs.
#' @param nmads Numeric scalar specifying the number of MADs to use when computing correction vectors.
#' @param mass.cap Integer scalar specifying the cap on the number of observations to use for center of mass calculations.
#' A value of 100,000 may be appropriate for speeding up correction of very large datasets.
#' If this is set to NULL, no cap is used.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param order Vector containing levels of \code{batch} in the desired merge order.
#' If \code{NULL}, a suitable merge order is automatically determined.
#' @param reference.policy String specifying the policy to use to choose the first reference batch.
#' This can be based on the largest batch (\code{"max-size"}), the most variable batch (\code{"max-variance"}), 
#' some combination of those two (\code{"max-rss"}) or the first specified input (\code{"input"}).
#' Only used for automatic merges, i.e., when \code{order=NULL}. 
#' @param approximate Logical scalar specifying whether to perform an approximate neighbor search.
#' Set to \code{FALSE} for historical reasons.
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
mnnCorrect.chan <- function(x, batch, k=15, nmads=3, mass.cap=NULL, order=NULL, reference.policy=c("max-rss", "max-size", "max-variance", "input"), approximate=FALSE, num.threads=1) {
    batch <- transform_factor(batch, n = ncol(x))

    if (!is.null(order)) {
        order <- match(order, batch$names)
        if (!identical(sort(order), seq_along(order))) {
            stop("'order' should contain unique values in 'batch'"); 
        }
        order <- order - 1L
    }

    if (is.null(mass.cap)) {
        mass.cap <- -1
    }

    output <- mnn_correct(
        x, 
        batch$index, 
        k=k, 
        nmads=nmads,
        mass_cap=mass.cap,
        nthreads=num.threads, 
        order=order, 
        ref_policy=match.arg(reference.policy), 
        approximate=approximate
    )

    output$merge.order <- batch$names[output$merge.order + 1L]
    output
}

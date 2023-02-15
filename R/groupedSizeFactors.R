#' Compute grouped size factors
#'
#' Correct for composition biases across groups of cells by computing size factors from per-group pseudo-bulk profiles.
#'
#' @param x A list of matrix data like that produced by \code{\link{initializeSparseMatrix}}.
#' @param group Vector or factor of length equal to the number of columns of \code{x},
#' containing the group assignment for each cell in \code{x}.
#' @param center Logical scalar specifying whether the size factors should be centered on output.
#' @param prior.count Numeric scalar specifying the prior count to add to the pseudo-bulk profiles.
#' Larger values improve stability at the cost of accurate removal of composition biases.
#' @param reference Identity of the group in \code{group} to be used as a reference.
#' If \code{NULL}, a suitable reference is automatically chosen.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return Numeric vector containing the size factor for each cell.
#' @author Aaron Lun
#' @examples
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' y <- initializeSparseMatrix(x)
#' 
#' # Obtain a sensible grouping with a short first-pass analysis.
#' normed <- logNormCounts.chan(y)
#' pcs <- runPCA.chan(norm)
#' clust <- clusterKmeans.chan(pcs$components, k=5)
#'
#' # Obtaining the grouped size factors.
#' groupedSizeFactors(y, clust$clusters)
#'
#' @export
groupedSizeFactors <- function(x, group, center=TRUE, prior.count=10, reference=NULL, num.threads=1) {
    g <- as.integer(factor(group)) - 1L
    grouped_size_factors(
        x, 
        g,
        center=center, 
        prior.count=prior.count, 
        reference=if (is.null(reference)) -1L else g[match(reference, group)], # index of the natural match
        nthreads=num.threads
    )
}

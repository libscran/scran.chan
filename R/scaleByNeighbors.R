#' Scale and combine multiple embeddings 
#'
#' Scale multiple embeddings (usually derived from different modalities across the same set of cells) so that their within-population variances are comparable,
#' and then combine them into a single embedding matrix for combined downstream analysis.
#'
#' @param x List of numeric matrices of PCs or other embeddings where rows are dimensions and columns are cells.
#' All entries should have the same number of columns.
#' @param k Integer scalar specifying the number of neighbors to use. 
#' @param approximate Logical scalar indicating whether an approximate neighbor search should be performed.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return List containing \code{scaling}, a vector of scaling factors to be aplied to each embedding;
#' and \code{combined}, a numeric matrix formed by scaling each entry of \code{x} and then \code{rbind}ing them together.
#'
#' @examples
#' pcs <- list(
#'     gene = matrix(rnorm(10000), ncol=200),
#'     protein = matrix(rnorm(1000, sd=3), ncol=200),
#'     guide = matrix(rnorm(2000, sd=5), ncol=200)
#' )
#'
#' out <- scaleByNeighbors.chan(pcs)
#' out$scaling
#' dim(out$combined)
#' @author Aaron Lun
#' @export
scaleByNeighbors.chan <- function(x, k = 20, approximate = TRUE, num.threads = 1) {
    ncols <- lapply(x, ncol)
    if (length(unique(ncols)) != 1L) {
        stop("all entries of 'x' should have the same number of columns")
    }

    scaling <- scale_by_neighbors(x, k, approximate, num.threads)
    names(scaling) <- names(x)
    for (i in seq_along(x)) {
        x[[i]] <- x[[i]] * scaling[i]
    }

    list(
         scaling = scaling,
         combined = do.call(rbind, x)
    )
}

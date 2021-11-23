#' Principal components analysis
#'
#' Perform an approximate principal components analysis using IRLBA.
#'
#' @param x A list containing a log-expression matrix like that produced by \code{\link{logNormCounts.chan}}.
#' @param num.comp Integer scalar specifying the number of top PCs to obtain.
#' @param subset Integer, logical or character vector specifying which features to use in the PCA (e.g., highly variable genes).
#' If \code{NULL}, all features in \code{x} are used.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' Alternatively \code{NULL} if all cells belong to the same batch.
#' @param batch.method String indicating how \code{batch} should be handled (if it is supplied).
#' \code{"block"} is equivalent to linear regression on \code{x} prior to PCA,
#' while \code{"weight"} will only weight each batch so that they contribute equally to the PCA.
#'
#' @return List containing \code{components}, containing the top principal components (note, columns correspond to cells);
#' and \code{prop.variance}, containing the proportion of variance explained by each component.
#'
#' @author Aaron Lun
#' @examples
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' y <- initializeSparseMatrix(x)
#' normed <- logNormCounts.chan(y)
#'
#' pcs <- runPCA.chan(normed)
#' dim(pcs$components)
#' barplot(pcs$prop.variance)
#'
#' # Only using features of interest.
#' subset <- 1:80
#' subbed <- runPCA.chan(normed, subset=subset)
#'
#' # Handling batch effects by blocking (i.e., regression)
#' # or by weighting to equalize contributions. 
#' b <- sample(1:3, ncol(x), replace=TRUE)
#' blocked <- runPCA.chan(normed, batch=b)
#' barplot(blocked$prop.variance)
#'
#' weighted <- runPCA.chan(normed, batch=b, batch.method="weight")
#' barplot(weighted$prop.variance)
#'
#' @export
runPCA.chan <- function(x, num.comp=50, subset=NULL, num.threads=1, batch=NULL, batch.method=c("block", "weight")) {
    if (!is.null(subset)) {
        subset <- to_logical(subset, n=tatami_dim(x$pointer)[1], names=x$rownames)
    }

    if (is.null(batch)) {
        run_pca(x$pointer, num.comp, subset, nthreads=num.threads)
    } else {
        batch.method <- match.arg(batch.method)
        batch <- transform_factor(batch, n = tatami_ncol(x))

        if (batch.method == "block") {
            run_blocked_pca(x$pointer, num.comp, batch$index, subset, nthreads=num.threads)
        } else {
            run_multibatch_pca(x$pointer, num.comp, batch$index, subset, nthreads=num.threads)
        }
    }
}

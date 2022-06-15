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
#' @param rotation Logical scalar indicating whether to report the rotation vectors.
#'
#' @return List containing:
#' \itemize{
#' \item \code{components}, a numeric matrix containing the top principal components.
#' Each row corresponds to a PC and each column corresponds to a cell.
#' \item \code{prop.variance}, containing the proportion of variance explained by each component.
#' \item \code{rotation}, a numeric matrix containing the rotation vectors.
#' Each column corresponds to a PC while each row corresponds to a feature.
#' If \code{subset} is provided, each row corresponds to a retained feature, ordered by its index in \code{x}.
#' }
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
runPCA.chan <- function(x, num.comp=50, subset=NULL, num.threads=1, batch=NULL, batch.method=c("block", "weight"), rotation=FALSE) {
    if (!is.null(subset)) {
        subset <- to_logical(subset, n=tatami_dim(x$pointer)[1], names=x$rownames)
    }

    if (is.null(batch)) {
        output <- run_pca(x$pointer, num.comp, subset, rotation=rotation, nthreads=num.threads)
    } else {
        batch.method <- match.arg(batch.method)
        batch <- transform_factor(batch, n = tatami_ncol(x))

        if (batch.method == "block") {
            output <- run_blocked_pca(x$pointer, num.comp, batch$index, subset, rotation=rotation, nthreads=num.threads)
        } else {
            output <- run_multibatch_pca(x$pointer, num.comp, batch$index, subset, rotation=rotation, nthreads=num.threads)
        }
    }

    if (rotation) {
        if (is.null(subset)) {
            rownames(output$rotation) <- x$rownames
        } else {
            rownames(output$rotation) <- x$rownames[subset]
        }
    }

    output
}

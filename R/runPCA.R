#' Principal components analysis
#'
#' Perform an approximate principal components analysis using IRLBA.
#'
#' @param x A list containing a log-expression matrix like that produced by \code{\link{logNormCounts.chan}}.
#' @param num.comp Integer scalar specifying the number of top PCs to obtain.
#' @param subset Logical vector specifying which features to use in the PCA (e.g., highly variable genes).
#' If \code{NULL}, all features in \code{x} are used.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return List containing \code{components}, containing the top principal components;
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
#' barplot(pcs$prop.variance)
#'
#' @export
runPCA.chan <- function(x, num.comp=50, subset=NULL, num.threads=1) {
    run_pca(x$pointer, num.comp, subset, nthreads=num.threads)
}

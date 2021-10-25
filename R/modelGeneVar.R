#' Model per-gene variance
#'
#' Model the per-gene variance in log-expression values, accounting for the mean-variance trend.
#'
#' @param x A list containing the log-expression matrix, typically produced by \code{\link{logNormCounts.chan}}.
#' @param span Numeric scalar containing the span to use for LOWESS smoothing of the trend.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return A data.frame of statistics for each gene in \code{x}.
#'
#' @author Aaron Lun
#' @examples
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' y <- initializeSparseMatrix(x)
#' normed <- logNormCounts.chan(y)
#'
#' df <- modelGeneVar.chan(normed)
#' df[1:10,]
#' plot(df$means, df$variances)
#' points(df$means, df$fitted, col="red", pch=16, cex=0.5)
#' 
#' @export
modelGeneVar.chan <- function(x, span = 0.3, num.threads = 1) {
    output <- model_gene_var(x$pointer, span, nthreads=num.threads)
    output <- do.call(data.frame, output)    
    rownames(output) <- x$rownames
    output
}

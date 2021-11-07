#' Model per-gene variance
#'
#' Model the per-gene variance in log-expression values, accounting for the mean-variance trend.
#'
#' @param x A list containing the log-expression matrix, typically produced by \code{\link{logNormCounts.chan}}.
#' @param span Numeric scalar containing the span to use for LOWESS smoothing of the trend.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' Alternatively \code{NULL} if all cells belong to the same batch.
#'
#' @return A list containing \code{statistics}, a data frame of statistics for each gene in \code{x}.
#'
#' If \code{batch} is supplied, this list will also contain \code{per.batch}, a list of data frames with statistics for each batch.
#' (In that case, \code{statistics} contains the averaged statistics for each gene across batches.)
#'
#' @author Aaron Lun
#' @examples
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' y <- initializeSparseMatrix(x)
#' normed <- logNormCounts.chan(y)
#'
#' dfl <- modelGeneVar.chan(normed)
#' df <- dfl$statistics
#' df[1:10,]
#' plot(df$means, df$variances)
#' points(df$means, df$fitted, col="red", pch=16, cex=0.5)
#' 
#' @export
modelGeneVar.chan <- function(x, span = 0.3, num.threads = 1, batch=NULL) {
    batch <- transform_factor(batch)
    output <- model_gene_var(x$pointer, span, nthreads=num.threads, batch=batch$index)

    rownames(output$statistics) <- x$rownames
    if (!is.null(batch$index)) {
        names(output$per.batch) <- batch$names
        for (i in seq_along(output$per.batch)) {
            rownames(output$per.batch[[i]]) <- x$rownames
        }
    }

    output
}

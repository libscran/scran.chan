#' Model per-gene variance
#'
#' Model the per-gene variance in log-expression values, accounting for the mean-variance trend.
#'
#' @param x A list containing the log-expression matrix, typically produced by \code{\link{logNormCounts.chan}}.
#' @param span Numeric scalar containing the span to use for LOWESS smoothing of the trend.
#' @param use.fixed.width Logical scalar specifying whether to apply a fixed-width constraint on the smoothing,
#' to avoid problems with large differences in density across the range of means.
#' @param fixed.width Numeric scalar specifying the window width when `use.fixed.width = TRUE`.
#' @param fixed.min.count Integer scalar specifying the minimum number of genes in a window when `use.fixed.width = TRUE`.
#' Windows with too few genes are expanded until they contain the required minimum. 
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
modelGeneVar.chan <- function(x, span = 0.3, use.fixed.width=FALSE, fixed.width=1, fixed.min.count=200, num.threads = 1, batch=NULL) {
    batch <- transform_factor(batch, n = tatami_ncol(x))
    output <- model_gene_var(
        x$pointer, 
        span=span, 
        use_fixed=use.fixed.width, 
        fixed_width=fixed.width, 
        min_count=fixed.min.count, 
        nthreads=num.threads, 
        batch=batch$index
    )

    rownames(output$statistics) <- x$rownames
    if (!is.null(batch$index)) {
        names(output$per.batch) <- batch$names
        for (i in seq_along(output$per.batch)) {
            rownames(output$per.batch[[i]]) <- x$rownames
        }
    }

    output
}

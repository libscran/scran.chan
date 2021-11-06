#' Compute per-cell QC metrics
#'
#' Calculate per-cell QC metrics from an initialized matrix.
#'
#' @param x List containing the output of \code{\link{initializeSparseMatrix}}.
#' @param subsets List of logical vectors specifying feature subsets of interest.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{sum}, a numeric vector containing the total sum for each cell.
#' \item \code{detected}, an integer vector containing the number of detected features per cell.
#' \item \code{subsets}, a list of numeric vectors containing the proportion of counts in each feature subset.
#' }
#' Each vector is of length equal to the number of cells.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#'
#' # Mocking up a subset:
#' sub <- list(Mito=rbinom(nrow(x), 1, 0.1) > 0)
#'
#' # Running the analysis:
#' y <- initializeSparseMatrix(x)
#' qc <- perCellQCMetrics.chan(y, sub)
#' summary(qc$sum)
#' summary(qc$detected)
#' summary(qc$subsets$Mito)
#'
#' @export
perCellQCMetrics.chan <- function(x, subsets, num.threads = 1) {
    nr <- tatami_dim(x$pointer)[1]

    # Converting to the expected logical vector.
    subsets <- as.list(subsets)
    subsets <- lapply(subsets, to_logical, n=nr, names=x$rownames)
 
    # Slapping on some names.
    output <- per_cell_qc_metrics(x$pointer, subsets, nthreads=num.threads)
    names(output$subsets) <- names(subsets)
    output
}

#' Compute per-cell ADT QC metrics
#'
#' Calculate per-cell QC metrics from an initialized matrix of ADT counts.
#'
#' @param x List containing the output of \code{\link{initializeSparseMatrix}},
#' corresponding to a matrix of ADT counts.
#' @param subsets List of logical vectors specifying tag subsets of interest, typically control tags like IgGs.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{sum}, a numeric vector containing the total ADT count for each cell.
#' \item \code{detected}, an integer vector containing the number of detected tags per cell.
#' \item \code{subsets}, a list of numeric vectors containing the total count of each control set. 
#' }
#' Each vector is of length equal to the number of cells.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#'
#' # Mocking up a control set.
#' sub <- list(IgG=rbinom(nrow(x), 1, 0.1) > 0)
#'
#' # Running the analysis:
#' y <- initializeSparseMatrix(x)
#' qc <- perCellAdtQcMetrics.chan(y, sub)
#' summary(qc$sum)
#' summary(qc$detected)
#' summary(qc$subsets$IgG)
#'
#' @export
perCellAdtQcMetrics.chan <- function(x, subsets, num.threads = 1) {
    nr <- tatami_dim(x$pointer)[1]

    # Converting to the expected logical vector.
    subsets <- as.list(subsets)
    subsets <- lapply(subsets, to_logical, n=nr, names=x$rownames)
 
    # Slapping on some names.
    output <- per_cell_adt_qc_metrics(x$pointer, subsets, nthreads=num.threads)
    names(output$subsets) <- names(subsets)
    output
}

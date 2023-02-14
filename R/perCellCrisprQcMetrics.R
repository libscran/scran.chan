#' Compute per-cell CRISPR QC metrics
#'
#' Calculate per-cell QC metrics from an initialized matrix of CRISPR guide counts.
#'
#' @param x List containing the output of \code{\link{initializeSparseMatrix}},
#' corresponding to a matrix of CRISPR guide counts.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{sum}, a numeric vector containing the total tag count for each cell.
#' \item \code{detected}, an integer vector containing the number of detected tags per cell.
#' \item \code{max_proportion}, a numeric vector specifying the proportion of counts assigned to the most abundant tag in each cell.
#' \item \code{max_index}, an integer vector containing the index of the most abundant tag in each cell.
#' }
#' Each vector is of length equal to the number of cells.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(20, 100, 0.1) * 100))
#'
#' y <- initializeSparseMatrix(x)
#' qc <- perCellCrisprQcMetrics.chan(y)
#' summary(qc$max_proportion)
#' table(qc$max_index)
#'
#' @export
perCellCrisprQcMetrics.chan <- function(x, num.threads = 1) {
    output <- per_cell_crispr_qc_metrics(x$pointer, nthreads=num.threads)
    output$max.index <- output$max.index + 1L
    output
}

#' Filter out low-quality cells
#'
#' Remove cells based on the filters defined from QC metrics.
#'
#' @param x A list containing the initialized matrix, as produced by \code{\link{initializeSparseMatrix}}.
#' @param discard Logical vector of length equal to the number of cells, specifying which cells in \code{x} should be removed.
#'
#' @return A copy of \code{x} with all the low-quality cells removed.
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
#' filters <- perCellQCFilters.chan(qc$sums, qc$detected, qc$subsets)
#' filtered <- filterCells.chan(y, filters$filters$overall)
#'
#' @export
filterCells.chan <- function(x, discard) {
    x$pointer <- filter_cells(x$pointer, discard)
    x$colnames <- x$colnames[!discard]
    x
}

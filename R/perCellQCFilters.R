#' Compute per-cell QC metrics
#'
#' Calculate per-cell QC metrics from an initialized matrix.
#'
#' @param sums Numeric vector containing the sum of counts for each cell.
#' @param detected Integer vector containing the total number of detected features for each cell.
#' @param subsets List of numeric vectors containing the proportion of counts assigned to each feature subset in each cell.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' Alternatively \code{NULL} if all cells belong to the same batch.
#' @param nmads Numeric scalar specifying the number of median absolute deviations to be used to detect outliers.
#'
#' @return A list containing:
#' \itemize{
#'     \item \code{filters}, a list containing:
#'     \itemize{
#'         \item \code{sums}, a logical vector indicating whether a cell was removed because its total count was too low.
#'         \item \code{detected}, a logical vector indicating whether a cell was removed because the number of detected features was too low.
#'         \item \code{subsets}, a list of logical vectors indicating whether a cell was removed because the proportion of counts in each feature subset was too high.
#'         \item \code{overall}, a logical vector indicating whether a cell was removed for any reason.
#'     }
#'     All logical vectors are of length equal to the number of cells.
#'     \item \code{thresholds}, a list containing:
#'     \itemize{
#'         \item \code{sums}, a numeric vector containing the minimum threshold on the total count (in each batch, if \code{batch} is not \code{NULL}).
#'         \item \code{detected}, a numeric vector containing the minimum threshold on the number of detected features.
#'         \item \code{subsets}, a list of numeric vectors containing the maximum threshold on each feature subset proportion.
#'     }
#'     Each numeric vector is of length equal to the number of blocks (default 1 if \code{block=NULL}).
#' }
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
#' str(filters)
#'
#' @export
perCellQCFilters.chan <- function(sums, detected, subsets, batch=NULL, nmads=3) {
    batch <- transform_factor(batch)
    filters <- per_cell_qc_filters(sums, detected, subsets, batch=batch$index, nmads=nmads)

    names(filters$filters$subsets) <- names(subsets)
    names(filters$thresholds$subsets) <- names(subsets)

    names(filters$thresholds$sum) <- batch$names
    names(filters$thresholds$detected) <- batch$names
    for (i in seq_along(subsets)) {
        names(filters$thresholds$subsets[[i]]) <- batch$names
    }

    filters
}

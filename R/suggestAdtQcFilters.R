#' Suggest filters for ADT QC metrics
#'
#' Suggest appropriate filters to be applied on the per-cell QC metrics for ADT data.
#'
#' @param detected Integer vector containing the total number of detected features for each cell.
#' @param subsets List of numeric vectors containing the total counts assigned to each feature subset in each cell.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' Alternatively \code{NULL} if all cells belong to the same batch.
#' @param min.detected.drop Numeric scalar specifying the minimum relative drop from the median number of detected features.
#' @param nmads Numeric scalar specifying the number of median absolute deviations to be used to detect outliers.
#'
#' @return A list containing:
#' \itemize{
#'     \item \code{filter}, a logical vector of length equal to the number of cells.
#'     True values indicate that a cell was considered to be low quality (for any reason) and removed.
#'     \item \code{thresholds}, a list containing:
#'     \itemize{
#'         \item \code{detected}, a numeric vector containing the minimum threshold on the number of detected features.
#'         \item \code{subsets}, a list of numeric vectors containing the maximum threshold on each feature subset total. 
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
#' sub <- list(IgG=rbinom(nrow(x), 1, 0.1) > 0)
#'
#' # Running the analysis:
#' y <- initializeSparseMatrix(x)
#' qc <- perCellAdtQcMetrics.chan(y, sub)
#' filters <- suggestAdtQcFilters.chan(qc$detected, qc$subsets)
#' str(filters)
#'
#' @export
suggestAdtQcFilters.chan <- function(detected, subsets, batch=NULL, min.detected.drop=0.1, nmads=3) {
    batch <- transform_factor(batch, n = length(detected))
    filters <- suggest_adt_qc_filters(detected, subsets, batch=batch$index, min_detected_drop=min.detected.drop, nmads=nmads)

    names(filters$thresholds$detected) <- batch$names
    names(filters$thresholds$subsets) <- names(subsets)
    for (i in seq_along(subsets)) {
        names(filters$thresholds$subsets[[i]]) <- batch$names
    }

    filters
}

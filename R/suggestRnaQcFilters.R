#' Suggest filters for RNA QC metrics
#'
#' Suggest appropriate filters to be applied on the per-cell QC metrics for RNA data.
#'
#' @param sums Numeric vector containing the sum of counts for each cell.
#' @param detected Integer vector containing the total number of detected features for each cell.
#' @param subsets List of numeric vectors containing the proportion of counts assigned to each feature subset in each cell.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' Alternatively \code{NULL} if all cells belong to the same batch.
#' @param nmads Numeric scalar specifying the number of median absolute deviations to be used to detect outliers.
#' @param ... Further arguments to pass to \code{suggestRnaQcFilters.chan}.
#'
#' @return A list containing:
#' \itemize{
#'     \item \code{filter}, a logical vector of length equal to the number of cells.
#'     True values indicate that a cell was considered to be low quality (for any reason) and removed.
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
#' qc <- perCellRnaQcMetrics.chan(y, sub)
#' filters <- suggestRnaQcFilters.chan(qc$sums, qc$detected, qc$subsets)
#' str(filters)
#'
#' @export
suggestRnaQcFilters.chan <- function(sums, detected, subsets, batch=NULL, nmads=3) {
    batch <- transform_factor(batch, n = length(sums))
    filters <- suggest_rna_qc_filters(sums, detected, subsets, batch=batch$index, nmads=nmads)

    names(filters$thresholds$sums) <- batch$names
    names(filters$thresholds$detected) <- batch$names

    names(filters$thresholds$subsets) <- names(subsets)
    for (i in seq_along(subsets)) {
        names(filters$thresholds$subsets[[i]]) <- batch$names
    }

    filters
}

#' @export
#' @rdname suggestRnaQcFilters.chan
perCellQCFilters.chan <- function(sums, detected, subsets, batch=NULL, ...) {
    res <- suggestRnaQcFilters.chan(sums, detected, subsets, batch=batch, ...)
    batch <- if (is.null(batch)) 1L else transform_factor(batch, n = length(sums))$index + 1L

    sum.thresh <- unname(res$thresholds$sums)[batch]
    by.sum <- sums < sum.thresh

    detected.thresh <- unname(res$thresholds$detected)[batch]
    by.detected <- detected < detected.thresh

    by.subset <- list()
    for (i in seq_along(subsets)) {
        sub.thresh <- unname(res$thresholds$subsets[[i]])[batch]
        by.subset[[i]] <- subsets[[i]] > sub.thresh
    }
    names(by.subset) <- names(subsets)

    res$filters <- list(
        sums = by.sum,
        detected = by.detected,
        subsets = by.subset,
        overall = res$filter
    )
    res$filter <- NULL

    res
}

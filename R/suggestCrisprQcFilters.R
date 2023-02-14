#' Suggest filters for CRISPR QC metrics
#'
#' Suggest appropriate filters to be applied on the per-cell QC metrics for CRISPR data.
#'
#' @param sums Numeric vector containing the sum of tag counts for each cell.
#' @param max.proportion Numeric vector containing the proportion of counts assigned to the most abundant tag within each cell.
#' @param batch Vector or factor of length equal to the number of cells, specifying the batch of origin for each cell.
#' Alternatively \code{NULL} if all cells belong to the same batch.
#' @param nmads Numeric scalar specifying the number of median absolute deviations to be used to detect outliers.
#'
#' @return A list containing:
#' \itemize{
#'     \item \code{filter}, a logical vector of length equal to the number of cells.
#'     True values indicate that a cell was considered to be low quality (for any reason) and removed.
#'     \item \code{thresholds}, a list containing:
#'     \itemize{
#'         \item \code{max.count}, a numeric vector containing the minimum threshold on the maximum count within each cell.
#'     }
#'     Each numeric vector is of length equal to the number of blocks (default 1 if \code{block=NULL}).
#' }
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(20, 100, 0.1) * 100))
#'
#' # Running the analysis:
#' y <- initializeSparseMatrix(x)
#' qc <- perCellCrisprQcMetrics.chan(y)
#' filters <- suggestCrisprQcFilters.chan(qc$sums, qc$max.proportion)
#' str(filters)
#'
#' @export
suggestCrisprQcFilters.chan <- function(sums, max.proportion, batch=NULL, nmads=3) {
    batch <- transform_factor(batch, n = length(sums))
    filters <- suggest_crispr_qc_filters(sums, max.proportion, batch=batch$index, nmads=nmads)
    names(filters$thresholds$max.count) <- batch$names
    filters
}

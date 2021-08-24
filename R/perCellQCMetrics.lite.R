#' Compute per-cell QC metrics
#'
#' Light-weight calculation of the per-cell QC metrics.
#'
#' @param x One of the supported matrix-like objects, containing genes as rows and cells as columns.
#' @param subsets List of logical vectors specifying feature subsets of interest.
#'
#' @return A list containing the total sum for each cell,
#' the number of detected features per cell,
#' and the proportion of counts in each feature subset.
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
#' qc <- diet.scran:::perCellQCMetrics.lite(x, sub)
#' summary(qc$sum)
#' summary(qc$detected)
#' summary(qc$subsets$Mito)
#'
#' @export
perCellQCMetrics.lite <- function(x, subsets=NULL) {
    output <- compute_qc_metrics(x, as.list(subsets))
    names(output) <- c("sum", "detected", "subsets")
    names(output$subsets) <- names(subsets)
    output
}

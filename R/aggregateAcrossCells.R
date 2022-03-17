#' Aggregate expression across cells
#'
#' Aggregate expression values across cells based on one or more grouping factors.
#' This is primarily used to create pseudo-bulk profiles for each cluster/sample combination.
#'
#' @param x A list of matrix data like that produced by \code{\link{initializeSparseMatrix}}.
#' Expression values are typically expected to be counts.
#' @param factors A list or data frame containing one or more grouping factors.
#' Each entry should be a factor of the same length as the number of cells in \code{x}.
#' @param num.threads Integer specifying the number of threads to be used for aggregation.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{sums}, a numeric matrix where each row corresponds to a gene and each column corresponds to a unique combination of grouping levels.
#' Each entry contains the summed expression across all cells with that combination. 
#' \item \code{detected}, an integer matrix where each row corresponds to a gene and each column corresponds to a unique combination of grouping levels.
#' Each entry contains the number of cells with detected expression in that combination.
#' \item \code{combinations}, a data frame describing the levels for each unique combination.
#' Rows of this data frame correspond to columns of \code{sums} and \code{detected},
#' while columns correspond to the factors in \code{factors}.
#' \item \code{counts}, the number of cells associated with each combination.
#' Each entry corresponds to a row of \code{combinations}.
#' \item \code{index}, an integer vector of length equal to the number of cells in \code{x}.
#' This specifies the combination in \code{combinations} to which each cell was assigned.
#' }
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' y <- initializeSparseMatrix(x)
#'
#' # Simple aggregation:
#' clusters <- sample(LETTERS, 100, replace=TRUE)
#' agg <- aggregateAcrossCells.chan(y, list(cluster=clusters))
#' str(agg)
#'
#' # Multi-factor aggregation
#' samples <- sample(1:5, 100, replace=TRUE)
#' agg2 <- aggregateAcrossCells.chan(y, list(cluster=clusters, sample=samples))
#' str(agg2)
#' 
#' @export
aggregateAcrossCells.chan <- function(x, factors, num.threads = 1) {
    f <- lapply(factors, factor) 
    f0 <- lapply(f, as.integer)

    output <- aggregate_across_cells(x$pointer, f0, num.threads)
    rownames(output$sums) <- rownames(output$detected) <- x$rownames 

    for (i in seq_along(output$combinations)) {
        current <- output$combinations[[i]]
        output$combinations[[i]] <- levels(f[[i]])[current]
    }

    names(output$combinations) <- names(factors)
    output$combinations <- data.frame(output$combinations)
    output$index <- output$index + 1L

    output
}

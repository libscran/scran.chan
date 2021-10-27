#' Score marker genes
#'
#' Score marker genes using a variety of effect sizes from pairwise comparisons.
#'
#' @param x A list of matrix data like that produced by \code{\link{logNormCounts.chan}}.
#' @param groups A vector specifying the group assignment for each cell in \code{x}.
#' 
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' y <- initializeSparseMatrix(x)
#' normed <- logNormCounts.chan(y)
#'
#' # Making up some groups.
#' g <- sample(4, ncol(x), replace=TRUE)
#' markers <- scoreMarkers.chan(normed, g)
#' names(markers)
#' head(markers[[1]])
#'
#' @export
scoreMarkers.chan <- function(x, groups) {
    f <- factor(groups)
    output <- score_markers(x$pointer, as.integer(f) - 1L)

    formatted <- vector("list", nlevels(f))
    for (i in seq_along(formatted)) {
        df <- data.frame(
            mean = output$means[[i]],
            detected = output$detected[[i]],
            cohen.min = output$cohen[[i]][[1]],
            cohen.mean = output$cohen[[i]][[2]],
            cohen.rank = output$cohen[[i]][[3]],
            auc.min = output$auc[[i]][[1]],
            auc.mean = output$auc[[i]][[2]],
            auc.rank = output$auc[[i]][[3]],
            row.names = x$rownames
        )

        formatted[[i]] <- df[order(df$cohen.rank),,drop=FALSE]
    }

    names(formatted) <- levels(f)
    formatted
}

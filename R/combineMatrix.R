#' Combine multiple initialized matrices
#'
#' @param ... Any number of lists of matrix data like that produced by \code{\link{initializeSparseMatrix}}.
#' All of them should have the same number of columns (for \code{rbindMatrix}) or columns (for \code{cbindMatrix}).
#'
#' @return A list containing the combined matrix data, along with their combined dimnames.
#'
#' @author Aaron Lun
#'
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(100, 10, 0.1) * 100))
#' rownames(x) <- sprintf("GENE_%i", seq_len(100))
#' colnames(x) <- LETTERS[1:10]
#'
#' y <- initializeSparseMatrix(x)
#' ry <- rbindMatrix(y, y)
#' cy <- cbindMatrix(y, y)
#'
#' @name combineMatrix
NULL

#' @export
#' @rdname combineMatrix
rbindMatrix <- function(...) {
    x <- list(...)
    list(
        pointer=combine_matrix(lapply(x, function(y) y$pointer), TRUE),
        rownames=.combine_names(x, 1, "rownames"),
        colnames=x[[1]]$colnames
    )
}

#' @export
#' @rdname combineMatrix
cbindMatrix <- function(...) {
    x <- list(...)
    list(
        pointer=combine_matrix(lapply(x, function(y) y$pointer), FALSE),
        rownames=x[[1]]$rownames,
        colnames=.combine_names(x, 2, "colnames")
    )
}

.combine_names <- function(x, dim, names) {
    names <- lapply(x, function(y) y[[names]])

    lost <- vapply(names, is.null, FUN.VALUE=TRUE)
    if (all(lost)) {
        return(NULL)
    }

    for (i in seq_along(names)) {
        if (lost[i]) {
            names[[i]] <- character(tatami_dim(x[[i]]$pointer)[dim])
        }
    }

    unlist(names)
}

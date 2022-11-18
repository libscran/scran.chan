prepare_indices <- function(indices, dim, names) {
    if (is.null(indices) || is.numeric(indices)) {
        indices
    } else if (is.logical(indices)) {
        if (length(indices) != dim) {
            stop("length of indices should be equal to the relevant dimension")
        }
        which(indices)
    } else if (is.character(indices)) {
        match(indices, names)
    } else {
        stop("unknown index type '", typeof (indices), "'") 
    }
}

#' Subset an initialized matrix
#'
#' @param x A list of matrix data like that produced by \code{\link{initializeSparseMatrix}}.
#' @param i,j Subsetting vectors along the rows (\code{i}) or columns (\code{j}).
#' Each vector may be numeric, in which case it should contain indices;
#' logical, specifying a filter of length equal to the number of rows/columns;
#' or character, containing the names of the rows/columns to retain.
#' This can also be set to \code{NULL} to not perform filtering on that dimension.
#'
#' @return A list containing the subsetted matrix data, along with correspondingly subsetted dimnames.
#'
#' @author Aaron Lun
#'
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(100, 10, 0.1) * 100))
#' rownames(x) <- sprintf("GENE_%i", seq_len(100))
#' colnames(x) <- LETTERS[1:10]
#' y <- initializeSparseMatrix(x)
#'
#' subsetMatrix(y, 1:5, c("A", "C", "E"))
#' subsetMatrix(y, 10:5, NULL)
#' @export
subsetMatrix <- function(x, i, j) {
    dims <- tatami_dim(x$pointer)
    i <- prepare_indices(i, dims[1], x$rownames)
    j <- prepare_indices(j, dims[2], x$colnames)

    x$pointer <- subset_matrix(x$pointer, i, j)
    if (!is.null(i)) {
        x$rownames <- x$rownames[i]
    }
    if (!is.null(j)) {
        x$colnames <- x$colnames[j]
    }

    x
}


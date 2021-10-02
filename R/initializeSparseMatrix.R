#' Initialize the count matrix
#'
#' Initialize a matrix object on the C++ side.
#' If the input is a \linkS4class{dgCMatrix}, initialization is performed without copying any data;
#' otherwise an efficient block realization is performed.
#'
#' @param x A matrix-like object.
#' This should contain non-negative integers.
#'
#' @return A list of internal components that can be used in other \pkg{scran.chan} functions.
#'
#' @details
#' Do not attempt to serialize this object; it contains a pointer to external memory, and will not be valid after a save/load cycle.
#' 
#' @examples
#' # Mocking up a count matrix:
#' x <- Matrix::rsparsematrix(1000, 100, 0.1)
#' y <- round(abs(x))
#'
#' stuff <- initializeSparseMatrix(y)
#' str(stuff)
#'
#' # Trying with a DelayedArray:
#' z <- DelayedArray::DelayedArray(y)
#' stuff2 <- initializeSparseMatrix(z)
#' str(stuff2)
#' 
#' @export
#' @importClassesFrom Matrix Matrix dgCMatrix
initializeSparseMatrix <- function(x, num.threads=1) {
    if (is(x, "dgCMatrix")) {
        ptr <- initialize_from_sparse(x@x, x@i, x@p, nrow(x), ncol(x))
        o <- seq_len(nrow(x)) - 1L
    } else {
        if (num.threads==1) {
            BPPARAM <- NULL
        } else {
            BPPARAM <- BiocParallel::MulticoreParam(num.threads)
        }

        m <- DelayedArray::blockApply(x, 
            grid=DelayedArray::colAutoGrid(x), 
            BPPARAM=BPPARAM, 
            as.sparse=NA, 
            FUN = function(block) {
                if (!is(block, "SparseArraySeed")) {
                    i <- which(block != 0, arr.ind=TRUE)
                    block <- DelayedArray::SparseArraySeed(dim(block), nzindex=i, nzdata=block[i], dimnames=dimnames(block))
                }
                block
            }
        )

        indices <- lapply(m, DelayedArray::nzindex)
        values <- lapply(m, DelayedArray::nzdata)
        ncols <- vapply(m, ncol, 0L)
        ptr <- initialize_from_blocks(indices, values, nrow(x), ncols, num.threads)
        o <- seq_len(nrow(x)) - 1L
    }

    list(pointer=ptr, order=o, rownames=rownames(x), colnames=colnames(x))
}

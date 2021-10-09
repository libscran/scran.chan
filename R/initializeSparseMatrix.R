#' Initialize the count matrix
#'
#' Initialize a matrix object on the C++ side.
#' If the input is a \linkS4class{dgCMatrix}, initialization is performed without copying any data;
#' otherwise an efficient block realization is performed.
#'
#' @param x A matrix-like object.
#' This should contain non-negative integers.
#' @param force.integer Logical scalar indicating whether double-precision \code{x} should be forced into integers.
#' @param no.sparse.copy Logical scalar indicating whether we should avoid a copy when \code{x} is a dgCMatrix.
#' This is more memory efficient if the data has already been loaded into memory.
#' If \code{TRUE}, any setting of \code{force.integer} is ignored.
#' @param num.threads Integer scalar specifying the number of threads to use when initializing \code{x}.
#' For dgCMatrix inputs, this is only relevant when \code{no.sparse.copy = FALSE}.
#'
#' @return A list of internal components that can be used in other \pkg{scran.chan} functions.
#'
#' @details
#' Do not attempt to serialize this object; it contains a pointer to external memory, and will not be valid after a save/load cycle.
#'
#' If \code{no.sparse.copy=TRUE}, the lifetime of the sparse matrix in \code{x} must be greater than the returned pointer. 
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
initializeSparseMatrix <- function(x, force.integer=TRUE, no.sparse.copy=TRUE, num.threads=1) {
    if (is(x, "dgCMatrix") && no.sparse.copy) {
        ptr <- initialize_from_dgCMatrix(x@x, x@i, x@p, nrow(x), ncol(x))

    } else {
        # We iterate across row-wise blocks and add them bit by bit to the
        # matrix so that we never have to load the entire matrix into memory.
        grid <- DelayedArray::rowAutoGrid(x)

        i <- 0L
        iterator <- function() {
            if (i >= length(grid)) {
                NULL
            } else {
                i <<- i + 1L
                grid[[i]]
            }
        }

        extractor <- function(vp) {
            # Doing the sort here, so we can exploit compute in multiple cores.
            output <- DelayedArray::read_block(x, vp, as.sparse=TRUE)
            nzi <- DelayedArray::nzindex(output)
            o <- order(nzi[,1], nzi[,2])
            list(row=nzi[o,1] - 1L, column=nzi[o,2] - 1L, value=DelayedArray::nzdata(output)[o], nrow=nrow(output))
        }

        NC <- ncol(x)
        is.int <- DelayedArray::type(x) == "integer" || force.integer
        ptr0 <- initialize_from_blocks(nrow(x), NC, is.int)
        reduction <- function(ignore, block) {
            add_new_block(ptr0, block$row, block$column, block$value, block$nrow, NC, is.int)
        }

        if (num.threads==1) {
            while (!is.null(g <- iterator())) {
                reduction(NULL, extractor(g))           
            }
        } else {
            BPPARAM <- BiocParallel::MulticoreParam(num.threads)
            BiocParallel::bpiterate(iterator, FUN=extractor, REDUCE=reduction, init=NULL, reduce.in.order=TRUE, BPPARAM=BPPARAM)
        }
        
        ptr <- finalize_all_blocks(ptr0, NC, is.int)
    }

    list(pointer=ptr, rownames=rownames(x), colnames=colnames(x))
}

#' Downsample cells based on their neighbors
#'
#' Implements a nearest neighbor-based downsampling scheme that preserves density while also guaranteeing representation of low-frequency subpopulations.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param k Integer scalar specifying the number of neighbors to use.
#' Larger values result in greater downsampling, roughly at a ratio of \code{k/2}.
#' @param approximate Logical scalar indicating whether an approximate neighbor search should be performed.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @return List containing:
#' \itemize{
#' \item \code{chosen}, an integer vector containing the column indices of \code{x} to retain.
#' \item \code{assigned}, integer vector of length equal to \code{ncol(x)},
#' containing the column index of the representative cell for each cell in \code{x}.
#' }
#'
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' down <- downsampleByNeighbors.chan(x, 10)
#' keep <- down$chosen
#' length(keep)
#'
#' # Visualizing what was kept.
#' plot(x[1,], x[2,], col=seq_len(ncol(x)) %in% keep + 1L)
#'
#' @export
downsampleByNeighbors.chan <- function(x, k, approximate=TRUE, num.threads=1) {
    out <- downsample_by_neighbors(x, k=k, approximate=approximate, nthreads=num.threads)
    out$chosen <- out$chosen + 1L
    out$assigned <- out$assigned + 1L
    out
}

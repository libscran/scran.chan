#' Run an approximate t-SNE
#'
#' Compute the t-stochastic neighbor embedding with assorted approximations for speed.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param perplexity Numeric scalar specifying the perplexity to use in the t-SNE algorithm.
#' @param interpolate Integer scalar specifying the number of grid points to use for interpolation.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @return A numeric matrix where rows are cells and columns are the two dimensions of the embedding.
#' 
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' embedding <- runTSNE.chan(x)
#' plot(embedding[,1], embedding[,2], col=iris[,5])
#' 
#' @export
runTSNE.chan <- function(x, perplexity=30, interpolate=200, num.threads=1) {
    neighbors <- build_nn_index(x)
    init <- initialize_tsne(neighbors, perplexity, interpolate, num.threads)
    output <- run_tsne(init, num.threads)
    t(output)
}
